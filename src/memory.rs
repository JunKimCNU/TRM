//! Optional process address-space ceiling, installed before scanning workers.
//!
//! RLIMIT_AS is a Linux virtual-address-space limit, not an RSS/cgroup memory
//! limit. It also covers stacks, shared mappings, and allocator reservations.
//! Reducing batches helps normal streaming inputs fit; it cannot bound a whole
//! FASTA record, candidate tables, or the allocator's own overhead. Allocator
//! exhaustion can terminate the process rather than return a recoverable error.

use serde::Serialize;

const MIB: u64 = 1024 * 1024;
const FIXED_RESERVE_BYTES: u64 = 16 * MIB;
const WORKER_RESERVE_BYTES: u64 = 4 * MIB;
const APPROXIMATE_BATCH_OVERHEAD: u64 = 8;
const MIN_PLANNED_BATCH_BASES: u64 = 4096;

#[derive(Clone, Debug, Serialize)]
pub struct MemoryPolicy {
    pub requested_limit_mib: Option<u64>,
    pub effective_virtual_limit_bytes: Option<u64>,
    pub startup_virtual_bytes: Option<u64>,
    pub requested_batch_bases: usize,
    pub effective_batch_bases: usize,
    pub enforcement: &'static str,
    pub batch_policy: &'static str,
}

/// Apply an explicitly requested ceiling before readers/workers are allocated.
///
/// A pre-existing lower soft/hard limit remains in force. Without an explicit
/// request this function does not inspect or change OS limits or batch sizing.
pub fn prepare(
    limit_mib: Option<u64>,
    threads: usize,
    requested_batch_bases: usize,
) -> Result<MemoryPolicy, String> {
    prepare_with_payload(limit_mib, threads, requested_batch_bases, 0)
}

/// Include explicitly retained per-base payload in the batch planning heuristic.
/// This does not change RLIMIT_AS semantics or guarantee that a scan fits.
pub fn prepare_with_payload(
    limit_mib: Option<u64>,
    threads: usize,
    requested_batch_bases: usize,
    extra_bytes_per_base: u64,
) -> Result<MemoryPolicy, String> {
    let Some(mib) = limit_mib else {
        return Ok(MemoryPolicy {
            requested_limit_mib: None,
            effective_virtual_limit_bytes: None,
            startup_virtual_bytes: None,
            requested_batch_bases,
            effective_batch_bases: requested_batch_bases,
            enforcement: "not_requested",
            batch_policy: "unchanged",
        });
    };
    let requested_bytes = mib
        .checked_mul(MIB)
        .filter(|&value| value > 0)
        .ok_or("--memory-limit-mib must be positive and fit in a byte count")?;
    install_limit(
        mib,
        requested_bytes,
        threads,
        requested_batch_bases,
        extra_bytes_per_base,
    )
}

/// This is deliberately a planning heuristic, not a memory-usage guarantee.
///
/// There are at most 2*t queued batches, t worker batches, and one producer
/// batch. A batch can overshoot by a complete record. Approximate overhead also
/// includes per-record strings and vectors; arbitrary IDs/short records can
/// exceed it. Leave additional space for workers, parser, tables, and output.
fn planned_batch_bases(
    ceiling_bytes: u64,
    startup_bytes: u64,
    threads: usize,
    requested_batch_bases: usize,
) -> Result<usize, String> {
    planned_batch_bases_with_payload(
        ceiling_bytes,
        startup_bytes,
        threads,
        requested_batch_bases,
        0,
    )
}

fn planned_batch_bases_with_payload(
    ceiling_bytes: u64,
    startup_bytes: u64,
    threads: usize,
    requested_batch_bases: usize,
    extra_bytes_per_base: u64,
) -> Result<usize, String> {
    if threads == 0 || requested_batch_bases == 0 {
        return Err("memory planning requires positive threads and batch bases".into());
    }
    let workers = u64::try_from(threads).map_err(|_| "thread count is too large")?;
    let batch_slots = workers
        .checked_mul(3)
        .and_then(|n| n.checked_add(1))
        .ok_or("thread count is too large for memory planning")?;
    let overhead = batch_slots
        .checked_mul(
            APPROXIMATE_BATCH_OVERHEAD
                .checked_add(extra_bytes_per_base)
                .ok_or("payload overhead exceeds memory planning range")?,
        )
        .ok_or("thread count is too large for memory planning")?;
    let reserve = workers
        .checked_mul(WORKER_RESERVE_BYTES)
        .and_then(|n| n.checked_add(FIXED_RESERVE_BYTES))
        .and_then(|n| n.checked_add(startup_bytes))
        .ok_or("memory planning byte count overflow")?;
    let floor = reserve
        .checked_add(
            MIN_PLANNED_BATCH_BASES
                .checked_mul(overhead)
                .ok_or("memory planning byte count overflow")?,
        )
        .ok_or("memory planning byte count overflow")?;
    if ceiling_bytes < floor {
        return Err(format!(
            "--memory-limit-mib is too small for this startup/batch plan: effective address-space \
             ceiling {ceiling_bytes} bytes, startup {startup_bytes} bytes, {threads} workers; \
             planning floor {} MiB (not a universal minimum or a success guarantee). \
             Increase the limit or reduce --threads",
            floor.div_ceil(MIB),
        ));
    }
    let budget = (ceiling_bytes - reserve) / overhead;
    let batch = budget.min(requested_batch_bases as u64);
    usize::try_from(batch).map_err(|_| "planned batch size exceeds this platform's usize".into())
}

#[cfg(target_os = "linux")]
fn parse_virtual_bytes(status: &str) -> Result<u64, String> {
    let line = status
        .lines()
        .find_map(|line| line.strip_prefix("VmSize:"))
        .ok_or("/proc/self/status has no VmSize; cannot plan an address-space limit")?;
    let mut fields = line.split_whitespace();
    let kib = fields
        .next()
        .ok_or("empty VmSize in /proc/self/status")?
        .parse::<u64>()
        .map_err(|_| "invalid VmSize in /proc/self/status")?;
    if fields.next() != Some("kB") || fields.next().is_some() {
        return Err("unexpected VmSize units in /proc/self/status".into());
    }
    kib.checked_mul(1024)
        .ok_or_else(|| "VmSize byte count overflow".into())
}

#[cfg(target_os = "linux")]
fn install_limit(
    requested_mib: u64,
    requested_bytes: u64,
    threads: usize,
    requested_batch_bases: usize,
    extra_bytes_per_base: u64,
) -> Result<MemoryPolicy, String> {
    let status = std::fs::read_to_string("/proc/self/status").map_err(|e| {
        format!("cannot read process address-space size for --memory-limit-mib: {e}")
    })?;
    let startup_bytes = parse_virtual_bytes(&status)?;
    let mut limits = libc::rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };
    // SAFETY: limits is valid writable memory for the OS rlimit ABI; RLIMIT_AS
    // is a supported Linux resource constant supplied by libc for this target.
    if unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut limits) } != 0 {
        return Err(format!(
            "cannot query RLIMIT_AS: {}",
            std::io::Error::last_os_error()
        ));
    }
    let requested = libc::rlim_t::try_from(requested_bytes)
        .map_err(|_| "--memory-limit-mib exceeds this platform's rlimit range")?;
    if requested == libc::RLIM_INFINITY {
        return Err("--memory-limit-mib cannot equal the unlimited-resource sentinel".into());
    }
    let effective = requested.min(limits.rlim_cur).min(limits.rlim_max);
    // rlim_t can have a narrower width on other Linux ABIs.
    #[allow(clippy::unnecessary_cast)]
    let effective_bytes = effective as u64;
    if effective_bytes <= startup_bytes {
        return Err(format!(
            "--memory-limit-mib effective address-space ceiling ({effective_bytes} bytes) is \
             at or below existing process VmSize ({startup_bytes} bytes); increase the limit"
        ));
    }
    let effective_batch_bases = if extra_bytes_per_base == 0 {
        planned_batch_bases(
            effective_bytes,
            startup_bytes,
            threads,
            requested_batch_bases,
        )?
    } else {
        planned_batch_bases_with_payload(
            effective_bytes,
            startup_bytes,
            threads,
            requested_batch_bases,
            extra_bytes_per_base,
        )?
    };
    limits.rlim_cur = effective;
    // SAFETY: limits is a fully initialized rlimit; only the soft address-space
    // limit is lowered. No global mutation of the caller's shell is performed.
    if unsafe { libc::setrlimit(libc::RLIMIT_AS, &limits) } != 0 {
        return Err(format!(
            "cannot apply --memory-limit-mib via RLIMIT_AS: {}",
            std::io::Error::last_os_error()
        ));
    }
    Ok(MemoryPolicy {
        requested_limit_mib: Some(requested_mib),
        effective_virtual_limit_bytes: Some(effective_bytes),
        startup_virtual_bytes: Some(startup_bytes),
        requested_batch_bases,
        effective_batch_bases,
        enforcement: "linux_rlimit_as_virtual_address_space",
        batch_policy: if effective_batch_bases < requested_batch_bases {
            "reduced_by_conservative_streaming_heuristic"
        } else {
            "requested_batch_fits_streaming_heuristic"
        },
    })
}

#[cfg(not(target_os = "linux"))]
fn install_limit(
    _requested_mib: u64,
    _requested_bytes: u64,
    _threads: usize,
    _requested_batch_bases: usize,
    _extra_bytes_per_base: u64,
) -> Result<MemoryPolicy, String> {
    Err(
        "--memory-limit-mib requires Linux RLIMIT_AS; no limit was applied. Use an OS/container \
         memory controller on this platform"
            .into(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_request_preserves_defaults() {
        let policy = prepare(None, 8, 4_000_000).unwrap();
        assert_eq!(policy.effective_batch_bases, 4_000_000);
        assert_eq!(policy.requested_limit_mib, None);
        assert_eq!(policy.effective_virtual_limit_bytes, None);
        assert_eq!(policy.enforcement, "not_requested");
    }

    #[test]
    fn invalid_byte_counts_fail_before_os_changes() {
        assert!(prepare(Some(0), 1, 4000).is_err());
        assert!(prepare(Some(u64::MAX), 1, 4000).is_err());
    }

    #[test]
    fn plans_only_reduce_batch_size() {
        assert_eq!(
            planned_batch_bases(1024 * MIB, 4 * MIB, 1, 1000).unwrap(),
            1000
        );
        let small = planned_batch_bases(64 * MIB, 4 * MIB, 4, 4_000_000).unwrap();
        let larger = planned_batch_bases(256 * MIB, 4 * MIB, 4, 4_000_000).unwrap();
        assert!(small >= MIN_PLANNED_BATCH_BASES as usize);
        assert!(small < larger && larger <= 4_000_000);
        assert!(planned_batch_bases(1024 * MIB, 4 * MIB, 8, 4_000_000).unwrap() >= larger);
    }

    #[test]
    fn impossible_plans_are_errors_not_wrapped_budgets() {
        assert!(planned_batch_bases(4 * MIB, 4 * MIB, 1, 1000).is_err());
        assert!(planned_batch_bases(24 * MIB, 4 * MIB, 8, 1000).is_err());
        assert!(planned_batch_bases(128 * MIB, 4 * MIB, 0, 1000).is_err());
        assert!(planned_batch_bases(128 * MIB, 4 * MIB, 1, 0).is_err());
        assert!(planned_batch_bases(u64::MAX, u64::MAX, 1, 1000).is_err());
    }

    #[test]
    fn more_workers_leave_less_batch_space() {
        let one = planned_batch_bases(256 * MIB, 4 * MIB, 1, usize::MAX).unwrap();
        let eight = planned_batch_bases(256 * MIB, 4 * MIB, 8, usize::MAX).unwrap();
        assert!(eight < one);
    }

    #[test]
    fn quality_payload_reduces_limited_batches_without_changing_default_plan() {
        let original = planned_batch_bases(128 * MIB, 4 * MIB, 4, 4_000_000).unwrap();
        assert_eq!(
            original,
            planned_batch_bases_with_payload(128 * MIB, 4 * MIB, 4, 4_000_000, 0).unwrap()
        );
        let quality =
            planned_batch_bases_with_payload(128 * MIB, 4 * MIB, 4, 4_000_000, 1).unwrap();
        assert!(quality < original);
        assert!(quality <= original * 8 / 9 + 1);
        assert!(
            planned_batch_bases_with_payload(128 * MIB, 4 * MIB, 4, 4_000_000, u64::MAX).is_err()
        );
        assert_eq!(
            prepare_with_payload(None, 4, 4_000_000, 1)
                .unwrap()
                .effective_batch_bases,
            4_000_000
        );
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn proc_vm_size_parser_requires_valid_bytes() {
        assert_eq!(
            parse_virtual_bytes("Name:\ttrm\nVmSize:\t 123 kB\n").unwrap(),
            125_952
        );
        assert!(parse_virtual_bytes("VmRSS:\t1 kB").is_err());
        assert!(parse_virtual_bytes("VmSize:\t1 MB").is_err());
        assert!(parse_virtual_bytes("VmSize:\tbogus kB").is_err());
        assert!(parse_virtual_bytes("VmSize:\t18446744073709551615 kB").is_err());
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn installed_limit_is_enforced_in_child() {
        const CHILD_FLAG: &str = "TRM_MEMORY_CAP_TEST_CHILD";
        if std::env::var_os(CHILD_FLAG).is_none() {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "memory::tests::installed_limit_is_enforced_in_child",
                    "--nocapture",
                ])
                .env(CHILD_FLAG, "1")
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "child limit test failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            return;
        }
        let vm =
            parse_virtual_bytes(&std::fs::read_to_string("/proc/self/status").unwrap()).unwrap();
        let mib = vm.div_ceil(MIB) + 96;
        let first = prepare(Some(mib), 1, 4_000_000).unwrap();
        assert_eq!(first.effective_virtual_limit_bytes, Some(mib * MIB));
        // A later request cannot raise a pre-existing stricter process limit.
        let second = prepare(Some(mib + 64), 1, 4_000_000).unwrap();
        assert_eq!(
            second.effective_virtual_limit_bytes,
            first.effective_virtual_limit_bytes
        );
        let len = usize::try_from(mib * MIB * 2).unwrap();
        // SAFETY: anonymous PROT_NONE mapping has no backing object, is never
        // dereferenced, and is released if it unexpectedly succeeds.
        let mapping = unsafe {
            libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_NONE,
                libc::MAP_PRIVATE | libc::MAP_ANONYMOUS,
                -1,
                0,
            )
        };
        let error = std::io::Error::last_os_error();
        if mapping != libc::MAP_FAILED {
            // SAFETY: mapping/len describe exactly the successful mmap above.
            unsafe {
                libc::munmap(mapping, len);
            }
            panic!("anonymous mapping larger than address-space ceiling unexpectedly succeeded");
        }
        assert_eq!(error.raw_os_error(), Some(libc::ENOMEM));
    }
}
