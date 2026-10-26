//! Opt-in producer wall-clock attribution at buffer and batch boundaries.
//! Times include scheduling waits; these are not CPU counters or disk bytes.
use serde::Serialize;
use std::io::{self, Read};
use std::sync::{Arc, Mutex};
use std::time::Instant;

#[derive(Default, Debug, Serialize)]
pub struct InputMetrics {
    pub path: String,
    pub gzip: bool,
    pub source_read_calls: u64,
    pub source_bytes: u64,
    pub source_read_wall_seconds: f64,
    pub decoded_stream_read_calls: u64,
    pub decoded_stream_bytes: u64,
    pub decoded_stream_read_wall_seconds: f64,
    pub source_read_nested_in_decode_wall_seconds: f64,
    pub decode_and_buffer_excluding_source_read_wall_seconds: f64,
}

pub type SharedInput = Arc<Mutex<InputMetrics>>;

pub fn input_metrics(path: &str) -> SharedInput {
    Arc::new(Mutex::new(InputMetrics {
        path: path.into(),
        ..Default::default()
    }))
}

struct TimedRead {
    source: Box<dyn Read + Send>,
    metrics: SharedInput,
    decoded: bool,
}

impl Read for TimedRead {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let before_source = if self.decoded {
            self.metrics
                .lock()
                .map_err(|_| io::Error::other("profile lock poisoned"))?
                .source_read_wall_seconds
        } else {
            0.0
        };
        let start = Instant::now();
        let result = self.source.read(buffer);
        let elapsed = start.elapsed().as_secs_f64();
        let mut metrics = self
            .metrics
            .lock()
            .map_err(|_| io::Error::other("profile lock poisoned"))?;
        let bytes = result.as_ref().copied().unwrap_or(0) as u64;
        if self.decoded {
            let nested_source = (metrics.source_read_wall_seconds - before_source).max(0.0);
            metrics.decoded_stream_read_calls += 1;
            metrics.decoded_stream_bytes += bytes;
            metrics.decoded_stream_read_wall_seconds += elapsed;
            metrics.source_read_nested_in_decode_wall_seconds += nested_source;
            metrics.decode_and_buffer_excluding_source_read_wall_seconds +=
                (elapsed - nested_source).max(0.0);
        } else {
            metrics.source_read_calls += 1;
            metrics.source_bytes += bytes;
            metrics.source_read_wall_seconds += elapsed;
        }
        result
    }
}

pub fn timed_source(source: Box<dyn Read + Send>, metrics: SharedInput) -> Box<dyn Read + Send> {
    Box::new(TimedRead {
        source,
        metrics,
        decoded: false,
    })
}

pub fn timed_decoded(source: Box<dyn Read + Send>, metrics: SharedInput) -> Box<dyn Read + Send> {
    Box::new(TimedRead {
        source,
        metrics,
        decoded: true,
    })
}

pub struct Producer {
    started: Instant,
    elapsed_seconds: f64,
    pub queue_send_wall_seconds: f64,
    pub batches_sent: u64,
    pub inputs: Vec<SharedInput>,
}

impl Producer {
    pub fn new() -> Self {
        Self {
            started: Instant::now(),
            elapsed_seconds: 0.0,
            queue_send_wall_seconds: 0.0,
            batches_sent: 0,
            inputs: Vec::new(),
        }
    }

    pub fn finish(&mut self) {
        self.elapsed_seconds = self.started.elapsed().as_secs_f64();
    }

    pub fn json(&self, records: u64, bases: u64) -> Result<serde_json::Value, String> {
        let mut inputs = Vec::new();
        let (mut io, mut decode) = (0.0, 0.0);
        for shared in &self.inputs {
            let metrics = shared.lock().map_err(|_| "profile lock poisoned")?;
            io += metrics.source_read_wall_seconds;
            decode += metrics.decode_and_buffer_excluding_source_read_wall_seconds;
            inputs.push(serde_json::to_value(&*metrics).map_err(|e| e.to_string())?);
        }
        Ok(serde_json::json!({
            "schema_version": 1,
            "producer_wall_seconds": self.elapsed_seconds,
            "source_read_wall_seconds": io,
            "decode_and_buffer_excluding_source_read_wall_seconds": decode,
            "queue_send_wall_seconds": self.queue_send_wall_seconds,
            "parse_mate_batch_and_control_residual_wall_seconds": (self.elapsed_seconds - io - decode - self.queue_send_wall_seconds).max(0.0),
            "batches_sent": self.batches_sent,
            "records": records,
            "bases": bases,
            "inputs": inputs,
            "method": "Opt-in buffer/batch wall clocks, not CPU time. Source reads include gzip magic reads, OS/page-cache/filesystem waits and scheduler delays; bytes are bytes returned to this reader, not physical disk IO. Decoded-stream calls include nested source reads; subtract their separately measured nested time for decode/buffering attribution. Residual groups parsing, mate checks, allocation, batch construction, input opening and instrumentation/control overhead. Producer overlaps workers; no worker time is measured. Prefix streams may read ahead. No per-record timer; no claim that IO/decode components are isolated causal costs."
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn nested_read_bytes_and_time_are_not_added_twice() {
        let metrics = input_metrics("memory");
        let inner = timed_source(Box::new(Cursor::new(vec![7u8; 101])), metrics.clone());
        let mut outer = timed_decoded(inner, metrics.clone());
        let mut data = Vec::new();
        outer.read_to_end(&mut data).unwrap();
        let m = metrics.lock().unwrap();
        assert_eq!(data.len(), 101);
        assert_eq!(m.source_bytes, 101);
        assert_eq!(m.decoded_stream_bytes, 101);
        assert_eq!(m.source_read_calls, m.decoded_stream_read_calls);
        assert!(m.decoded_stream_read_wall_seconds >= m.source_read_nested_in_decode_wall_seconds);
        assert_eq!(
            m.source_read_wall_seconds,
            m.source_read_nested_in_decode_wall_seconds
        );
    }
}
