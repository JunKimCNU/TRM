"""Plan threads for the current native scan and Python batch runner.

This is a launch preflight, not an OS quota or a cross-runner scheduler.
Keep these execution-only options outside native arguments and cache keys.
"""
from __future__ import annotations

import json
import sys


def add_arguments(parser):
    parser.add_argument('--thread-budget', type=int,
                        help='reject before launch if planned total threads exceed this positive limit; separate from worker cpu-budget')
    parser.add_argument('--thread-reserve', type=int, default=0,
                        help='nonnegative threads reserved for other project jobs/helpers (default: 0; not automatically measured)')
    parser.add_argument('--resource-plan-only', action='store_true',
                        help='print scheduling resource plan as JSON and exit without scanner launch or output creation; not an input-integrity check')


def resource_plan(*, jobs, workers, selected_samples, timed, cpu_affinity,
                  thread_budget=None, thread_reserve=0):
    for name, value in (('jobs', jobs), ('workers', workers), ('cpu_affinity', cpu_affinity)):
        if value < 1:
            raise ValueError(f'{name} must be positive')
    if selected_samples < 0:
        raise ValueError('selected_samples must be nonnegative')
    if thread_budget is not None and thread_budget < 1:
        raise ValueError('--thread-budget must be positive')
    if thread_reserve < 0:
        raise ValueError('--thread-reserve must be nonnegative')
    active = min(jobs, selected_samples)
    python_threads = 1 + active
    native_threads = active * (workers + 1)
    time_threads = active * int(bool(timed))
    runner = python_threads + native_threads + time_threads
    planned = runner + thread_reserve
    limit = thread_budget if thread_budget is not None else cpu_affinity
    return {
        'schema_version': 1,
        'configured_jobs': jobs, 'workers_per_scan': workers,
        'selected_samples': selected_samples, 'effective_jobs': active,
        'python_threads': python_threads, 'native_threads': native_threads,
        'gnu_time_threads': time_threads, 'record_resources': bool(timed),
        'runner_threads': runner, 'thread_reserve': thread_reserve,
        'planned_total_threads': planned, 'thread_budget': thread_budget,
        'cpu_affinity': cpu_affinity, 'comparison_limit': limit,
        'budget_mode': 'reject' if thread_budget is not None else 'warn',
        'within_budget': planned <= limit,
        'exceeds_cpu_affinity': planned > cpu_affinity,
        'formula': '1 + min(jobs, selected_samples) * (workers + 2 + GNU_time) + thread_reserve',
        'scope': 'Current native worker/main topology plus Python main/pool and GNU time; selected cached samples still count in case rescan is needed.',
        'limitations': 'Planning estimate, not measured use or an OS limit. Other runners, custom binary helpers and observers must fit the explicit reserve. No global resource arbitration or memory enforcement.',
    }


def preflight(parser, args, *, selected_samples, timed, cpu_affinity):
    """Attach the plan; reject explicit overflow; return True for plan-only."""
    try:
        plan = resource_plan(jobs=args.jobs, workers=args.threads,
                             selected_samples=selected_samples, timed=timed,
                             cpu_affinity=cpu_affinity, thread_budget=args.thread_budget,
                             thread_reserve=args.thread_reserve)
    except ValueError as error:
        parser.error(str(error))
    args._resource_plan = plan
    if args.resource_plan_only:
        print(json.dumps(plan, indent=2), flush=True)
    total = plan['planned_total_threads']
    if plan['budget_mode'] == 'reject' and not plan['within_budget']:
        parser.error(f'planned total threads ({total}) exceeds --thread-budget ({args.thread_budget}); '
                     'includes producer, Python coordinator/pool, GNU time and --thread-reserve. '
                     'Lower --jobs/--threads or allocate a sufficient explicit budget; no automatic resizing.')
    if plan['exceeds_cpu_affinity']:
        print(f'WARNING: planned total threads ({total}) exceeds CPU affinity count ({cpu_affinity}); '
              '--cpu-budget counts workers only. Use --thread-budget with an explicit --thread-reserve '
              'to check the full planned allocation before launch.', file=sys.stderr)
    elif not args.resource_plan_only:
        print(f'Resource plan: {total} threads including reserve {args.thread_reserve}; '
              f"{plan['budget_mode']} limit {plan['comparison_limit']} (not an OS quota)", file=sys.stderr)
    return args.resource_plan_only
