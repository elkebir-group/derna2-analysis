#!/usr/bin/env python3
"""measure.py OUT.json [--timeout S] [--stdin FILE] [--stdout FILE] -- CMD...

Run CMD once and record what it cost: user-space instructions and cycles over
every thread and child (perf_event_open, inherit), wall-clock, user and system
CPU time, peak RSS, and how it ended.

The counters are opened by this process, so each measured command needs its
own wrapper: a counter opened with `inherit` counts every child forked after
it, and a driver that opened counters for several concurrent cells would charge
each cell for all of them.

Instructions retired do not move with host load the way wall and CPU time do.
Cycles do (cache and memory contention), which is why both are kept.

Exit categories are kept distinct, because a censored cell is a lower bound and
not a missing row:
  ok       exit status 0
  timeout  killed by this wrapper at --timeout
  oom      killed by SIGKILL that this wrapper did not send (the kernel OOM
           killer), or the command reported an allocation failure
  crash    any other signal or non-zero exit
"""

import argparse
import ctypes
import json
import os
import re
import signal
import struct
import sys
import time

PERF_TYPE_HARDWARE = 0
COUNTERS = {"instructions": 1, "cycles": 0}
NR_PERF_EVENT_OPEN = 298  # x86_64
# flags bitfield: disabled(0) inherit(1) exclude_kernel(5) exclude_hv(6) enable_on_exec(12)
FLAGS = (1 << 0) | (1 << 1) | (1 << 5) | (1 << 6) | (1 << 12)
ALLOC_FAILURE = re.compile(r"memory allocation of \d+ bytes failed|std::bad_alloc|Cannot allocate memory")

libc = ctypes.CDLL(None, use_errno=True)


def open_counter(config):
    """One hardware counter on this process and its future children, or None
    when the kernel refuses (the run is still timed, just not counted)."""
    attr = struct.pack("=IIQQQQQIIQ", PERF_TYPE_HARDWARE, 64, config, 0, 0, 0, FLAGS, 0, 0, 0)
    buf = ctypes.create_string_buffer(attr, 64)
    fd = libc.syscall(NR_PERF_EVENT_OPEN, buf, 0, -1, -1, 0)
    return fd if fd >= 0 else None


def read_counter(fd):
    return struct.unpack("Q", os.read(fd, 8))[0] if fd is not None else None


def spawn(cmd, stdin, stdout, stderr):
    """Fork and exec CMD in its own process group, so a timeout kills every
    process it started and not just the leader."""
    pid = os.fork()
    if pid == 0:
        os.setpgid(0, 0)
        for path, target, flags in ((stdin, 0, os.O_RDONLY), (stdout, 1, os.O_WRONLY | os.O_CREAT | os.O_TRUNC),
                                    (stderr, 2, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)):
            if path:
                fd = os.open(path, flags, 0o644)
                os.dup2(fd, target)
                os.close(fd)
        try:
            os.execvp(cmd[0], cmd)
        except OSError as err:
            os.write(2, f"exec failed: {err}\n".encode())
            os._exit(127)
    return pid


def wait_with_timeout(pid, timeout):
    """wait4 the child, killing its process group at the deadline.
    Returns (status, rusage, timed_out)."""
    deadline = time.monotonic() + timeout if timeout else None
    while True:
        done, status, usage = os.wait4(pid, os.WNOHANG)
        if done:
            return status, usage, False
        if deadline and time.monotonic() >= deadline:
            os.killpg(pid, signal.SIGKILL)
            _, status, usage = os.wait4(pid, 0)
            return status, usage, True
        time.sleep(0.05)


def classify(status, timed_out, stderr_path):
    if timed_out:
        return "timeout"
    if os.WIFSIGNALED(status):
        return "oom" if os.WTERMSIG(status) == signal.SIGKILL else "crash"
    if os.WEXITSTATUS(status) == 0:
        return "ok"
    if stderr_path and os.path.exists(stderr_path):
        with open(stderr_path, errors="replace") as f:
            if ALLOC_FAILURE.search(f.read()):
                return "oom"
    return "crash"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--timeout", type=float, default=0)
    ap.add_argument("--stdin")
    ap.add_argument("--stdout")
    ap.add_argument("--stderr")
    argv = sys.argv[1:]
    split = argv.index("--") if "--" in argv else len(argv)
    a = ap.parse_args(argv[:split])
    cmd = argv[split + 1:]
    assert cmd, "no command after --"
    fds = {name: open_counter(cfg) for name, cfg in COUNTERS.items()}
    t0 = time.time()
    pid = spawn(cmd, a.stdin, a.stdout, a.stderr)
    status, usage, timed_out = wait_with_timeout(pid, a.timeout)
    wall = time.time() - t0
    rec = {name: read_counter(fd) for name, fd in fds.items()}
    rec.update(wall_s=wall, user_s=usage.ru_utime, sys_s=usage.ru_stime, maxrss_kb=usage.ru_maxrss,
               exit=classify(status, timed_out, a.stderr), raw_status=status,
               cpus=sorted(os.sched_getaffinity(0)), cmd=cmd)
    with open(a.out, "w") as f:
        json.dump(rec, f, indent=1)
    sys.exit(0 if rec["exit"] == "ok" else 1)


main()
