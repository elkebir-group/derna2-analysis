"""icount.py OUT.json -- CMD...: run CMD and count its user-space instructions
and cycles over every thread and child (perf_event_open, inherit), plus wall.

Instructions retired do not move with host load the way wall and CPU time do,
so an A/B on a busy host can compare work. Cycles do move (cache and memory
contention), which is why both are kept.
"""
import ctypes
import json
import os
import struct
import sys
import time

PERF_TYPE_HARDWARE = 0
COUNTERS = {"instructions": 1, "cycles": 0}
NR_PERF_EVENT_OPEN = 298  # x86_64
# flags bitfield: disabled(0) inherit(1) exclude_kernel(5) exclude_hv(6) enable_on_exec(12)
FLAGS = (1 << 0) | (1 << 1) | (1 << 5) | (1 << 6) | (1 << 12)

libc = ctypes.CDLL(None, use_errno=True)


def open_counter(config):
    # perf_event_attr, PERF_ATTR_SIZE_VER0 (64 bytes)
    attr = struct.pack("=IIQQQQQIIQ", PERF_TYPE_HARDWARE, 64, config, 0, 0, 0, FLAGS, 0, 0, 0)
    buf = ctypes.create_string_buffer(attr, 64)
    fd = libc.syscall(NR_PERF_EVENT_OPEN, buf, 0, -1, -1, 0)
    if fd < 0:
        raise OSError(ctypes.get_errno(), "perf_event_open")
    return fd


def main():
    out, sep, cmd = sys.argv[1], sys.argv[2], sys.argv[3:]
    assert sep == "--" and cmd
    fds = {name: open_counter(cfg) for name, cfg in COUNTERS.items()}
    t0 = time.time()
    pid = os.fork()
    if pid == 0:
        os.execvp(cmd[0], cmd)
    _, status, usage = os.wait4(pid, 0)
    wall = time.time() - t0
    rec = {name: struct.unpack("Q", os.read(fd, 8))[0] for name, fd in fds.items()}
    rec.update(wall_s=wall, user_s=usage.ru_utime, sys_s=usage.ru_stime,
               maxrss_kb=usage.ru_maxrss, status=os.waitstatus_to_exitcode(status), cmd=cmd)
    with open(out, "w") as f:
        json.dump(rec, f, indent=1)
    print(f"instructions {rec['instructions'] / 1e9:.1f} G  cycles {rec['cycles'] / 1e9:.1f} G  "
          f"wall {wall:.1f} s  user {usage.ru_utime:.1f} s", file=sys.stderr)
    sys.exit(rec["status"])


main()
