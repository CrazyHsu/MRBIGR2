"""Worker entry for ``mrbigr.mcp.longjob``.

Spawned by ``start_or_poll`` as::

    python -m mrbigr.mcp.longjob_worker <job_file>

Reads the job JSON, dispatches according to ``runner.type``:

* ``shell``  — execve a subprocess with ``runner.cmd`` (list[str]) in
  ``runner.cwd`` (optional). Captures returncode.
* ``python`` — ``importlib.import_module(runner.module)`` then call
  ``getattr(module, runner.function)(**runner.kwargs)``. Treats a non-None
  return as success.

On completion, updates the job file with ``status``, ``completed_at``,
``returncode``, and (for python runners) ``result``. Exits non-zero on
failure so the parent's ``proc.poll()`` reflects it.
"""
from __future__ import annotations

import importlib
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _run_shell(runner: dict) -> int:
    cmd = runner.get("cmd")
    if not cmd:
        raise ValueError("shell runner requires 'cmd'")
    cwd = runner.get("cwd")
    if isinstance(cmd, str):
        # shell=True path
        proc = subprocess.run(cmd, shell=True, cwd=cwd)
    else:
        proc = subprocess.run(list(cmd), cwd=cwd)
    return int(proc.returncode)


def _run_python(runner: dict) -> tuple[int, object]:
    module_name = runner.get("module")
    function_name = runner.get("function")
    kwargs = dict(runner.get("kwargs") or {})
    if not module_name or not function_name:
        raise ValueError("python runner requires 'module' and 'function'")
    module = importlib.import_module(module_name)
    fn = getattr(module, function_name)
    result = fn(**kwargs)
    # Explicit failure signals only. Equality checks against DataFrame/array-like
    # results can be ambiguous, so avoid tuple membership here.
    if result is None or result is False:
        return 1, result
    return 0, result


def main(job_file: str) -> int:
    job_path = Path(job_file).resolve()
    job = _read(job_path)
    runner = job.get("runner") or {}
    rtype = runner.get("type")

    job["worker_pid"] = os.getpid()
    job["worker_started_at"] = time.time()
    _write(job_path, job)

    try:
        if rtype == "shell":
            rc = _run_shell(runner)
            job["returncode"] = rc
            job["status"] = "completed" if rc == 0 else "failed"
        elif rtype == "python":
            rc, result = _run_python(runner)
            job["returncode"] = rc
            job["status"] = "completed" if rc == 0 else "failed"
            if rc == 0:
                # Persist serializable result (e.g. list of output files).
                try:
                    json.dumps(result)
                    job["result"] = result
                except (TypeError, ValueError):
                    job["result"] = str(result)
        else:
            job["status"] = "failed"
            job["returncode"] = 2
            job["error"] = f"unknown runner.type: {rtype!r}"
    except Exception as exc:  # noqa: BLE001 - persist worker failure for polling.
        job["status"] = "failed"
        job["returncode"] = 1
        job["error"] = repr(exc)
        job["traceback"] = traceback.format_exc()
    finally:
        job["completed_at"] = time.time()
        _write(job_path, job)

    return int(job.get("returncode") or 0)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python -m mrbigr.mcp.longjob_worker <job_file>", file=sys.stderr)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
