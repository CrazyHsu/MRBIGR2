"""Unit tests for ``mrbigr.mcp.longjob``.

Avoids depending on real GEMMA / PLINK binaries — uses small shell sleep +
file-write fakes to validate the framework's spawn/poll/wait/dedup contract
end to end.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import pytest

# Make ``mrbigr`` importable when the test runs without ``pip install -e .``.
REPO_SRC = Path(__file__).resolve().parent.parent / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from mrbigr.mcp import longjob  # noqa: E402
from mrbigr.mcp.longjob import (  # noqa: E402
    JobSpec,
    WAIT_FOR_JOB_MAX_SECONDS,
    poll_job,
    register_wait_for_job,
    start_or_poll,
    wait_for_job,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _shell_spec(tmp_path: Path, *, kind: str, sleep_s: int, marker_kind: str):
    """Build a JobSpec whose worker is `sleep N && write file` for testing."""
    out_file = tmp_path / f"{kind}.out"
    log_file = tmp_path / f"{kind}.log"
    if marker_kind == "log_contains":
        cmd = f"sleep {sleep_s} && echo 'total computation time' > {log_file} && echo done > {out_file}"
        marker = {"type": "log_contains", "pattern": "total computation time", "log_files": [str(log_file)]}
    elif marker_kind == "returncode_zero":
        cmd = f"sleep {sleep_s} && echo done > {out_file}"
        marker = {"type": "returncode_zero"}
    else:  # file_nonempty
        cmd = f"sleep {sleep_s} && echo done > {out_file}"
        marker = {"type": "file_nonempty"}
    return JobSpec(
        kind=kind,
        runner={"type": "shell", "cmd": ["/bin/sh", "-c", cmd]},
        output_dir=str(tmp_path),
        output_name=kind,
        expected_files=[str(out_file)],
        completion_marker=marker,
        eta_seconds=sleep_s,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_first_call_running_then_completed_via_wait(tmp_path: Path) -> None:
    spec = _shell_spec(tmp_path, kind="fake_a", sleep_s=2, marker_kind="file_nonempty")

    r1 = start_or_poll(spec)
    assert r1["status"] == "running"
    assert r1["kind"] == "fake_a"
    assert r1["do_not_bypass"] is True
    assert r1["recommended_poll_seconds"] >= 30

    res = wait_for_job(r1["job_file"], max_wait_seconds=15)
    assert res["status"] == "completed"
    assert res["result_files"] == [str(tmp_path / "fake_a.out")]
    assert (tmp_path / "fake_a.out").exists()


def test_repeat_call_dedup_does_not_spawn_second_worker(tmp_path: Path) -> None:
    spec = _shell_spec(tmp_path, kind="fake_b", sleep_s=3, marker_kind="file_nonempty")

    r1 = start_or_poll(spec)
    assert r1["status"] == "running"
    pid1 = r1["pid"]

    r2 = start_or_poll(spec)
    assert r2["status"] == "running"
    assert r2["pid"] == pid1, "second call must surface the same worker, not start a new one"

    res = wait_for_job(r1["job_file"], max_wait_seconds=15)
    assert res["status"] == "completed"


def test_poll_job_does_not_complete_while_worker_is_still_running(tmp_path: Path) -> None:
    out_file = tmp_path / "early_marker.out"
    spec = JobSpec(
        kind="early_marker",
        runner={"type": "shell", "cmd": ["/bin/sh", "-c", f"echo done > {out_file} && sleep 3"]},
        output_dir=str(tmp_path),
        output_name="early_marker",
        expected_files=[str(out_file)],
        completion_marker={"type": "file_nonempty"},
        eta_seconds=3,
    )

    r1 = start_or_poll(spec)
    deadline = time.time() + 5
    while not out_file.exists() and time.time() < deadline:
        time.sleep(0.05)

    assert out_file.exists()
    assert poll_job(r1["job_file"])["status"] == "running"
    assert wait_for_job(r1["job_file"], max_wait_seconds=10)["status"] == "completed"


def test_completed_call_is_idempotent(tmp_path: Path) -> None:
    spec = _shell_spec(tmp_path, kind="fake_c", sleep_s=1, marker_kind="file_nonempty")
    r1 = start_or_poll(spec)
    wait_for_job(r1["job_file"], max_wait_seconds=10)
    r3 = start_or_poll(spec)
    assert r3["status"] == "completed"
    assert r3["result_files"] == [str(tmp_path / "fake_c.out")]


def test_log_contains_marker(tmp_path: Path) -> None:
    spec = _shell_spec(tmp_path, kind="fake_log", sleep_s=1, marker_kind="log_contains")
    r1 = start_or_poll(spec)
    res = wait_for_job(r1["job_file"], max_wait_seconds=10)
    assert res["status"] == "completed", res


def test_failure_persists_returncode(tmp_path: Path) -> None:
    out_file = tmp_path / "fail.out"
    spec = JobSpec(
        kind="fake_fail",
        runner={"type": "shell", "cmd": ["/bin/sh", "-c", "exit 7"]},
        output_dir=str(tmp_path),
        output_name="fake_fail",
        expected_files=[str(out_file)],
        completion_marker={"type": "file_nonempty"},
        eta_seconds=1,
    )
    r1 = start_or_poll(spec)
    res = wait_for_job(r1["job_file"], max_wait_seconds=10)
    assert res["status"] == "failed"
    assert res["returncode"] == 7
    assert not out_file.exists()


def test_python_runner_dispatch(tmp_path: Path) -> None:
    """Python runner path: import a stdlib helper and verify completion via
    file_nonempty marker."""
    helper_dir = tmp_path / "helper_pkg"
    helper_dir.mkdir()
    helper_path = helper_dir / "_lj_helper.py"
    helper_path.write_text(
        "def make_file(path):\n"
        "    import pathlib\n"
        "    pathlib.Path(path).write_text('ok\\n')\n"
        "    return path\n",
        encoding="utf-8",
    )

    out_file = tmp_path / "py.out"
    spec = JobSpec(
        kind="fake_py",
        runner={
            "type": "python",
            "module": "_lj_helper",
            "function": "make_file",
            "kwargs": {"path": str(out_file)},
        },
        output_dir=str(tmp_path),
        output_name="fake_py",
        expected_files=[str(out_file)],
        completion_marker={"type": "file_nonempty"},
        eta_seconds=1,
    )

    # Worker runs in a subprocess — it needs ``helper_dir`` on PYTHONPATH.
    old_pp = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = str(helper_dir) + (os.pathsep + old_pp if old_pp else "")
    try:
        r1 = start_or_poll(spec)
        res = wait_for_job(r1["job_file"], max_wait_seconds=10)
    finally:
        if old_pp:
            os.environ["PYTHONPATH"] = old_pp
        else:
            os.environ.pop("PYTHONPATH", None)

    assert res["status"] == "completed", res
    assert out_file.exists()


def test_register_wait_for_job_attaches_tool_named_wait_for_job() -> None:
    class FakeMCP:
        def __init__(self):
            self.tools = {}

        def tool(self):
            def deco(fn):
                self.tools[fn.__name__] = fn
                return fn
            return deco

    mcp = FakeMCP()
    register_wait_for_job(mcp)
    assert "wait_for_job" in mcp.tools


def test_wait_for_job_caps_max_wait() -> None:
    """Block bound is hard-capped at WAIT_FOR_JOB_MAX_SECONDS even if a higher
    number is passed."""
    fake = "/nonexistent/job.json"
    started = time.time()
    res = wait_for_job(fake, max_wait_seconds=10_000)
    elapsed = time.time() - started
    assert res["status"] == "unknown"
    # Unknown returns immediately (no sleep loop), but make sure we didn't go
    # near the requested 10000 s anyway.
    assert elapsed < WAIT_FOR_JOB_MAX_SECONDS + 5


def test_corrupt_job_file_recovers(tmp_path: Path) -> None:
    """If a stale .mcp_job.json exists with a dead PID and no outputs, the
    next start_or_poll should mark it failed and surface the failure."""
    spec = _shell_spec(tmp_path, kind="fake_stale", sleep_s=1, marker_kind="file_nonempty")
    job_file = Path(spec.output_dir) / f"{spec.output_name}.{spec.kind}.mcp_job.json"
    job_file.parent.mkdir(parents=True, exist_ok=True)
    job_file.write_text(json.dumps({
        "job_id": "fake_stale:stale",
        "kind": "fake_stale",
        "status": "running",
        "pid": 999_999_999,  # very unlikely to exist
        "expected_files": spec.expected_files,
        "completion_marker": spec.completion_marker,
        "started_at": time.time() - 600,
    }), encoding="utf-8")

    res = start_or_poll(spec)
    assert res["status"] == "failed"


def test_stale_live_proc_without_pid_is_failed(tmp_path: Path, monkeypatch) -> None:
    """A stale Popen handle alone must not keep a job permanently running."""
    spec = _shell_spec(tmp_path, kind="fake_stale_live", sleep_s=1, marker_kind="file_nonempty")
    job_file = Path(spec.output_dir).resolve() / f"{spec.output_name}.{spec.kind}.mcp_job.json"
    job_file.parent.mkdir(parents=True, exist_ok=True)
    job_file.write_text(json.dumps({
        "job_id": "fake_stale_live:stale",
        "kind": "fake_stale_live",
        "status": "running",
        "pid": 999_999_999,
        "expected_files": spec.expected_files,
        "completion_marker": spec.completion_marker,
        "started_at": time.time() - 600,
    }), encoding="utf-8")

    class StaleProc:
        def poll(self):
            return None

    longjob._LIVE_PROCS[str(job_file)] = StaleProc()
    monkeypatch.setattr(longjob, "_pid_running", lambda pid: False)
    try:
        res = poll_job(str(job_file))
    finally:
        longjob._LIVE_PROCS.pop(str(job_file), None)

    assert res["status"] == "failed"
