"""Shared long-running job framework for MRBIGR2 MCP servers.

Any tool whose underlying work takes more than a few seconds (gemma, plink,
ClusterONE, MR loops, large-GTF parsing, ...) should be wrapped via
``start_or_poll(JobSpec(...))``. The first call spawns a detached worker and
returns ``status="running"`` instantly; repeating the same call (same
``kind`` + ``output_dir`` + ``output_name``) polls the existing job.

Companion tool ``wait_for_job(job_file, max_wait_seconds)`` blocks
server-side up to 180 s and returns the latest status — registered in each
server via ``register_wait_for_job(mcp)``.

Worker entry: ``python -m mrbigr.mcp.longjob_worker <job_file>``. The job
JSON tells the worker whether to spawn a subprocess (``runner.type ==
"shell"``) or call into a Python adapter (``runner.type == "python"``).

Goal: every wrapped MCP tool returns the same response shape so all coding
agents (Claude Code, codex, opencode, cursor, aider, ...) see one contract.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

# Hard cap on the server-side blocking helper. Beyond ~180 s the MCP client's
# request timeout kicks in (codex CLI default 60 s requires the user to bump
# request_timeout_sec; see plan doc fourth layer).
WAIT_FOR_JOB_MAX_SECONDS = 180

# Live subprocess.Popen handles, keyed by the job_file path. Lets the same
# server process recognise its own children even before they finish writing
# the completion marker.
_LIVE_PROCS: dict[str, subprocess.Popen] = {}


# ---------------------------------------------------------------------------
# JobSpec
# ---------------------------------------------------------------------------


@dataclass
class JobSpec:
    """Description of a long-running job — fully JSON-serializable."""

    kind: str  # e.g. "gwas_lmm", "kinship", "mr_analysis"
    runner: dict  # see _spawn_worker for accepted shapes
    output_dir: str
    output_name: str
    expected_files: list[str]
    completion_marker: dict  # {"type": "file_nonempty"} | {"type": "log_contains", ...}
    eta_seconds: int = 0
    do_not_bypass_hint: str = (
        "Re-call this same tool with the same arguments to poll. "
        "Do not bypass via shell gemma/plink/java; the server already runs "
        "it detached and dedup-keys by (kind, output_dir, output_name)."
    )
    # Free-form extras callers want surfaced in the job JSON / responses.
    extras: dict = field(default_factory=dict)

    def as_json(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Filesystem / process helpers (migrated from gwas_mcp + geno_mcp)
# ---------------------------------------------------------------------------


def _job_file_for(spec: JobSpec) -> Path:
    out_dir = Path(spec.output_dir).resolve()
    return out_dir / f"{spec.output_name}.{spec.kind}.mcp_job.json"


def _pid_running(pid) -> bool:
    try:
        pid_int = int(pid)
        os.kill(pid_int, 0)
    except (OSError, ValueError, TypeError):
        return False
    stat_path = Path("/proc") / str(pid_int) / "stat"
    try:
        stat_fields = stat_path.read_text(errors="ignore").split()
    except OSError:
        return True
    if len(stat_fields) > 2 and stat_fields[2] == "Z":
        return False
    return True


def _read_job(job_file: Path) -> Optional[dict]:
    try:
        return json.loads(job_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_job(job_file: Path, data: dict) -> None:
    job_file.parent.mkdir(parents=True, exist_ok=True)
    job_file.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Completion markers
# ---------------------------------------------------------------------------


def _completion_check(spec_or_job) -> bool:
    """Return True if the job's expected output exists and matches its marker.

    Accepts either a JobSpec or a job dict (as loaded from the job file).
    """
    if isinstance(spec_or_job, JobSpec):
        expected = list(spec_or_job.expected_files)
        marker = dict(spec_or_job.completion_marker or {})
    else:
        expected = list(spec_or_job.get("expected_files") or [])
        marker = dict(spec_or_job.get("completion_marker") or {})

    if not expected:
        return False
    for path in expected:
        p = Path(path)
        if not (p.is_file() and p.stat().st_size > 0):
            return False

    mtype = (marker or {}).get("type", "file_nonempty")
    if mtype == "file_nonempty":
        return True
    if mtype == "log_contains":
        log_files = marker.get("log_files") or [
            str(Path(p).with_suffix("").with_suffix(".log.txt"))
            for p in expected
        ]
        pattern = marker.get("pattern", "")
        for lf in log_files:
            try:
                text = Path(lf).read_text(errors="ignore")
            except OSError:
                return False
            if pattern and pattern not in text:
                return False
        return True
    if mtype == "log_regex":
        log_files = marker.get("log_files") or [
            str(Path(p).with_suffix("").with_suffix(".log.txt"))
            for p in expected
        ]
        pattern = marker.get("pattern", "")
        rx = re.compile(pattern)
        for lf in log_files:
            try:
                text = Path(lf).read_text(errors="ignore")
            except OSError:
                return False
            if not rx.search(text):
                return False
        return True
    if mtype == "returncode_zero":
        # Marker only meaningful after process exits; checked elsewhere.
        return True
    return True


# ---------------------------------------------------------------------------
# Spawn / poll
# ---------------------------------------------------------------------------


def _mrbigr_src_parent() -> str:
    """Path that needs to be on PYTHONPATH for `import mrbigr` to work in
    the worker subprocess. Resolved from this file's location."""
    return str(Path(__file__).resolve().parent.parent.parent)  # .../src


def _spawn_worker(spec: JobSpec, job_file: Path, stdout_file: Path, stderr_file: Path) -> subprocess.Popen:
    """Spawn the worker subprocess. Returns Popen handle (already started)."""
    cmd = [sys.executable, "-m", "mrbigr.mcp.longjob_worker", str(job_file)]
    env = os.environ.copy()
    src_parent = _mrbigr_src_parent()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = src_parent + (os.pathsep + existing if existing else "")
    with stdout_file.open("w", encoding="utf-8") as so, stderr_file.open("w", encoding="utf-8") as se:
        proc = subprocess.Popen(
            cmd,
            stdout=so,
            stderr=se,
            text=True,
            start_new_session=True,
            env=env,
        )
    return proc


def _running_response(spec: JobSpec, job: dict) -> dict:
    return {
        "status": "running",
        "kind": spec.kind,
        "job_id": job.get("job_id"),
        "job_file": str(_job_file_for(spec)),
        "pid": job.get("pid"),
        "output_dir": str(Path(spec.output_dir).resolve()),
        "output_name": spec.output_name,
        "expected_files": list(spec.expected_files),
        "started_at": job.get("started_at"),
        "eta_seconds": spec.eta_seconds,
        "recommended_poll_seconds": _recommended_poll(spec.eta_seconds),
        "do_not_bypass": True,
        "hint": spec.do_not_bypass_hint,
        "stdout_file": job.get("stdout_file"),
        "stderr_file": job.get("stderr_file"),
    }


def _completed_response(spec: JobSpec, job: dict) -> dict:
    response = {
        "status": "completed",
        "kind": spec.kind,
        "job_id": job.get("job_id"),
        "job_file": str(_job_file_for(spec)),
        "output_dir": str(Path(spec.output_dir).resolve()),
        "output_name": spec.output_name,
        "result_files": list(spec.expected_files),
        "expected_files": list(spec.expected_files),
        "started_at": job.get("started_at"),
        "completed_at": job.get("completed_at"),
        "stdout_file": job.get("stdout_file"),
        "stderr_file": job.get("stderr_file"),
    }
    if "result" in job:
        response["result"] = job.get("result")
    return response


def _failed_response(spec: JobSpec, job: dict, returncode) -> dict:
    return {
        "status": "failed",
        "kind": spec.kind,
        "job_id": job.get("job_id"),
        "job_file": str(_job_file_for(spec)),
        "output_dir": str(Path(spec.output_dir).resolve()),
        "output_name": spec.output_name,
        "expected_files": list(spec.expected_files),
        "returncode": returncode,
        "started_at": job.get("started_at"),
        "completed_at": job.get("completed_at"),
        "error": job.get("error"),
        "stdout_file": job.get("stdout_file"),
        "stderr_file": job.get("stderr_file"),
    }


def _recommended_poll(eta_seconds: int) -> int:
    if not eta_seconds or eta_seconds <= 0:
        return 30
    return max(30, min(eta_seconds // 4, WAIT_FOR_JOB_MAX_SECONDS))


def start_or_poll(spec: JobSpec) -> dict:
    """First call: spawn detached worker, return ``status='running'`` instantly.
    Repeated call with the same dedup key: poll job state. Same response shape
    is reused across all wrapped MCP tools.
    """
    output_dir = Path(spec.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    spec.output_dir = str(output_dir)

    job_file = _job_file_for(spec)

    # Already complete? (idempotent re-call after success)
    if _completion_check(spec):
        existing = _read_job(job_file) or {}
        existing.setdefault("job_id", f"{spec.kind}:{output_dir}:{spec.output_name}")
        existing.setdefault("started_at", time.time())
        existing.setdefault("completed_at", time.time())
        return _completed_response(spec, existing)

    job = _read_job(job_file)
    if job:
        if job.get("status") == "failed":
            return _failed_response(spec, job, job.get("returncode"))
        if job.get("status") == "completed":
            marker_type = dict(job.get("completion_marker") or {}).get("type")
            if _completion_check(spec) or (marker_type == "returncode_zero" and job.get("returncode") == 0):
                return _completed_response(spec, job)
            job["status"] = "failed"
            job.setdefault("completed_at", time.time())
            job.setdefault("returncode", 1)
            job.setdefault("error", "job marked completed but expected output marker is missing")
            _write_job(job_file, job)
            return _failed_response(spec, job, job.get("returncode"))
        proc = _LIVE_PROCS.get(str(job_file))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and _pid_running(job.get("pid")):
            return _running_response(spec, job)

        # Process is gone. Re-check completion in case of race.
        if _completion_check(spec):
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_job(job_file, job)
            return _completed_response(spec, job)

        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_job(job_file, job)
        return _failed_response(spec, job, returncode)

    # New job: spawn worker.
    stdout_file = output_dir / f"{spec.output_name}.{spec.kind}.mcp_stdout.log"
    stderr_file = output_dir / f"{spec.output_name}.{spec.kind}.mcp_stderr.log"

    job_id = f"{spec.kind}:{output_dir}:{spec.output_name}"
    job_payload = {
        "job_id": job_id,
        "kind": spec.kind,
        "status": "running",
        "runner": spec.runner,
        "output_dir": str(output_dir),
        "output_name": spec.output_name,
        "expected_files": list(spec.expected_files),
        "completion_marker": dict(spec.completion_marker or {}),
        "do_not_bypass_hint": spec.do_not_bypass_hint,
        "eta_seconds": spec.eta_seconds,
        "extras": dict(spec.extras or {}),
        "started_at": time.time(),
        "stdout_file": str(stdout_file),
        "stderr_file": str(stderr_file),
    }
    _write_job(job_file, job_payload)

    proc = _spawn_worker(spec, job_file, stdout_file, stderr_file)
    job_payload["pid"] = proc.pid
    _LIVE_PROCS[str(job_file)] = proc
    _write_job(job_file, job_payload)
    return _running_response(spec, job_payload)


def poll_job(job_file: str | Path) -> dict:
    """Poll a job file directly (no JobSpec needed). Used by wait_for_job."""
    job_file = Path(job_file)
    job = _read_job(job_file)
    if not job:
        return {"status": "unknown", "job_file": str(job_file), "hint": "job file not found"}

    spec_for_check = JobSpec(
        kind=job.get("kind", "unknown"),
        runner=job.get("runner", {}),
        output_dir=job.get("output_dir", "."),
        output_name=job.get("output_name", ""),
        expected_files=list(job.get("expected_files") or []),
        completion_marker=dict(job.get("completion_marker") or {}),
        eta_seconds=int(job.get("eta_seconds") or 0),
        do_not_bypass_hint=job.get("do_not_bypass_hint", ""),
    )

    if job.get("status") == "failed":
        return _failed_response(spec_for_check, job, job.get("returncode"))
    if job.get("status") == "completed":
        marker_type = dict(job.get("completion_marker") or {}).get("type")
        if _completion_check(spec_for_check) or (marker_type == "returncode_zero" and job.get("returncode") == 0):
            return _completed_response(spec_for_check, job)
        job["status"] = "failed"
        job.setdefault("completed_at", time.time())
        job.setdefault("returncode", 1)
        job.setdefault("error", "job marked completed but expected output marker is missing")
        _write_job(job_file, job)
        return _failed_response(spec_for_check, job, job.get("returncode"))

    proc = _LIVE_PROCS.get(str(job_file))
    returncode = proc.poll() if proc is not None else None
    if returncode is None and _pid_running(job.get("pid")):
        return _running_response(spec_for_check, job)

    if _completion_check(spec_for_check):
        if job.get("status") != "completed":
            job["status"] = "completed"
            job.setdefault("completed_at", time.time())
            job["returncode"] = 0 if returncode is None else returncode
            _write_job(job_file, job)
        return _completed_response(spec_for_check, job)

    job["status"] = "failed"
    job["completed_at"] = time.time()
    job["returncode"] = returncode
    _write_job(job_file, job)
    return _failed_response(spec_for_check, job, returncode)


# ---------------------------------------------------------------------------
# wait_for_job — registered as an MCP tool by every server
# ---------------------------------------------------------------------------


def wait_for_job(job_file: str | Path, max_wait_seconds: int = WAIT_FOR_JOB_MAX_SECONDS) -> dict:
    """Block server-side up to ``max_wait_seconds`` (hard-capped at
    ``WAIT_FOR_JOB_MAX_SECONDS``=180), then return latest poll result. Polls
    every 2 s. Requires the calling MCP client's request timeout to be
    >= max_wait_seconds + 60 s — see plan doc."""
    try:
        max_wait = int(max_wait_seconds)
    except (TypeError, ValueError):
        max_wait = WAIT_FOR_JOB_MAX_SECONDS
    max_wait = max(0, min(max_wait, WAIT_FOR_JOB_MAX_SECONDS))

    deadline = time.time() + max_wait
    while True:
        result = poll_job(job_file)
        if result.get("status") in ("completed", "failed", "unknown"):
            return result
        if time.time() >= deadline:
            return result
        time.sleep(2)


def register_wait_for_job(mcp) -> None:
    """Attach the universal ``wait_for_job`` tool to the given FastMCP."""

    @mcp.tool()
    def wait_for_job(job_file: str, max_wait_seconds: int = WAIT_FOR_JOB_MAX_SECONDS):  # noqa: D401
        """Block server-side up to ``max_wait_seconds`` (default 180, capped
        at 180), then return the latest status of the job at ``job_file``.
        Use this instead of busy-polling the originating tool. Returns the
        same response shape as the originating tool's poll.

        IMPORTANT — your MCP client request timeout must be >= 240 s for
        this to be safe. See MRBIGR2 plan doc, fourth layer.
        """
        return _wait_for_job_impl(job_file, max_wait_seconds=max_wait_seconds)


# Internal alias so the tool function above can call our public helper
# without colliding with its own name.
_wait_for_job_impl = wait_for_job


# ---------------------------------------------------------------------------
# Public re-exports
# ---------------------------------------------------------------------------

__all__ = [
    "JobSpec",
    "start_or_poll",
    "poll_job",
    "wait_for_job",
    "register_wait_for_job",
    "WAIT_FOR_JOB_MAX_SECONDS",
]
