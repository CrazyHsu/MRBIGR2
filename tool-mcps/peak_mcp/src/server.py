"""peak_mcp — peak & haplotype visualization MCP (4 tools)."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path


def _bootstrap_syspath() -> None:
    here = Path(__file__).resolve()
    for candidate in here.parents:
        src_dir = candidate / "src" / "mrbigr"
        if src_dir.is_dir():
            src_parent = str(src_dir.parent)
            if src_parent not in sys.path:
                sys.path.insert(0, src_parent)
            return


_bootstrap_syspath()

from mrbigr.core import peak  # noqa: E402

_BOXPLOT_PROCS = {}


def _pid_running(pid):
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError, TypeError):
        return False
    return True


def _read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _boxplot_paths(output_dir):
    output_dir = Path(output_dir).resolve()
    return {
        "job": output_dir / "plot_qtl_boxplot.mcp_job.json",
        "result": output_dir / "plot_qtl_boxplot.mcp_result.json",
        "qtl_csv": output_dir / "plot_qtl_boxplot.mcp_qtl.csv",
        "stdout": output_dir / "plot_qtl_boxplot.mcp_stdout.log",
        "stderr": output_dir / "plot_qtl_boxplot.mcp_stderr.log",
    }


def _boxplot_completed(output_dir, result_file):
    result = _read_json(result_file)
    if not result:
        return None
    output_files = result.get("output_files", [])
    if output_files and all(Path(path).is_file() and Path(path).stat().st_size > 0 for path in output_files):
        return {
            "status": "completed",
            "result": output_files,
            "output_dir": str(Path(output_dir).resolve()),
            "summary_file": str(Path(output_dir).resolve() / "qtl_boxplot_summary.csv"),
        }
    return None


def _coerce_qtl_df(qtl_df):
    import pandas as pd
    if isinstance(qtl_df, dict):
        return pd.DataFrame(qtl_df)
    if isinstance(qtl_df, (str, Path)):
        return pd.read_csv(qtl_df)
    if isinstance(qtl_df, pd.DataFrame):
        return qtl_df
    return pd.DataFrame(qtl_df)


def _plot_qtl_boxplot_status(pheno_file, geno_prefix, qtl_df, output_dir=None, test_method="t-test"):
    if output_dir is None:
        output_dir = "boxplot_output"
    output_dir = str(Path(output_dir).resolve())
    paths = _boxplot_paths(output_dir)
    completed = _boxplot_completed(output_dir, paths["result"])
    if completed:
        return completed

    job = _read_json(paths["job"])
    if job:
        proc = _BOXPLOT_PROCS.get(str(paths["job"]))
        returncode = proc.poll() if proc is not None else None
        if returncode is None and (proc is not None or _pid_running(job.get("pid"))):
            return {
                "status": "running",
                "job_id": job.get("job_id"),
                "pid": job.get("pid"),
                "output_dir": output_dir,
                "started_at": job.get("started_at"),
            }
        completed = _boxplot_completed(output_dir, paths["result"])
        if completed:
            job["status"] = "completed"
            job["completed_at"] = time.time()
            job["returncode"] = 0 if returncode is None else returncode
            _write_json(paths["job"], job)
            return completed
        job["status"] = "failed"
        job["completed_at"] = time.time()
        job["returncode"] = returncode
        _write_json(paths["job"], job)
        return {
            "status": "failed",
            "result": False,
            "output_dir": output_dir,
            "returncode": returncode,
            "stdout_file": str(paths["stdout"]),
            "stderr_file": str(paths["stderr"]),
        }

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    qtl_data = _coerce_qtl_df(qtl_df)
    qtl_data.to_csv(paths["qtl_csv"], index=False)
    helper = (
        "import json, sys; "
        "from pathlib import Path; "
        "from mrbigr.core import peak; "
        "outputs = peak.plot_qtl_boxplot(sys.argv[1], sys.argv[2], sys.argv[3], output_dir=sys.argv[4], test_method=sys.argv[5]); "
        "Path(sys.argv[6]).write_text(json.dumps({'output_files': outputs}, ensure_ascii=False), encoding='utf-8')"
    )
    env = os.environ.copy()
    sys_path = [path for path in sys.path if path]
    if env.get("PYTHONPATH"):
        sys_path.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(sys_path)
    stdout_fh = paths["stdout"].open("w", encoding="utf-8")
    stderr_fh = paths["stderr"].open("w", encoding="utf-8")
    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            helper,
            pheno_file,
            geno_prefix,
            str(paths["qtl_csv"]),
            output_dir,
            str(test_method),
            str(paths["result"]),
        ],
        stdout=stdout_fh,
        stderr=stderr_fh,
        text=True,
        start_new_session=True,
        env=env,
    )
    stdout_fh.close()
    stderr_fh.close()
    job = {
        "job_id": f"plot_qtl_boxplot:{output_dir}",
        "status": "running",
        "pid": proc.pid,
        "pheno_file": pheno_file,
        "geno_prefix": geno_prefix,
        "output_dir": output_dir,
        "test_method": test_method,
        "started_at": time.time(),
    }
    _BOXPLOT_PROCS[str(paths["job"])] = proc
    _write_json(paths["job"], job)
    return {
        "status": "running",
        "job_id": job["job_id"],
        "pid": proc.pid,
        "output_dir": output_dir,
        "started_at": job["started_at"],
    }


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df, output_dir=None, test_method='t-test'):
        """Generate boxplots for QTL regions.

        Long-running mode: first call starts the plot job and returns
        status='running'. Repeating the same call polls until completed.
        """
        return _plot_qtl_boxplot_status(
            pheno_file, geno_prefix, qtl_df,
            output_dir=output_dir, test_method=test_method,
        )

    @mcp.tool()
    def plot_grouped_boxplot(data_dict, output_file=None, test_method='t-test'):
        """Generate grouped boxplot."""
        return peak.plot_grouped_boxplot(data_dict, output_file=output_file, test_method=test_method)

    @mcp.tool()
    def haplotype_test(pheno_df, geno_df, snp_id, test_method='t-test'):
        """Perform statistical test for haplotype effect."""
        result = peak.haplotype_test(pheno_df, geno_df, snp_id, test_method=test_method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix):
        """Generate multi-trait Manhattan plot with QTL regions."""
        return peak.multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("peak_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
