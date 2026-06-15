"""net_mcp — network analysis MCP (2 tools)."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
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

from mrbigr.core import net  # noqa: E402
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _stable_hash(value) -> str:
    try:
        payload = json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)
    except TypeError:
        payload = repr(value)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


def _job_root() -> Path:
    cwd = Path.cwd()
    if os.access(cwd, os.W_OK):
        return cwd / ".mrbigr_mcp_jobs"
    return Path(tempfile.gettempdir()) / "mrbigr2_mcp_jobs"


def _module_paths(edge_weight):
    if isinstance(edge_weight, (str, Path)):
        edge_path = Path(edge_weight).expanduser().resolve()
        prefix = str(edge_path).replace(".csv", "").replace(".edge_list", "")
        output_csv = Path(f"{prefix}.cluster_one.result.csv").resolve()
        return str(edge_path), output_csv, output_csv.parent, output_csv.stem

    import pandas as pd

    job_dir = (_job_root() / "net_module_identify" / _stable_hash(edge_weight)).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)
    edge_csv = job_dir / "edge_weight.csv"
    edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
    if not isinstance(edge_df, pd.DataFrame):
        edge_df = pd.DataFrame(edge_df)
    edge_df.to_csv(edge_csv, index=False)
    return str(edge_csv), job_dir / "module_identify.csv", job_dir, "module_identify"


def _hub_paths(edge_weight, cluster_one_res):
    """Resolve (edge_csv, cluster_csv, output_cluster_csv, output_hub_csv,
    output_dir, output_name) for the hub_identify long-job.

    Inputs may be CSV paths, DataFrames, or dicts. In-memory inputs are
    serialized into a stable, content-hashed job dir under ``_job_root()``.
    """
    import pandas as pd

    job_dir = (_job_root() / "net_hub_identify" / _stable_hash([edge_weight, cluster_one_res])).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(edge_weight, (str, Path)):
        edge_csv = str(Path(edge_weight).expanduser().resolve())
    else:
        edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
        if not isinstance(edge_df, pd.DataFrame):
            edge_df = pd.DataFrame(edge_df)
        edge_csv = str(job_dir / "edge_weight.csv")
        edge_df.to_csv(edge_csv, index=False)

    if isinstance(cluster_one_res, (str, Path)):
        cluster_csv = str(Path(cluster_one_res).expanduser().resolve())
    else:
        cluster_df = pd.DataFrame(cluster_one_res) if isinstance(cluster_one_res, dict) else cluster_one_res
        if not isinstance(cluster_df, pd.DataFrame):
            cluster_df = pd.DataFrame(cluster_df)
        cluster_csv = str(job_dir / "cluster_one.csv")
        cluster_df.to_csv(cluster_csv, index=False)

    output_cluster_csv = job_dir / "hub_cluster.csv"
    output_hub_csv = job_dir / "hub_genes.csv"
    return edge_csv, cluster_csv, str(output_cluster_csv), str(output_hub_csv), str(job_dir), "hub_identify"


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def module_identify(edge_weight, module_size=3):
        """Identify network modules using ClusterONE or a NetworkX fallback.

        Long-running: first call starts the ClusterONE/NetworkX job; repeated
        calls with the same inputs poll until completed.
        """
        edge_csv, output_csv, output_dir, output_name = _module_paths(edge_weight)
        spec = JobSpec(
            kind="module_identify",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "net_module_identify_runner",
                "kwargs": {
                    "edge_weight_csv": edge_csv,
                    "output_csv": str(output_csv),
                    "module_size": int(module_size),
                },
            },
            output_dir=str(output_dir),
            output_name=output_name,
            expected_files=[str(output_csv)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed" and output_csv.is_file():
            import pandas as pd
            result["result"] = pd.read_csv(output_csv).to_dict()
            result["output_file"] = str(output_csv)
        return result

    @mcp.tool()
    def hub_identify(edge_weight, cluster_one_res):
        """Identify hub genes in each network module.

        Long-running: first call returns status='running'; re-call with the
        same arguments to poll, or use wait_for_job(job_file). Results are
        written to <job_dir>/hub_cluster.csv and <job_dir>/hub_genes.csv.
        """
        edge_csv, cluster_csv, out_cluster_csv, out_hub_csv, output_dir, output_name = _hub_paths(
            edge_weight, cluster_one_res
        )
        spec = JobSpec(
            kind="hub_identify",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "net_hub_identify_runner",
                "kwargs": {
                    "edge_weight_csv": edge_csv,
                    "cluster_one_csv": cluster_csv,
                    "output_cluster_csv": out_cluster_csv,
                    "output_hub_csv": out_hub_csv,
                },
            },
            output_dir=output_dir,
            output_name=output_name,
            expected_files=[out_cluster_csv],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed":
            import pandas as pd
            cluster_df = pd.read_csv(out_cluster_csv) if Path(out_cluster_csv).is_file() else None
            hub_df = pd.read_csv(out_hub_csv) if Path(out_hub_csv).is_file() else None
            result["cluster_result"] = cluster_df.to_dict() if cluster_df is not None else None
            result["hub_result"] = hub_df.to_dict() if hub_df is not None else None
            result["output_cluster_file"] = out_cluster_csv
            result["output_hub_file"] = out_hub_csv
        return result

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("net_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
