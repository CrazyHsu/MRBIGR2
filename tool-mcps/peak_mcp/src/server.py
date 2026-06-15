"""peak_mcp — peak & haplotype visualization MCP (4 tools)."""
from __future__ import annotations

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

from mrbigr.core import peak  # noqa: E402
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _coerce_qtl_df(qtl_df):
    import pandas as pd
    if isinstance(qtl_df, dict):
        return pd.DataFrame(qtl_df)
    if isinstance(qtl_df, (str, Path)):
        return pd.read_csv(qtl_df)
    if isinstance(qtl_df, pd.DataFrame):
        return qtl_df
    return pd.DataFrame(qtl_df)


def _read_json_file(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _default_boxplot_dir() -> Path:
    if os.access(Path.cwd(), os.W_OK):
        return Path("boxplot_output").resolve()
    return Path(tempfile.gettempdir()) / "mrbigr2_mcp_jobs" / "boxplot_output"


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df, output_dir=None, test_method='t-test'):
        """Generate boxplots for QTL regions.

        Long-running mode: first call starts the plot job and returns
        status='running'. Repeating the same call polls until completed.
        """
        if output_dir is None:
            output_dir = str(_default_boxplot_dir())
        output_dir = str(Path(output_dir).expanduser().resolve())
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        qtl_csv = Path(output_dir) / "plot_qtl_boxplot.mcp_qtl.csv"
        result_json = Path(output_dir) / "plot_qtl_boxplot.mcp_result.json"
        qtl_data = _coerce_qtl_df(qtl_df)
        qtl_data.to_csv(qtl_csv, index=False)
        spec = JobSpec(
            kind="plot_qtl_boxplot",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "peak_boxplot_runner",
                "kwargs": {
                    "pheno_file": pheno_file,
                    "geno_prefix": geno_prefix,
                    "qtl_csv": str(qtl_csv),
                    "output_dir": output_dir,
                    "test_method": test_method,
                    "result_json": str(result_json),
                },
            },
            output_dir=output_dir,
            output_name="plot_qtl_boxplot",
            expected_files=[str(result_json)],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed":
            payload = _read_json_file(result_json)
            result["result"] = payload.get("output_files", [])
            result["output_dir"] = output_dir
            result["summary_file"] = payload.get("summary_file")
            result["n_plots"] = payload.get("n_plots")
            result["reason"] = payload.get("reason")
        return result

    @mcp.tool()
    def plot_grouped_boxplot(data_dict, output_file=None, test_method='t-test'):
        """Generate grouped boxplot."""
        from mrbigr.core._argjson import maybe_json_loads
        return peak.plot_grouped_boxplot(maybe_json_loads(data_dict), output_file=output_file, test_method=test_method)

    @mcp.tool()
    def haplotype_test(pheno_df, geno_df, snp_id, test_method='t-test'):
        """Perform statistical test for haplotype effect."""
        from mrbigr.core._argjson import coerce_table
        result = peak.haplotype_test(coerce_table(pheno_df), coerce_table(geno_df), snp_id, test_method=test_method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix, file_format='png'):
        """Generate multi-trait Manhattan plot with QTL regions."""
        output_prefix = str(Path(output_prefix).expanduser().resolve())
        output_file = f"{output_prefix}_multi_trait.{file_format}"
        spec = JobSpec(
            kind="multi_trait_qtl_plot",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "core_call_runner",
                "kwargs": {
                    "module": "mrbigr.core.peak",
                    "function": "multi_trait_qtl_plot",
                    "kwargs": {
                        "gwas_dir": gwas_dir,
                        "qtl_file": qtl_file,
                        "output_prefix": output_prefix,
                        "file_format": file_format,
                    },
                },
            },
            output_dir=str(Path(output_prefix).parent),
            output_name=Path(output_prefix).name,
            expected_files=[output_file],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=300,
        )
        result = start_or_poll(spec)
        if result.get("status") == "completed":
            result["result"] = output_file
            result["output_file"] = output_file
        return result

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("peak_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
