"""mr_mcp — Mendelian randomization MCP (6 tools).

``run_mr_analysis`` and ``run_qtl_target_analysis`` are wrapped via
:mod:`mrbigr.mcp.longjob` because internal MR loops (and any underlying
gemma/plink invocations) can take minutes on real datasets. Both write
the result to a CSV under ``output_dir / output_name.mr.csv`` (or
``.qtl_target.csv``) for file-based completion detection.
"""
from __future__ import annotations

import sys
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

from mrbigr.core import mr  # noqa: E402
from mrbigr.mcp.longjob import JobSpec, register_wait_for_job, start_or_poll  # noqa: E402


def _gwas_frame(value):
    """Return an MR-ready GWAS DataFrame from a path or MCP JSON object."""
    import pandas as pd

    if isinstance(value, pd.DataFrame):
        df = value.copy()
    elif isinstance(value, (str, Path)):
        df = pd.read_csv(value, sep=None, engine="python")
    else:
        df = pd.DataFrame(value)

    rename = {}
    if "snp" in df.columns and "rs" not in df.columns:
        rename["snp"] = "rs"
    if "SNP" in df.columns and "rs" not in df.columns:
        rename["SNP"] = "rs"
    if "p_wald" in df.columns and "pvalue" not in df.columns:
        rename["p_wald"] = "pvalue"
    if "p_score" in df.columns and "pvalue" not in df.columns:
        rename["p_score"] = "pvalue"
    return df.rename(columns=rename)


def _instrument_snps(instrument_file):
    """Read SNP ids from a PLINK .clumped file or simple SNP table."""
    import pandas as pd

    if not instrument_file:
        return None
    inst = pd.read_csv(instrument_file, sep=r"\s+", engine="python")
    if not any(col in inst.columns for col in ("SNP", "rs", "snp")) and inst.shape[1] == 1:
        inst = pd.read_csv(instrument_file, sep=None, engine="python")
    for col in ("SNP", "rs", "snp"):
        if col in inst.columns:
            return inst[col].dropna().astype(str).tolist()
    if inst.shape[1] == 1:
        return inst.iloc[:, 0].dropna().astype(str).tolist()
    raise ValueError(f"instrument file missing SNP column: {instrument_file}")


def _filter_instruments(df, instrument_file):
    snps = _instrument_snps(instrument_file)
    if not snps:
        return df
    if "rs" not in df.columns:
        raise ValueError("GWAS table missing rs/SNP column for instrument filtering")
    keep = set(snps)
    out = df[df["rs"].astype(str).isin(keep)].copy()
    if out.empty:
        raise ValueError(f"no GWAS rows matched instruments from {instrument_file}")
    return out


def _persist_gwas(value, out_dir: Path, name: str, instrument_file=None) -> str:
    """Coerce ``value`` to an MR-ready DataFrame and write it as CSV."""
    df = _filter_instruments(_gwas_frame(value), instrument_file)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.csv"
    df.to_csv(path, index=False)
    return str(path)


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_mr_analysis(exposure_gwas, outcome_gwas, method='ivw', output_dir=None, output_name='mr', instrument_file=None):
        """Mendelian Randomization analysis (IVW, MR-Egger).

        Long-running: first call returns status='running' instantly; re-call
        with the same arguments to poll, or use wait_for_job(job_file). Result
        DataFrame is written to ``output_dir / output_name.mr.csv``.
        """
        out_dir = Path(output_dir or Path.cwd() / "mr_output").resolve()
        exp_csv = _persist_gwas(exposure_gwas, out_dir, f"{output_name}_exposure", instrument_file=instrument_file)
        out_csv = _persist_gwas(outcome_gwas, out_dir, f"{output_name}_outcome", instrument_file=instrument_file)
        result_csv = str(out_dir / f"{output_name}.mr.csv")
        spec = JobSpec(
            kind="mr_analysis",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "mr_analysis_runner",
                "kwargs": {
                    "exposure_gwas_csv": exp_csv,
                    "outcome_gwas_csv": out_csv,
                    "output_csv": result_csv,
                    "method": method,
                },
            },
            output_dir=str(out_dir),
            output_name=output_name,
            expected_files=[result_csv],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
            extras={"instrument_file": str(instrument_file) if instrument_file else None},
        )
        return start_or_poll(spec)

    @mcp.tool()
    def calculate_mr_causal_estimate(exposure_gwas, outcome_gwas, snp_col='rs', beta_col='beta', se_col='se', method='ivw', instrument_file=None):
        """Calculate MR causal estimate between two traits."""
        exp_df = _filter_instruments(_gwas_frame(exposure_gwas), instrument_file)
        out_df = _filter_instruments(_gwas_frame(outcome_gwas), instrument_file)
        result = mr.mr_causal_estimate(
            exp_df, out_df, snp_col=snp_col, beta_col=beta_col, se_col=se_col, method=method
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def test_mr_pleiotropy(exposure_gwas, outcome_gwas, instrument_file=None):
        """Test for horizontal pleiotropy (MR-Egger intercept)."""
        exp_df = _filter_instruments(_gwas_frame(exposure_gwas), instrument_file)
        out_df = _filter_instruments(_gwas_frame(outcome_gwas), instrument_file)
        return mr.test_pleiotropy(exp_df, out_df)

    @mcp.tool()
    def test_mr_heterogeneity(exposure_gwas, outcome_gwas, instrument_file=None):
        """Test for heterogeneity using IVW method."""
        exp_df = _filter_instruments(_gwas_frame(exposure_gwas), instrument_file)
        out_df = _filter_instruments(_gwas_frame(outcome_gwas), instrument_file)
        return mr.heterogeneity_test(exp_df, out_df)

    @mcp.tool()
    def run_qtl_target_analysis(qtl_df, tf_genes, target_genes, window=500000, output_dir=None, output_name='qtl_target'):
        """QTL targeting analysis - find which TFs target which genes.

        Long-running for large gene sets: first call returns status='running';
        re-call with the same arguments to poll, or use wait_for_job. Result
        DataFrame written to ``output_dir / output_name.qtl_target.csv``.
        """
        import pandas as pd
        out_dir = Path(output_dir or Path.cwd() / "qtl_target_output").resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        if isinstance(qtl_data, (str, Path)):
            qtl_csv = str(qtl_data)
        else:
            qtl_csv = str(out_dir / f"{output_name}_qtl.csv")
            (qtl_data if hasattr(qtl_data, "to_csv") else pd.DataFrame(qtl_data)).to_csv(qtl_csv, index=False)
        result_csv = str(out_dir / f"{output_name}.qtl_target.csv")
        spec = JobSpec(
            kind="qtl_target",
            runner={
                "type": "python",
                "module": "mrbigr.mcp.runners",
                "function": "qtl_target_runner",
                "kwargs": {
                    "qtl_csv": qtl_csv,
                    "tf_genes": list(tf_genes) if not isinstance(tf_genes, str) else tf_genes,
                    "target_genes": list(target_genes) if not isinstance(target_genes, str) else target_genes,
                    "output_csv": result_csv,
                    "window": int(window),
                },
            },
            output_dir=str(out_dir),
            output_name=output_name,
            expected_files=[result_csv],
            completion_marker={"type": "file_nonempty"},
            eta_seconds=600,
        )
        return start_or_poll(spec)

    @mcp.tool()
    def format_qtl_for_mr(qtl_file, output_file=None):
        """Format QTL file for MR analysis."""
        result = mr.format_qtl_for_mr(qtl_file, output_file=output_file)
        return result.to_dict() if result is not None else None

    register_wait_for_job(mcp)


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("mr_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
