"""mr_mcp — Mendelian randomization MCP (6 tools)."""
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


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def run_mr_analysis(exposure_gwas, outcome_gwas, method='ivw'):
        """Mendelian Randomization analysis (IVW, MR-Egger)."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        result = mr.mr_analysis(exp_df, out_df, method=method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def calculate_mr_causal_estimate(exposure_gwas, outcome_gwas, snp_col='rs', beta_col='beta', se_col='se', method='ivw'):
        """Calculate MR causal estimate between two traits."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        result = mr.mr_causal_estimate(
            exp_df, out_df, snp_col=snp_col, beta_col=beta_col, se_col=se_col, method=method
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def test_mr_pleiotropy(exposure_gwas, outcome_gwas):
        """Test for horizontal pleiotropy (MR-Egger intercept)."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        return mr.test_pleiotropy(exp_df, out_df)

    @mcp.tool()
    def test_mr_heterogeneity(exposure_gwas, outcome_gwas):
        """Test for heterogeneity using IVW method."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        return mr.heterogeneity_test(exp_df, out_df)

    @mcp.tool()
    def run_qtl_target_analysis(qtl_df, tf_genes, target_genes, window=500000):
        """QTL targeting analysis - find which TFs target which genes."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        result = mr.qtl_target_analysis(qtl_data, tf_genes, target_genes, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def format_qtl_for_mr(qtl_file, output_file=None):
        """Format QTL file for MR analysis."""
        result = mr.format_qtl_for_mr(qtl_file, output_file=output_file)
        return result.to_dict() if result is not None else None


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("mr_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
