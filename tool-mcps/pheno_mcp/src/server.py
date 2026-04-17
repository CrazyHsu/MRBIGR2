"""pheno_mcp — phenotype preprocessing MCP (8 tools)."""
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

from mrbigr.core import pheno  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def filter_abundance(d, abundance):
        """Filter features by minimum abundance threshold."""
        import pandas as pd
        return pheno.abundance_filter(pd.DataFrame(d), abundance).to_dict()

    @mcp.tool()
    def filter_missing(d, missing_ratio):
        """Filter features by missing ratio threshold."""
        import pandas as pd
        return pheno.missing_filter(pd.DataFrame(d), missing_ratio).to_dict()

    @mcp.tool()
    def scale_phenotype(d, method='zscore'):
        """Scale phenotype data (log2, log10, zscore, minmax, robust, boxcox)."""
        import pandas as pd
        return pheno.scale_wrapper(pd.DataFrame(d), method).to_dict()

    @mcp.tool()
    def impute_phenotype(d, method='mean'):
        """Impute missing phenotype values."""
        import pandas as pd
        return pheno.pheno_imputer(pd.DataFrame(d), method=method).to_dict()

    @mcp.tool()
    def remove_outliers(d, method='zscore'):
        """Remove outliers from phenotype data."""
        import pandas as pd
        return pheno.outlier(pd.DataFrame(d), method=method).to_dict()

    @mcp.tool()
    def run_blup(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUP via REML mixed model."""
        import pandas as pd
        return pheno.blup(pd.DataFrame(d), method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def run_blue(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUE via OLS fixed-effects model."""
        import pandas as pd
        return pheno.blue(pd.DataFrame(d), method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def correct_trait(pc, y):
        """Correct phenotype using principal components."""
        import pandas as pd
        import numpy as np
        return pheno.trait_correct(pd.DataFrame(pc), np.array(y)).tolist()


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("pheno_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
