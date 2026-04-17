"""gwas_mcp — association analysis MCP (12 tools)."""
from __future__ import annotations

import os
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

from mrbigr.core import gwas, vis  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    def _auto_plot(result_files):
        plots = []
        for f in (result_files or []):
            if not f or not os.path.exists(f):
                continue
            stem = f.replace('.assoc.txt', '').replace('.assoc.linear', '')
            mp = vis.manhattan_plot(f, output_file=f"{stem}_manhattan.png")
            qp = vis.qq_plot(f, output_file=f"{stem}_qq.png")
            plots.append({"file": f, "manhattan": mp, "qq": qp})
        return plots

    @mcp.tool()
    def run_gwas_lm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False):
        """Linear Model GWAS using GEMMA. Set auto_plot=True to generate Manhattan and QQ plots."""
        import pandas as pd
        result_files = gwas.gwas_lm(pd.DataFrame(phe), geno_prefix, output_name=output_name, output_dir=output_dir)
        if auto_plot and result_files:
            return {"gwas_results": result_files, "plots": _auto_plot(result_files)}
        return result_files

    @mcp.tool()
    def run_gwas_lmm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False):
        """Linear Mixed Model GWAS using GEMMA. Set auto_plot=True to generate Manhattan and QQ plots."""
        import pandas as pd
        result_files = gwas.gwas_lmm(pd.DataFrame(phe), geno_prefix, output_name=output_name, output_dir=output_dir)
        if auto_plot and result_files:
            return {"gwas_results": result_files, "plots": _auto_plot(result_files)}
        return result_files

    @mcp.tool()
    def run_gwas_plink(phe, geno_prefix, pheno_col=None, output_name="gwas", auto_plot=False):
        """GWAS using PLINK linear regression. Set auto_plot=True to generate Manhattan and QQ plots."""
        import pandas as pd
        result_file = gwas.gwas_plink(pd.DataFrame(phe), geno_prefix, pheno_col=pheno_col, output_name=output_name)
        if auto_plot and result_file:
            return {"gwas_results": [result_file], "plots": _auto_plot([result_file])}
        return result_file

    @mcp.tool()
    def add_rs_id_to_vcf(vcf_file, output_file=None):
        """Add rs IDs to VCF file if SNP IDs are missing."""
        return gwas.add_rs_id_to_vcf(vcf_file, output_file=output_file)

    @mcp.tool()
    def ensure_rs_id(geno_prefix):
        """Ensure PLINK bim file has rs IDs for SNPs."""
        return gwas.ensure_rs_id(geno_prefix)

    @mcp.tool()
    def run_simple_gwas(Y, G, cov=None):
        """Simple linear regression GWAS (pure Python, no external tools)."""
        import numpy as np
        cov_arr = np.array(cov) if cov is not None else None
        result = gwas.simple_gwas(np.array(Y), np.array(G), cov=cov_arr)
        return result.to_dict()

    @mcp.tool()
    def generate_clump(gwas_dir):
        """Generate clump input files from GWAS results."""
        return gwas.generate_clump_input(gwas_dir)

    @mcp.tool()
    def run_gwas_clump(geno_prefix, p1=0.001, p2=0.05):
        """GWAS result clumping using PLINK."""
        return gwas.gwas_clump(geno_prefix, p1=p1, p2=p2)

    @mcp.tool()
    def get_top_snps(gwas_file, n=10):
        """Get top SNPs from GWAS results."""
        return gwas.get_top_snps(gwas_file, n=n).to_dict()

    @mcp.tool()
    def calculate_lambda(pvalues):
        """Calculate genomic inflation factor (lambda)."""
        import numpy as np
        return gwas.calculate_lambda(np.array(pvalues))

    @mcp.tool()
    def qq_plot_data(pvalues):
        """Generate QQ plot data from p-values."""
        import numpy as np
        return gwas.qq_plot(np.array(pvalues))

    @mcp.tool()
    def manhattan_plot_data(gwas_results):
        """Generate Manhattan plot data from GWAS results."""
        import pandas as pd
        return gwas.manhattan_plot(pd.DataFrame(gwas_results)).to_dict()


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("gwas_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
