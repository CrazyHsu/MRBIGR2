#!/usr/bin/env python3
"""
MRBIGR2 CLI — file-based command-line interface to all MRBIGR2 tools.

Designed so that Claude Code (or any shell) can call MRBIGR2 tools immediately
after installation, without waiting for the MCP server to be loaded on next restart.

Usage:
    python mrbigr_cli.py <tool> [options]
    python mrbigr_cli.py list              # show all available tools
    python mrbigr_cli.py help <tool>       # show help for a tool

Examples:
    python mrbigr_cli.py gwas_lmm --phe pheno.csv --geno data/chr_HAMP
    python mrbigr_cli.py snp_qc --input data/geno --output output/geno_qc --maf 0.05
    python mrbigr_cli.py plot_manhattan --gwas_file output/result.assoc.txt
    python mrbigr_cli.py get_top_snps --gwas_file output/result.assoc.txt --n 20
"""
import sys
import os
import json
import argparse
import importlib
import warnings

warnings.filterwarnings("ignore")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import pandas as pd
import numpy as np


class _LazyModule:
    """Import tool modules only when a command actually uses them."""

    def __init__(self, module_name):
        self._module_name = module_name
        self._module = None

    def _load(self):
        if self._module is None:
            self._module = importlib.import_module(self._module_name)
        return self._module

    def __getattr__(self, name):
        return getattr(self._load(), name)


pheno = _LazyModule("mrbigr.core.pheno")
geno = _LazyModule("mrbigr.core.geno")
gwas = _LazyModule("mrbigr.core.gwas")
vis = _LazyModule("mrbigr.core.vis")
anno = _LazyModule("mrbigr.core.anno")
qtl = _LazyModule("mrbigr.core.qtl")
multi = _LazyModule("mrbigr.core.parallel")  # internal helper, CLI-only surface
peak = _LazyModule("mrbigr.core.peak")
mr = _LazyModule("mrbigr.core.mr")
go = _LazyModule("mrbigr.core.go")
net = _LazyModule("mrbigr.core.net")


def _load_phe(path, cols=None):
    """Load phenotype file (CSV or gzipped CSV), return DataFrame with sample index."""
    df = pd.read_csv(path, compression='infer', index_col=0)
    df.index = df.index.astype(str)
    if cols is not None:
        cols_list = [c.strip() for c in cols.split(",")]
        df = df[cols_list]
    return df


def _load_phe_raw(path):
    """Load phenotype file without forcing the first column to be the index."""
    return pd.read_csv(path, compression='infer')


def _print_result(result, label="Result"):
    """Pretty-print a result for Claude to read."""
    if result is None:
        print(f"{label}: None")
    elif isinstance(result, pd.DataFrame):
        print(f"{label}: DataFrame ({result.shape[0]} rows x {result.shape[1]} cols)")
        print(result.to_string(max_rows=30, max_cols=15))
    elif isinstance(result, dict):
        print(json.dumps(result, indent=2, default=str, ensure_ascii=False))
    elif isinstance(result, (list, tuple)):
        for i, item in enumerate(result):
            print(f"  [{i}] {item}")
    else:
        print(f"{label}: {result}")


# =============================================================================
# Tool registry: each entry is (function, description, argparse-setup-function)
# =============================================================================

TOOLS = {}


def tool(name, desc):
    """Decorator to register a CLI tool."""
    def decorator(fn):
        TOOLS[name] = {"fn": fn, "desc": desc}
        return fn
    return decorator


# ======================== GENO ========================

@tool("snp_qc", "SNP quality control using PLINK")
def cmd_snp_qc(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output PLINK prefix")
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--missing_rate", type=float, default=0.2)
    p.add_argument("--mind", type=float, default=0.2)
    a = p.parse_args(args)
    result = geno.snp_qc(a.input, a.output, maf=a.maf, missing_rate=a.missing_rate, mind=a.mind)
    _print_result(result, "snp_qc")

@tool("subset_genotype", "Subset PLINK genotype data by chromosome or random proportion")
def cmd_subset_genotype(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output PLINK prefix")
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--chromosomes", help="Chromosome list, e.g. 1,2,3")
    mode.add_argument("--proportion", type=float, help="Random SNP retention proportion in (0,1]")
    p.add_argument("--seed", type=int, default=42, help="Random seed for proportion-based subset")
    a = p.parse_args(args)
    result = geno.subset_plink(a.input, a.output, chromosomes=a.chromosomes, proportion=a.proportion, seed=a.seed)
    _print_result(result, "subset_genotype")

@tool("genotype_pca", "PCA analysis for genotype data")
def cmd_genotype_pca(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--n_components", type=int, default=10)
    p.add_argument("--max_snps", type=int, default=20000)
    p.add_argument("--output", default=None, help="Save PCA to CSV")
    a = p.parse_args(args)
    pc_df, var_ratio = geno.calculate_pca(a.input, n_components=a.n_components, max_snps=a.max_snps)
    print(f"PCA: {pc_df.shape[0]} samples x {len(var_ratio)} components")
    print(f"Variance explained: {[f'{v:.4f}' for v in var_ratio]}")
    if a.output:
        pc_df.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        print(pc_df.head(10).to_string())

@tool("calculate_kinship", "Calculate kinship matrix using GEMMA")
def cmd_calculate_kinship(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output prefix")
    a = p.parse_args(args)
    result = geno.calculate_kinship(a.input, a.output)
    _print_result(result, "kinship")

@tool("calculate_ibd", "Calculate Identity by Descent (IBD) matrix")
def cmd_calculate_ibd(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output prefix")
    a = p.parse_args(args)
    result = geno.calculate_ibd(a.input, a.output)
    _print_result(result, "IBD")

@tool("snp_pruning", "LD-based SNP pruning using PLINK")
def cmd_snp_pruning(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--window", type=int, default=50)
    p.add_argument("--shift", type=int, default=5)
    p.add_argument("--r2", type=float, default=0.5)
    p.add_argument("--maf", type=float, default=0.05)
    a = p.parse_args(args)
    result = geno.snp_pruning(a.input, a.output, window=a.window, shift=a.shift, r2=a.r2, maf=a.maf)
    _print_result(result, "snp_pruning")

@tool("snp_clumping", "Genotype-based LD clumping (greedy, like bigsnpr::snp_clumping)")
def cmd_snp_clumping(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output prefix")
    p.add_argument("--r2", type=float, default=0.5)
    p.add_argument("--maf", type=float, default=0.05)
    p.add_argument("--window_kb", type=int, default=250, help="LD window in kb")
    a = p.parse_args(args)
    result = geno.snp_clumping(a.input, a.output, r2=a.r2, maf=a.maf, window_kb=a.window_kb)
    _print_result(result, "snp_clumping")

@tool("snp_impute", "Impute missing genotype values")
def cmd_snp_impute(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--method", default="mean")
    a = p.parse_args(args)
    result = geno.snp_impute(a.input, a.output, method=a.method)
    _print_result(result, "snp_impute")

@tool("snp_stats", "Get basic SNP statistics")
def cmd_snp_stats(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    a = p.parse_args(args)
    result = geno.get_snp_stats(a.input)
    _print_result(result, "snp_stats")

@tool("convert_vcf", "Convert VCF to PLINK format")
def cmd_convert_vcf(args):
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(args)
    result = geno.vcf_to_plink(a.vcf, a.output)
    _print_result(result, "convert_vcf")

@tool("convert_hapmap", "Convert HapMap to PLINK format")
def cmd_convert_hapmap(args):
    p = argparse.ArgumentParser()
    p.add_argument("--hapmap", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(args)
    result = geno.hapmap_to_plink(a.hapmap, a.output)
    _print_result(result, "convert_hapmap")

@tool("plink_to_vcf", "Convert PLINK to VCF format")
def cmd_plink_to_vcf(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(args)
    result = geno.plink_to_vcf(a.input, a.output)
    _print_result(result, "plink_to_vcf")

@tool("build_tree", "Build phylogenetic tree using FastTree")
def cmd_build_tree(args):
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Input PLINK prefix")
    p.add_argument("--output", required=True, help="Output prefix")
    a = p.parse_args(args)
    result = geno.build_phylogenetic_tree(a.input, a.output)
    _print_result(result, "build_tree")


# ======================== PHENO ========================

@tool("filter_abundance", "Filter features by minimum abundance threshold")
def cmd_filter_abundance(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--abundance", type=float, required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe)
    result = pheno.abundance_filter(df, a.abundance)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "filter_abundance")

@tool("filter_missing", "Filter features by missing ratio threshold")
def cmd_filter_missing(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--missing_ratio", type=float, required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe)
    result = pheno.missing_filter(df, a.missing_ratio)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "filter_missing")

@tool("scale_phenotype", "Scale phenotype data (log2, log10, zscore, minmax, robust, boxcox)")
def cmd_scale_phenotype(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--method", default="zscore")
    p.add_argument("--output", default=None)
    p.add_argument("--cols", default=None, help="Comma-separated column names")
    a = p.parse_args(args)
    df = _load_phe(a.phe, cols=a.cols)
    result = pheno.scale_wrapper(df, a.method)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "scale_phenotype")

@tool("impute_phenotype", "Impute missing phenotype values")
def cmd_impute_phenotype(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--method", default="mean")
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe)
    result = pheno.pheno_imputer(df, method=a.method)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "impute_phenotype")

@tool("remove_outliers", "Remove outliers from phenotype data")
def cmd_remove_outliers(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--method", default="zscore")
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe)
    result = pheno.outlier(df, method=a.method)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "remove_outliers")

@tool("blup", "BLUP via REML mixed model (matrix, long, or legacy)")
def cmd_blup(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--method", choices=["matrix", "long", "legacy"], default="matrix")
    p.add_argument("--y_col", default="y")
    p.add_argument("--geno_col", default="genotype")
    p.add_argument("--env_col", default="env")
    p.add_argument("--rep_col", default=None)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe) if a.method in ("matrix", "legacy") else _load_phe_raw(a.phe)
    result = pheno.blup(df, method=a.method, y_col=a.y_col, geno_col=a.geno_col, env_col=a.env_col, rep_col=a.rep_col)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "blup")

@tool("blue", "BLUE via OLS fixed-effects model (matrix, long, or legacy)")
def cmd_blue(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--method", choices=["matrix", "long", "legacy"], default="matrix")
    p.add_argument("--y_col", default="y")
    p.add_argument("--geno_col", default="genotype")
    p.add_argument("--env_col", default="env")
    p.add_argument("--rep_col", default=None)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    df = _load_phe(a.phe) if a.method in ("matrix", "legacy") else _load_phe_raw(a.phe)
    result = pheno.blue(df, method=a.method, y_col=a.y_col, geno_col=a.geno_col, env_col=a.env_col, rep_col=a.rep_col)
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "blue")

@tool("trait_correct", "Correct phenotype for population structure using PCA")
def cmd_trait_correct(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True, help="Phenotype CSV file")
    p.add_argument("--pc", required=True, help="PCA CSV file (from genotype_pca)")
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    phe_df = _load_phe(a.phe)
    pc_df = pd.read_csv(a.pc)
    if 'sample' in pc_df.columns:
        pc_df = pc_df.set_index('sample')
    else:
        pc_df = pc_df.set_index(pc_df.columns[0])
    pc_df.index = pc_df.index.astype(str)
    pc_cols = [c for c in pc_df.columns if c.startswith('PC')]
    pc_df = pc_df[pc_cols].astype(float)
    common = phe_df.index.intersection(pc_df.index)
    result = pheno.trait_correct(pc_df.loc[common], phe_df.loc[common])
    if a.output:
        result.to_csv(a.output)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "trait_correct")


# ======================== GWAS ========================

@tool("gwas_lm", "Linear Model GWAS using GEMMA")
def cmd_gwas_lm(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True, help="Phenotype CSV file")
    p.add_argument("--geno", required=True, help="PLINK genotype prefix")
    p.add_argument("--cols", default=None, help="Comma-separated phenotype column names")
    p.add_argument("--output_name", default=None)
    a = p.parse_args(args)
    phe_df = _load_phe(a.phe, cols=a.cols)
    results = gwas.gwas_lm(phe_df, a.geno, output_name=a.output_name)
    _print_result(results, "gwas_lm output files")

@tool("gwas_lmm", "Linear Mixed Model GWAS using GEMMA (recommended)")
def cmd_gwas_lmm(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True, help="Phenotype CSV file")
    p.add_argument("--geno", required=True, help="PLINK genotype prefix")
    p.add_argument("--cols", default=None, help="Comma-separated phenotype column names")
    p.add_argument("--output_name", default=None)
    p.add_argument("--threads", type=int, default=None, help="Max parallel GEMMA processes")
    a = p.parse_args(args)
    phe_df = _load_phe(a.phe, cols=a.cols)
    results = gwas.gwas_lmm(phe_df, a.geno, output_name=a.output_name, num_threads=a.threads)
    _print_result(results, "gwas_lmm output files")

@tool("gwas_plink", "GWAS using PLINK linear regression")
def cmd_gwas_plink(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--geno", required=True)
    p.add_argument("--pheno_col", default=None)
    p.add_argument("--output_name", default="gwas")
    a = p.parse_args(args)
    phe_df = _load_phe(a.phe)
    result = gwas.gwas_plink(phe_df, a.geno, pheno_col=a.pheno_col, output_name=a.output_name)
    _print_result(result, "gwas_plink")

@tool("get_top_snps", "Get top SNPs from GWAS results")
def cmd_get_top_snps(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True, help="GWAS assoc.txt file")
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--pvalue_cutoff", type=float, default=None)
    a = p.parse_args(args)
    result = gwas.get_top_snps(a.gwas_file, n=a.n, pvalue_cutoff=a.pvalue_cutoff)
    _print_result(result, "top_snps")

@tool("calculate_lambda", "Calculate genomic inflation factor (lambda)")
def cmd_calculate_lambda(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--pval_col", default=None)
    a = p.parse_args(args)
    df = pd.read_csv(a.gwas_file, sep='\t')
    pval_col = a.pval_col
    if pval_col is None:
        for col in ['p_wald', 'p_lrt', 'p_score', 'P', 'pvalue']:
            if col in df.columns:
                pval_col = col
                break
    lam = gwas.calculate_lambda(df[pval_col].dropna().values)
    print(f"Genomic inflation factor (lambda): {lam:.4f}")

@tool("ensure_rs_id", "Ensure PLINK bim file has rs IDs for SNPs")
def cmd_ensure_rs_id(args):
    p = argparse.ArgumentParser()
    p.add_argument("--geno", required=True, help="PLINK prefix")
    a = p.parse_args(args)
    result = gwas.ensure_rs_id(a.geno)
    print(f"rs IDs added: {result}")

@tool("add_rs_id_to_vcf", "Add rs IDs to VCF file if missing")
def cmd_add_rs_id_to_vcf(args):
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = gwas.add_rs_id_to_vcf(a.vcf, output_file=a.output)
    print(f"Output: {result}")

@tool("generate_clump", "Generate clump input files from GWAS results")
def cmd_generate_clump(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_dir", required=True)
    a = p.parse_args(args)
    result = gwas.generate_clump_input(a.gwas_dir)
    print(f"Clump input dir: {result}")

@tool("plink_clump", "PLINK clumping for GWAS results")
def cmd_plink_clump(args):
    p = argparse.ArgumentParser()
    p.add_argument("--geno", required=True)
    p.add_argument("--p1", type=float, default=0.05)
    p.add_argument("--p2", type=float, default=0.001)
    a = p.parse_args(args)
    result = gwas.plink_clump(a.geno, p1=a.p1, p2=a.p2)
    print(f"Clump result dir: {result}")


# ======================== VIS ========================

@tool("plot_manhattan", "Generate Manhattan plot from GWAS results")
def cmd_plot_manhattan(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--significance", type=float, default=5e-8)
    p.add_argument("--suggest", type=float, default=1e-5)
    a = p.parse_args(args)
    result = vis.manhattan_plot(a.gwas_file, output_file=a.output,
                                significance=a.significance, suggest=a.suggest)
    _print_result(result, "manhattan_plot")

@tool("plot_qq", "Generate Q-Q plot from GWAS results")
def cmd_plot_qq(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = vis.qq_plot(a.gwas_file, output_file=a.output)
    _print_result(result, "qq_plot")

@tool("plot_pca", "Generate PCA scatter plot")
def cmd_plot_pca(args):
    p = argparse.ArgumentParser()
    p.add_argument("--pca_file", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--pc_x", default="PC1")
    p.add_argument("--pc_y", default="PC2")
    a = p.parse_args(args)
    result = vis.pca_plot(a.pca_file, output_file=a.output, pc_x=a.pc_x, pc_y=a.pc_y)
    _print_result(result, "pca_plot")

@tool("plot_ld_heatmap", "Generate LD heatmap from genotype data")
def cmd_plot_ld_heatmap(args):
    p = argparse.ArgumentParser()
    p.add_argument("--geno", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--max_snps", type=int, default=500)
    a = p.parse_args(args)
    result = vis.ld_heatmap(a.geno, output_file=a.output, max_snps=a.max_snps)
    _print_result(result, "ld_heatmap")

@tool("plot_pheno_hist", "Generate histogram for phenotype distribution")
def cmd_plot_pheno_hist(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = vis.phenotype_hist(a.phe, output_file=a.output)
    _print_result(result, "pheno_hist")

@tool("plot_pheno_boxplot", "Generate boxplot for phenotypes")
def cmd_plot_pheno_boxplot(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = vis.phenotype_boxplot(a.phe, output_file=a.output)
    _print_result(result, "pheno_boxplot")

@tool("plot_pheno_correlation", "Generate phenotype correlation heatmap")
def cmd_plot_pheno_correlation(args):
    p = argparse.ArgumentParser()
    p.add_argument("--phe", required=True)
    p.add_argument("--output", default=None)
    p.add_argument("--method", default="pearson")
    a = p.parse_args(args)
    result = vis.phenotype_correlation(a.phe, output_file=a.output, method=a.method)
    _print_result(result, "pheno_correlation")

@tool("gwas_summary", "Generate comprehensive GWAS summary plots (Manhattan + QQ)")
def cmd_gwas_summary(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_dir", required=True)
    p.add_argument("--output_prefix", default="gwas_summary")
    a = p.parse_args(args)
    result = vis.gwas_summary_plot(a.gwas_dir, output_prefix=a.output_prefix)
    _print_result(result, "gwas_summary")


# ======================== ANNO ========================

@tool("parse_gtf", "Parse GTF annotation file")
def cmd_parse_gtf(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = anno.parse_gtf(a.gtf)
    if result is not None and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "parse_gtf")

@tool("annotate_snps", "Annotate SNPs using GTF file")
def cmd_annotate_snps(args):
    p = argparse.ArgumentParser()
    p.add_argument("--snps", required=True, help="Comma-separated SNP list (chr:pos) or file")
    p.add_argument("--gtf", required=True)
    p.add_argument("--window", type=int, default=5000)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    if os.path.isfile(a.snps):
        snp_list = pd.read_csv(a.snps, header=None)[0].tolist()
    else:
        snp_list = a.snps.split(",")
    result = anno.annotate_snps_simple(snp_list, a.gtf, window=a.window)
    if result is not None and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "annotate_snps")

@tool("get_genes_in_region", "Get all genes in a genomic region")
def cmd_get_genes_in_region(args):
    p = argparse.ArgumentParser()
    p.add_argument("--chr", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--gtf", required=True)
    a = p.parse_args(args)
    result = anno.get_genes_in_region(getattr(a, 'chr'), a.start, a.end, a.gtf)
    _print_result(result, "genes_in_region")

@tool("qtl_annotation", "Annotate QTL regions with genes")
def cmd_qtl_annotation(args):
    p = argparse.ArgumentParser()
    p.add_argument("--qtl_file", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = anno.qtl_annotation(a.qtl_file, a.annotation)
    if result is not None and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        _print_result(result, "qtl_annotation")


# ======================== QTL ========================

@tool("detect_qtl", "Detect QTL regions from GWAS results")
def cmd_detect_qtl(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--p1", type=float, default=1e-7)
    p.add_argument("--p2", type=float, default=1e-5)
    p.add_argument("--window", type=int, default=500000)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = qtl.detect_qtl(a.gwas_file, p1=a.p1, p2=a.p2, window=a.window)
    if result is not None and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output} ({result.shape})")
    else:
        _print_result(result, "detect_qtl")

@tool("get_lead_snp", "Get lead SNP in a genomic region")
def cmd_get_lead_snp(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--chr", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    a = p.parse_args(args)
    result = qtl.get_lead_snp(a.gwas_file, getattr(a, 'chr'), a.start, a.end)
    _print_result(result, "lead_snp")

@tool("identify_peak_snps", "Identify peak SNPs from GWAS results")
def cmd_identify_peak_snps(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_dir", required=True)
    p.add_argument("--p_threshold", type=float, default=1e-5)
    p.add_argument("--max_peaks", type=int, default=50)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = qtl.identify_peak_snps(a.gwas_dir, p_threshold=a.p_threshold, max_peaks=a.max_peaks)
    if result is not None and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        _print_result(result, "peak_snps")

@tool("plot_qtl_region", "Plot GWAS results for a QTL region")
def cmd_plot_qtl_region(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_file", required=True)
    p.add_argument("--chr", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = qtl.plot_qtl_region(a.gwas_file, getattr(a, 'chr'), a.start, a.end, output_file=a.output)
    _print_result(result, "plot_qtl_region")

@tool("qtl_summary", "Generate summary of all QTL regions")
def cmd_qtl_summary(args):
    p = argparse.ArgumentParser()
    p.add_argument("--gwas_dir", required=True)
    p.add_argument("--output_prefix", default="qtl_summary")
    a = p.parse_args(args)
    result = qtl.qtl_summary(a.gwas_dir, output_prefix=a.output_prefix)
    _print_result(result, "qtl_summary")


# ======================== MR ========================

@tool("mr_analysis", "Mendelian Randomization analysis (IVW, MR-Egger)")
def cmd_mr_analysis(args):
    p = argparse.ArgumentParser()
    p.add_argument("--exposure", required=True, help="Exposure GWAS file")
    p.add_argument("--outcome", required=True, help="Outcome GWAS file")
    p.add_argument("--method", default="ivw")
    a = p.parse_args(args)
    exp_df = pd.read_csv(a.exposure, sep='\t')
    out_df = pd.read_csv(a.outcome, sep='\t')
    result = mr.mr_analysis(exp_df, out_df, method=a.method)
    _print_result(result, "mr_analysis")

@tool("mr_causal_estimate", "Calculate MR causal estimate between two traits")
def cmd_mr_causal_estimate(args):
    p = argparse.ArgumentParser()
    p.add_argument("--exposure", required=True)
    p.add_argument("--outcome", required=True)
    a = p.parse_args(args)
    exp_df = pd.read_csv(a.exposure, sep='\t')
    out_df = pd.read_csv(a.outcome, sep='\t')
    result = mr.mr_causal_estimate(exp_df, out_df)
    _print_result(result, "mr_causal_estimate")

@tool("format_qtl_for_mr", "Format QTL file for MR analysis")
def cmd_format_qtl_for_mr(args):
    p = argparse.ArgumentParser()
    p.add_argument("--qtl_file", required=True)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    result = mr.format_qtl_for_mr(a.qtl_file, output_file=a.output)
    _print_result(result, "format_qtl_for_mr")


# ======================== GO ========================

@tool("go_enrichment", "GO enrichment analysis (BP, MF, CC)")
def cmd_go_enrichment(args):
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True, help="Comma-separated gene list or file")
    p.add_argument("--organism", default="Mouse")
    p.add_argument("--pvalue_cutoff", type=float, default=0.05)
    p.add_argument("--qvalue_cutoff", type=float, default=0.05)
    p.add_argument("--mode", default="local", choices=["auto", "online", "local"])
    p.add_argument("--gene_sets_file", default=None, help="Local two-column gene2GO file or GMT file")
    p.add_argument("--background_genes", default=None, help="Comma-separated background gene list or file")
    p.add_argument("--go_obo_file", default=None, help="Optional local go-basic.obo path")
    p.add_argument("--obo_cache_dir", default=None, help="Directory for cached/downloaded go-basic.obo")
    p.add_argument("--no_auto_download_obo", action="store_true", help="Disable auto-download of go-basic.obo")
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    if os.path.isfile(a.genes):
        genes = pd.read_csv(a.genes, header=None)[0].tolist()
    else:
        genes = a.genes.split(",")
    if a.background_genes:
        if os.path.isfile(a.background_genes):
            background_genes = pd.read_csv(a.background_genes, header=None)[0].tolist()
        else:
            background_genes = a.background_genes.split(",")
    else:
        background_genes = None
    result = go.go_enrich(
        genes,
        organism=a.organism,
        pvalue_cutoff=a.pvalue_cutoff,
        qvalue_cutoff=a.qvalue_cutoff,
        mode=a.mode,
        gene_sets_file=a.gene_sets_file,
        background_genes=background_genes,
        go_obo_file=a.go_obo_file,
        obo_cache_dir=a.obo_cache_dir,
        auto_download_obo=not a.no_auto_download_obo,
    )
    if isinstance(result, pd.DataFrame) and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        _print_result(result, "go_enrichment")

@tool("kegg_enrichment", "KEGG pathway enrichment analysis")
def cmd_kegg_enrichment(args):
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True)
    p.add_argument("--organism", default="mmu")
    p.add_argument("--pvalue_cutoff", type=float, default=0.05)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    if os.path.isfile(a.genes):
        genes = pd.read_csv(a.genes, header=None)[0].tolist()
    else:
        genes = a.genes.split(",")
    result = go.kegg_enrich(genes, organism=a.organism, pvalue_cutoff=a.pvalue_cutoff)
    if isinstance(result, pd.DataFrame) and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        _print_result(result, "kegg_enrichment")


# ======================== NET ========================

@tool("module_identify", "Identify network modules using ClusterONE")
def cmd_module_identify(args):
    p = argparse.ArgumentParser()
    p.add_argument("--edges", required=True, help="Edge weight CSV file")
    p.add_argument("--module_size", type=int, default=5)
    p.add_argument("--output", default=None)
    a = p.parse_args(args)
    edges = pd.read_csv(a.edges)
    result = net.module_identify(edges, module_size=a.module_size)
    if isinstance(result, pd.DataFrame) and a.output:
        result.to_csv(a.output, index=False)
        print(f"Saved to {a.output}")
    else:
        _print_result(result, "module_identify")

@tool("hub_identify", "Identify hub genes from network")
def cmd_hub_identify(args):
    p = argparse.ArgumentParser()
    p.add_argument("--edges", required=True)
    p.add_argument("--clusters", required=True)
    a = p.parse_args(args)
    edges = pd.read_csv(a.edges)
    clusters = pd.read_csv(a.clusters)
    result = net.hub_identify(edges, clusters)
    _print_result(result, "hub_identify")


# ======================== MULTI ========================

@tool("get_optimal_threads", "Get optimal number of threads")
def cmd_get_optimal_threads(args):
    p = argparse.ArgumentParser()
    p.add_argument("--max", type=int, default=None)
    a = p.parse_args(args)
    result = multi.get_optimal_threads(a.max)
    print(f"Optimal threads: {result}")


# =============================================================================
# Main entry point
# =============================================================================

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    cmd = sys.argv[1]

    if cmd == "list":
        print("Available MRBIGR2 CLI tools:\n")
        categories = {
            "Genotype": ["snp_qc", "subset_genotype", "genotype_pca", "calculate_kinship", "calculate_ibd",
                         "snp_pruning", "snp_clumping", "snp_impute", "snp_stats",
                         "convert_vcf", "convert_hapmap", "plink_to_vcf", "build_tree"],
            "Phenotype": ["filter_abundance", "filter_missing", "scale_phenotype",
                          "impute_phenotype", "remove_outliers", "blup", "blue", "trait_correct"],
            "GWAS": ["gwas_lm", "gwas_lmm", "gwas_plink", "get_top_snps",
                     "calculate_lambda", "ensure_rs_id", "add_rs_id_to_vcf",
                     "generate_clump", "plink_clump"],
            "Visualization": ["plot_manhattan", "plot_qq", "plot_pca",
                              "plot_ld_heatmap", "plot_pheno_hist",
                              "plot_pheno_boxplot", "plot_pheno_correlation",
                              "gwas_summary"],
            "Annotation": ["parse_gtf", "annotate_snps", "get_genes_in_region",
                           "qtl_annotation"],
            "QTL": ["detect_qtl", "get_lead_snp", "identify_peak_snps",
                    "plot_qtl_region", "qtl_summary"],
            "MR": ["mr_analysis", "mr_causal_estimate", "format_qtl_for_mr"],
            "GO/KEGG": ["go_enrichment", "kegg_enrichment"],
            "Network": ["module_identify", "hub_identify"],
            "Utility": ["get_optimal_threads"],
        }
        for cat, names in categories.items():
            print(f"  {cat}:")
            for n in names:
                if n in TOOLS:
                    print(f"    {n:30s} {TOOLS[n]['desc']}")
            print()
        print(f"Total: {len(TOOLS)} tools")
        print(f"\nUsage: python {sys.argv[0]} <tool> [options]")
        print(f"Help:  python {sys.argv[0]} <tool> --help")
        sys.exit(0)

    if cmd == "help":
        if len(sys.argv) < 3:
            print("Usage: python mrbigr_cli.py help <tool>")
            sys.exit(1)
        tool_name = sys.argv[2]
        if tool_name not in TOOLS:
            print(f"Unknown tool: {tool_name}")
            sys.exit(1)
        TOOLS[tool_name]["fn"](["--help"])
        sys.exit(0)

    if cmd not in TOOLS:
        print(f"Unknown tool: {cmd}")
        print(f"Run 'python {sys.argv[0]} list' to see all available tools.")
        sys.exit(1)

    TOOLS[cmd]["fn"](sys.argv[2:])


if __name__ == "__main__":
    main()
