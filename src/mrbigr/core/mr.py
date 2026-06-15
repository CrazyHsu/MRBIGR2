#!/usr/bin/env python3
"""
MRMCP - Mendelian Randomization Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/mr.py

Functions:
- Basic MR analysis (IVW-like)
- MR with Linear Mixed Model (MLM)
- Target type analysis (QTL targeting)
"""
import pandas as pd
import numpy as np
from scipy.stats import chi2
from scipy.stats import norm
from sklearn.linear_model import LinearRegression
import os
import subprocess
import glob
import shutil
import warnings

warnings.filterwarnings("ignore")

# External tools
import os as mr_os
from .paths import repo_root as _repo_root
SCRIPT_DIR = str(_repo_root())
PLINK_BIN = os.path.join(SCRIPT_DIR, "utils", "plink")
GEMMA_BIN = os.path.join(SCRIPT_DIR, "utils", "gemma.linux")


def _ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def _coerce_gwas_df(data):
    from ._argjson import maybe_json_loads
    data = maybe_json_loads(data)
    if isinstance(data, pd.DataFrame):
        return data.copy()
    if isinstance(data, dict):
        return pd.DataFrame(data)
    return pd.DataFrame(data)


def _coerce_gene_df(data):
    from ._argjson import maybe_json_loads
    data = maybe_json_loads(data)
    if isinstance(data, list):
        return pd.DataFrame({'gene_id': data})

    df = data.copy() if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
    if 'geneid' in df.columns and 'gene_id' not in df.columns:
        df = df.rename(columns={'geneid': 'gene_id'})

    if 'position' in df.columns and not {'chr', 'start', 'end'}.issubset(df.columns):
        position = df['position'].astype(str).str.extract(r'(?P<chr>[^:]+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])')
        df = pd.concat([df, position], axis=1)
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')

    return df


# ========== Core MR Functions ==========

def lm_res(y, X):
    """Linear regression for SNP-effect calculation.
    
    Args:
        y: Outcome variable (array, shape n)
        X: Exposure variable (array, shape n or n x m)
    
    Returns:
        Series with effect, se, rsq (for 1D X) or DataFrame (for 2D X)
    """
    # Handle 2D array (multiple exposures)
    if X.ndim > 1:
        results = []
        for i in range(X.shape[1]):
            result = lm_res(y, X[:, i])
            results.append(result)
        return pd.DataFrame(results)
    
    # 1D case
    idx = ~(np.isnan(X) | np.isnan(y))
    X_clean = X[idx].reshape(-1, 1)
    y_clean = y[idx]
    
    if len(X_clean) < 3:
        return pd.Series({'effect': np.nan, 'se': np.nan, 'rsq': np.nan})
    
    lm = LinearRegression().fit(X_clean, y_clean)
    sigma2 = np.sum((y_clean - lm.predict(X_clean))**2) / (X_clean.shape[0] - 1)
    
    try:
        XtX_inv = np.linalg.pinv(np.dot(X_clean.T, X_clean))
        se = np.sqrt(np.diag(XtX_inv)[-1] * sigma2)
    except:
        se = np.nan
    
    effect = lm.coef_[-1] if len(lm.coef_) > 0 else np.nan
    
    try:
        rsq = lm.score(X_clean, y_clean)
    except:
        rsq = 0
    
    return pd.Series({'effect': effect, 'se': se, 'rsq': rsq})


def mr_analysis(exposure_beta, exposure_se=None, outcome_beta=None, outcome_se=None,
                method='ivw', snp_col='rs', beta_col='beta', se_col='se'):
    """Mendelian Randomization analysis.
    
    Args:
        exposure_beta: Beta coefficients for exposure, or exposure GWAS DataFrame
        exposure_se: Standard errors for exposure, or outcome GWAS DataFrame
        outcome_beta: Beta coefficients for outcome
        outcome_se: Standard errors for outcome
        method: MR method ('ivw', 'egger')
    
    Returns:
        DataFrame with MR results
    """
    # Allow direct GWAS DataFrame input for CLI/MCP wrappers.
    if isinstance(exposure_beta, pd.DataFrame) and isinstance(exposure_se, pd.DataFrame) \
            and outcome_beta is None and outcome_se is None:
        merged = pd.merge(
            exposure_beta[[snp_col, beta_col, se_col]],
            exposure_se[[snp_col, beta_col, se_col]],
            on=snp_col,
            suffixes=('_exp', '_out')
        )
        if len(merged) < 3:
            return pd.DataFrame()
        result = mr_analysis(
            merged[f'{beta_col}_exp'].values,
            merged[f'{se_col}_exp'].values,
            merged[f'{beta_col}_out'].values,
            merged[f'{se_col}_out'].values,
            method=method
        )
        result['n_snps'] = len(merged)
        return result

    exposure_beta = np.asarray(exposure_beta, dtype=float)
    exposure_se = np.asarray(exposure_se, dtype=float)
    outcome_beta = np.asarray(outcome_beta, dtype=float)
    outcome_se = np.asarray(outcome_se, dtype=float)

    idx = ~(np.isnan(exposure_beta) | np.isnan(exposure_se) | np.isnan(outcome_beta) | np.isnan(outcome_se))
    exposure_beta = exposure_beta[idx]
    exposure_se = exposure_se[idx]
    outcome_beta = outcome_beta[idx]
    outcome_se = outcome_se[idx]

    if len(exposure_beta) < 1:
        return pd.DataFrame()

    if str(method).lower() in {'weighted_median', 'weighted-median', 'wm'}:
        valid = (exposure_beta != 0) & (outcome_se > 0) & (exposure_se >= 0)
        if valid.sum() < 1:
            return pd.DataFrame()
        bx = exposure_beta[valid]
        bx_se = exposure_se[valid]
        by = outcome_beta[valid]
        by_se = outcome_se[valid]
        ratios = by / bx
        ratio_var = (by_se ** 2 / bx ** 2) + ((by ** 2) * (bx_se ** 2) / (bx ** 4))
        ratio_var = np.where(ratio_var <= 0, np.nan, ratio_var)
        valid_ratio = ~np.isnan(ratios) & ~np.isnan(ratio_var)
        if valid_ratio.sum() < 1:
            return pd.DataFrame()
        ratios = ratios[valid_ratio]
        weights = 1 / ratio_var[valid_ratio]
        order = np.argsort(ratios)
        ratios = ratios[order]
        weights = weights[order]
        cum_weights = np.cumsum(weights) / np.sum(weights)
        wm_effect = ratios[np.searchsorted(cum_weights, 0.5)]
        wm_se = np.sqrt(1 / np.sum(weights))
        TMR = (wm_effect / wm_se) ** 2 if wm_se > 0 else np.nan
        pvalue = 1 - chi2.cdf(TMR, 1) if not np.isnan(TMR) else np.nan
        return pd.DataFrame({
            'method': ['Weighted median'],
            'effect': [wm_effect],
            'se': [wm_se],
            'TMR': [TMR],
            'pvalue': [pvalue]
        })

    # Calculate IVW estimate
    # Effect = sum(bx*by/sx^2) / sum(bx^2/sx^2)
    weights = 1 / (exposure_se ** 2)
    
    # IVW method
    ivw_effect = np.sum(exposure_beta * outcome_beta * weights) / np.sum(exposure_beta ** 2 * weights)
    
    # Variance and SE
    ivw_var = 1 / np.sum(exposure_beta ** 2 * weights)
    ivw_se = np.sqrt(ivw_var)
    
    # Wald statistic
    TMR = (ivw_effect ** 2) / ivw_var
    pvalue = 1 - chi2.cdf(TMR, 1)
    
    # MR-Egger (if enough SNPs)
    if method == 'egger' and len(exposure_beta) > 2:
        # Simple Egger: regress outcome_beta on exposure_beta
        idx = ~(np.isnan(exposure_beta) | np.isnan(outcome_beta))
        if idx.sum() > 2:
            egger_model = LinearRegression().fit(
                exposure_beta[idx].reshape(-1, 1), 
                outcome_beta[idx]
            )
            egger_intercept = egger_model.intercept_
            egger_slope = egger_model.coef_[0]
            
            # Calculate Egger SE
            residuals = outcome_beta[idx] - egger_intercept - egger_slope * exposure_beta[idx]
            egger_var = np.sum(residuals**2) / (idx.sum() - 2)
            egger_se = np.sqrt(egger_var / np.sum((exposure_beta[idx] - exposure_beta[idx].mean())**2))
            
            return pd.DataFrame({
                'method': ['IVW', 'MR-Egger'],
                'effect': [ivw_effect, egger_slope],
                'se': [ivw_se, egger_se],
                'intercept': [np.nan, egger_intercept],
                'TMR': [TMR, egger_slope**2 / egger_se**2],
                'pvalue': [pvalue, 1 - chi2.cdf(egger_slope**2 / egger_se**2, 1)]
            })
    
    return pd.DataFrame({
        'method': ['IVW'],
        'effect': [ivw_effect],
        'se': [ivw_se],
        'TMR': [TMR],
        'pvalue': [pvalue]
    })


def mr_causal_estimate(exposure_gwas, outcome_gwas, snp_col='rs', beta_col='beta', se_col='se', method='ivw'):
    """Calculate MR causal estimate between two traits.
    
    Args:
        exposure_gwas: GWAS summary for exposure trait
        outcome_gwas: GWAS summary for outcome trait
        snp_col: SNP column name
        beta_col: Beta column name
        se_col: SE column name
    
    Returns:
        DataFrame with MR results
    """
    # Merge on SNP
    merged = pd.merge(
        exposure_gwas[[snp_col, beta_col, se_col]],
        outcome_gwas[[snp_col, beta_col, se_col]],
        on=snp_col,
        suffixes=('_exp', '_out')
    )
    
    if len(merged) < 3:
        return pd.DataFrame()
    
    # Run MR analysis
    result = mr_analysis(
        merged[f'{beta_col}_exp'].values,
        merged[f'{se_col}_exp'].values,
        merged[f'{beta_col}_out'].values,
        merged[f'{se_col}_out'].values,
        method=method
    )
    
    result['n_snps'] = len(merged)
    
    return result


# ========== MR with MLM (GEMMA) ==========

def prepare_mr_geno(exposure_qtl, outcome_qtl, genotype_prefix, output_dir):
    """Prepare genotype files for MR analysis.
    
    Args:
        exposure_qtl: QTL for exposure trait
        outcome_qtl: QTL for outcome trait
        genotype_prefix: PLINK genotype prefix
        output_dir: Output directory
    
    Returns:
        Dictionary with output paths
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Get all unique SNPs
    all_snps = set(exposure_qtl['SNP'].unique()) | set(outcome_qtl['SNP'].unique())
    
    # Write SNP list
    snp_file = f"{output_dir}/mr_snps.txt"
    with open(snp_file, 'w') as f:
        for snp in all_snps:
            f.write(f"{snp}\n")
    
    # Extract SNPs using PLINK
    output_prefix = f"{output_dir}/mr_geno"
    cmd = f"{PLINK_BIN} --bfile {genotype_prefix} --extract {snp_file} --make-bed --out {output_prefix} --allow-extra-chr"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    return {
        'bed': output_prefix,
        'snps': list(all_snps)
    }


def run_mr_gemma(exposure_qtl, outcome_qtl, genotype_prefix, output_dir, threads=1):
    """Run MR analysis using GEMMA for effect estimation.
    
    Args:
        exposure_qtl: QTL DataFrame for exposure
        outcome_qtl: QTL DataFrame for outcome
        genotype_prefix: PLINK genotype prefix
        output_dir: Output directory
        threads: Number of threads
    
    Returns:
        MR results DataFrame
    """
    os.makedirs(output_dir, exist_ok=True)
    gemma_output_dir = os.path.join(output_dir, 'output')
    os.makedirs(gemma_output_dir, exist_ok=True)
    
    # Prepare genotype
    geno_info = prepare_mr_geno(exposure_qtl, outcome_qtl, genotype_prefix, output_dir)
    if geno_info is None:
        return None
    
    geno_prefix = geno_info['bed']
    geno_name = os.path.basename(geno_prefix)
    
    # Calculate relatedness matrix
    kinship_cmd = f"{GEMMA_BIN} -bfile {geno_prefix} -gk 1 -o {geno_name}"
    subprocess.run(kinship_cmd, shell=True, capture_output=True, cwd=output_dir)
    
    kinship_file = os.path.join(gemma_output_dir, f"{geno_name}.cXX.txt")
    
    # Run GWAS for each exposure and outcome
    results = []
    
    # Exposure traits
    for exp_trait in exposure_qtl['phe_name'].unique():
        exp_snps = exposure_qtl[exposure_qtl['phe_name'] == exp_trait]['SNP']
        
        # Create phenotype file for this exposure
        fam_file = geno_prefix + '.fam'
        fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                         names=['fam', 'id', 'pat', 'mat', 'sex', 'pheno'])
        fam[5] = 1
        fam.to_csv(fam_file, sep=' ', header=False, index=False)
        
        # Run GEMMA
        cmd = f"{GEMMA_BIN} -bfile {geno_prefix} -k {kinship_file} -lmm -n 1 -o {exp_trait}_exp"
        subprocess.run(cmd, shell=True, capture_output=True, cwd=output_dir)
        
        # Read results
        assoc_file = os.path.join(gemma_output_dir, f"{exp_trait}_exp.assoc.txt")
        if os.path.exists(assoc_file):
            exp_effect = pd.read_csv(assoc_file, sep='\t')
            exp_effect.index = exp_effect['rs']
            exp_effect = exp_effect[['beta', 'se']]
        else:
            continue
        
        # Outcome traits
        for out_trait in outcome_qtl['phe_name'].unique():
            out_snps = outcome_qtl[outcome_qtl['phe_name'] == out_trait]['SNP']
            
            # Update phenotype
            # (Simplified - would need proper phenotype matching)
            
            cmd = f"{GEMMA_BIN} -bfile {geno_prefix} -k {kinship_file} -lmm -n 2 -o {out_trait}_out"
            subprocess.run(cmd, shell=True, capture_output=True, cwd=output_dir)
            
            out_file = os.path.join(gemma_output_dir, f"{out_trait}_out.assoc.txt")
            if os.path.exists(out_file):
                out_effect = pd.read_csv(out_file, sep='\t')
                out_effect.index = out_effect['rs']
                out_effect = out_effect[['beta', 'se']]
            else:
                continue
            
            # Calculate MR for overlapping SNPs
            common_snps = list(set(exp_snps) & set(out_snps) & set(exp_effect.index) & set(out_effect.index))
            
            if len(common_snps) < 3:
                continue
            
            # MR calculation
            exp_b = exp_effect.loc[common_snps, 'beta']
            exp_s = exp_effect.loc[common_snps, 'se']
            out_b = out_effect.loc[common_snps, 'beta']
            out_s = out_effect.loc[common_snps, 'se']
            
            mr_result = mr_analysis(exp_b.values, exp_s.values, out_b.values, out_s.values)
            mr_result['exposure'] = exp_trait
            mr_result['outcome'] = out_trait
            mr_result['n_snps'] = len(common_snps)
            
            results.append(mr_result)
    
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()


# ========== QTL Targeting Analysis ==========

def qtl_target_analysis(qtl_df, tf_genes, target_genes, window=500000):
    """QTL targeting analysis - find which TFs target which genes.
    
    Args:
        qtl_df: QTL DataFrame with CHR, SNP, phe_name columns
        tf_genes: DataFrame with TF gene IDs
        target_genes: DataFrame with target gene IDs
    
    Returns:
        DataFrame with targeting relationships
    """
    # Try to import pyranges
    try:
        import pyranges as pr
        HAS_PYRANGES = True
    except ImportError:
        HAS_PYRANGES = False
    
    qtl_df = qtl_df.copy() if isinstance(qtl_df, pd.DataFrame) else pd.DataFrame(qtl_df)
    tf_genes = _coerce_gene_df(tf_genes)
    target_genes = _coerce_gene_df(target_genes)

    target_col = 'gene_id' if 'gene_id' in target_genes.columns else target_genes.columns[0]
    target_ids = set(target_genes[target_col].dropna().astype(str)) if len(target_genes.columns) > 0 else set()
    if target_ids and 'phe_name' in qtl_df.columns:
        qtl_df = qtl_df[qtl_df['phe_name'].astype(str).isin(target_ids)]

    if not HAS_PYRANGES:
        results = []
        for _, qtl in qtl_df.iterrows():
            qtl_pos = int(qtl['SNP'].split('_')[-1]) if '_' in qtl['SNP'] else 0
            qtl_chr = str(qtl['CHR'])

            for _, tf in tf_genes.iterrows():
                if 'start' in tf.index and 'end' in tf.index and pd.notna(tf['start']) and pd.notna(tf['end']):
                    if (str(tf.get('chr', '')) == qtl_chr and 
                        tf['start'] - window <= qtl_pos <= tf['end'] + window):
                        results.append({
                            'tf_gene': tf.get('gene_id', tf.get('gene_name', 'unknown')),
                            'target_gene': qtl['phe_name'],
                            'snp': qtl['SNP'],
                            'qtl_chr': qtl_chr,
                            'qtl_pos': qtl_pos
                        })
        return pd.DataFrame(results)
    
    # Use pyranges for more accurate overlap
    tf_ranges = tf_genes[['chr', 'start', 'end', 'gene_id']].dropna().copy()
    tf_ranges.columns = ['Chromosome', 'Start', 'End', 'gene_id']
    tf_pr = pr.PyRanges(tf_ranges)
    
    # Convert QTL to ranges
    qtl_ranges = qtl_df[['CHR', 'SNP', 'phe_name']].copy()
    qtl_ranges['Start'] = qtl_ranges['SNP'].apply(lambda x: int(x.split('_')[-1]) - window if '_' in str(x) else 0)
    qtl_ranges['End'] = qtl_ranges['SNP'].apply(lambda x: int(x.split('_')[-1]) + window if '_' in str(x) else 0)
    qtl_ranges['phe_name'] = qtl_df['phe_name']
    qtl_ranges.columns = ['Chromosome', 'SNP', 'phe_name', 'Start', 'End']
    qtl_pr = pr.PyRanges(qtl_ranges[['Chromosome', 'Start', 'End', 'phe_name']])
    
    # Join
    overlap = tf_pr.join(qtl_pr)
    
    results = []
    for k in sorted(overlap.dfs.keys()):
        df = overlap.dfs[k]
        for _, row in df.iterrows():
            results.append({
                'tf_gene': row.get('gene_id', 'unknown'),
                'target_gene': row.get('phe_name', 'unknown')
            })
    
    return pd.DataFrame(results)


# ========== pleiotropy and heterogeneity ==========

def test_pleiotropy(exposure_gwas, outcome_gwas):
    """Test for horizontal pleiotropy (MR-Egger intercept).
    
    Args:
        exposure_gwas: Exposure GWAS summary
        outcome_gwas: Outcome GWAS summary
    
    Returns:
        Dictionary with test results
    """
    merged = pd.merge(
        exposure_gwas[['rs', 'beta', 'se']],
        outcome_gwas[['rs', 'beta', 'se']],
        on='rs',
        suffixes=('_exp', '_out')
    )
    
    if len(merged) < 3:
        return {'n_snps': 0, 'intercept': np.nan, 'pvalue': np.nan}
    
    # MR-Egger
    idx = ~(merged['beta_exp'].isna() | merged['beta_out'].isna())
    X = merged.loc[idx, 'beta_exp'].values.reshape(-1, 1)
    y = merged.loc[idx, 'beta_out'].values
    
    model = LinearRegression().fit(X, y)
    intercept = model.intercept_
    slope = model.coef_[0]
    
    # Calculate SE of intercept
    residuals = y - intercept - slope * X.flatten()
    n = len(y)
    s_sq = np.sum(residuals**2) / (n - 2)
    mean_x = X.mean()
    intercept_se = np.sqrt(s_sq * (1/n + mean_x**2 / np.sum((X - mean_x)**2)))
    
    # Z-test for intercept
    z = intercept / intercept_se
    pvalue = 2 * norm.cdf(-abs(z))
    
    return {
        'n_snps': len(merged),
        'intercept': intercept,
        'intercept_se': intercept_se,
        'pvalue': pvalue
    }


def heterogeneity_test(exposure_gwas, outcome_gwas):
    """Test for heterogeneity using IVW method.
    
    Args:
        exposure_gwas: Exposure GWAS summary
        outcome_gwas: Outcome GWAS summary
    
    Returns:
        Dictionary with heterogeneity statistics
    """
    merged = pd.merge(
        exposure_gwas[['rs', 'beta', 'se']],
        outcome_gwas[['rs', 'beta', 'se']],
        on='rs',
        suffixes=('_exp', '_out')
    )
    
    if len(merged) < 3:
        return {'q_statistic': np.nan, 'q_df': 0, 'q_pvalue': np.nan}
    
    # IVW estimate
    weights = 1 / (merged['se_exp'] ** 2)
    ivw_effect = np.sum(merged['beta_exp'] * merged['beta_out'] * weights) / np.sum(merged['beta_exp'] ** 2 * weights)
    
    # Q statistic
    Q = np.sum(weights * (merged['beta_out'] - ivw_effect * merged['beta_exp']) ** 2)
    Q_df = len(merged) - 1
    Q_pvalue = 1 - chi2.cdf(Q, Q_df)
    
    return {
        'q_statistic': Q,
        'q_df': Q_df,
        'q_pvalue': Q_pvalue
    }


# ========== Utility Functions ==========

def read_gwas_summary(gwas_file):
    """Read GWAS summary statistics.
    
    Args:
        gwas_file: GWAS result file (GEMMA format)
    
    Returns:
        DataFrame with rs, beta, se, pvalue
    """
    df = pd.read_csv(gwas_file, sep='\t')
    return df[['rs', 'beta', 'se', 'p_wald']].rename(columns={'p_wald': 'pvalue'})


def _infer_qtl_trait_name(qtl_file):
    """Infer a single-trait name from a QTL file path when no ``phe_name``
    column is present. GWAS wrappers commonly write
    ``<output_name>_<trait>.qtl.csv``; use the final underscore token as the
    trait label for downstream MR/QTL-target tools."""
    if isinstance(qtl_file, (str, os.PathLike)):
        name = os.path.basename(str(qtl_file))
        if name.endswith(".csv"):
            name = name[:-4]
        if name.endswith(".qtl"):
            name = name[:-4]
        if "_" in name:
            return name.rsplit("_", 1)[-1]
        return name
    return "trait"


def format_qtl_for_mr(qtl_file, output_file=None):
    """Format QTL file for MR analysis.
    
    Args:
        qtl_file: QTL result CSV
    
    Returns:
        DataFrame with required columns
    """
    df = qtl_file.copy() if isinstance(qtl_file, pd.DataFrame) else pd.read_csv(qtl_file)
    if 'phe_name' not in df.columns:
        df['phe_name'] = _infer_qtl_trait_name(qtl_file)
    # Ensure required columns exist
    required = ['CHR', 'SNP', 'phe_name', 'P']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"QTL file missing required column: {col}")
    if output_file is not None:
        _ensure_parent_dir(output_file)
        df.to_csv(output_file, index=False)
    return df


def heterogeneity_test_simple(betas, ses):
    """Simple heterogeneity test with arrays.
    
    Args:
        betas: Array of beta coefficients
        ses: Array of standard errors
    
    Returns:
        Dictionary with heterogeneity statistics
    """
    betas = np.array(betas)
    ses = np.array(ses)
    
    if len(betas) < 3 or len(ses) < 3:
        return {'q_statistic': np.nan, 'q_df': 0, 'q_pvalue': np.nan}
    
    weights = 1 / (ses ** 2)
    ivw_effect = np.sum(betas * weights) / np.sum(weights)
    Q = np.sum(weights * (betas - ivw_effect) ** 2)
    q_df = len(betas) - 1
    q_pvalue = 1 - chi2.cdf(Q, q_df)
    
    return {'q_statistic': Q, 'q_df': q_df, 'q_pvalue': q_pvalue}


def test_pleiotropy_simple(exposure_beta, outcome_beta, geno):
    """Simple pleiotropy test with arrays.
    
    Args:
        exposure_beta: Exposure effect sizes
        outcome_beta: Outcome effect sizes
        geno: Genotype matrix
    
    Returns:
        Dictionary with pleiotropy test results
    """
    from sklearn.linear_model import LinearRegression
    
    if len(exposure_beta) != len(outcome_beta):
        return {'n_snps': 0, 'intercept': np.nan, 'pvalue': np.nan}
    
    idx = ~(np.isnan(exposure_beta) | np.isnan(outcome_beta))
    if idx.sum() < 3:
        return {'n_snps': idx.sum(), 'intercept': np.nan, 'pvalue': np.nan}
    
    X = exposure_beta[idx].reshape(-1, 1)
    y = outcome_beta[idx]
    
    model = LinearRegression().fit(X, y)
    intercept = model.intercept_
    slope = model.coef_[0]
    
    residuals = y - intercept - slope * X.flatten()
    var = np.sum(residuals ** 2) / (idx.sum() - 2)
    se = np.sqrt(var / np.sum((X - X.mean()) ** 2))
    
    TMR = intercept ** 2 / (se ** 2) if se > 0 else 0
    pvalue = 1 - chi2.cdf(TMR, 1)
    
    return {'n_snps': idx.sum(), 'intercept': intercept, 'pvalue': pvalue, 'slope': slope}


__all__ = [
    'lm_res', 'mr_analysis', 'mr_causal_estimate',
    'prepare_mr_geno', 'run_mr_gemma',
    'qtl_target_analysis',
    'test_pleiotropy', 'heterogeneity_test',
    'read_gwas_summary', 'format_qtl_for_mr'
]
