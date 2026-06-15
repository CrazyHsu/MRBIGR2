#!/usr/bin/env python3
"""
GWASMCP - GWAS Analysis Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/gwas.py
"""
import pandas as pd
import numpy as np
import os
import subprocess
import glob
import shutil
from scipy import stats
from scipy.stats import chi2
import statsmodels.api as sm
import warnings
from . import parallel as multi
from .paths import repo_root as _repo_root

warnings.filterwarnings("ignore")


# GEMMA binary path
import os as gwas_os
SCRIPT_DIR = str(_repo_root())
GEMMA_BIN = os.path.join(SCRIPT_DIR, "utils", "gemma.linux")
PLINK_BIN = os.path.join(SCRIPT_DIR, "utils", "plink")


def _resolve_num_threads(num_threads):
    """Use multi.py as the single thread-count policy source."""
    return multi.get_optimal_threads(num_threads)


def _run_shell_commands(cmds, num_threads, label):
    """Run shell commands sequentially or via multi.parallel_run."""
    n_parallel = min(_resolve_num_threads(num_threads), len(cmds))
    if n_parallel <= 1 or len(cmds) <= 1:
        return [subprocess.run(cmd, shell=True, capture_output=True, text=True).returncode for cmd in cmds]

    print(f"[info] running {len(cmds)} {label} jobs in parallel (max {n_parallel} concurrent)")
    return multi.parallel_run(multi.run_cmd, [(cmd,) for cmd in cmds], num_threads=n_parallel)


def read_plink_prefix(geno_prefix):
    """Read PLINK genotype files."""
    try:
        bed_file = geno_prefix + '.bed'
        bim_file = geno_prefix + '.bim'
        fam_file = geno_prefix + '.fam'
        
        # Read bim file
        snps = pd.read_csv(bim_file, sep="\t", header=None, 
                           names=["chr", "snp", "cm", "pos", "a1", "a2"])
        
        # Read fam file
        fam = pd.read_csv(fam_file, sep=r"\s+", header=None,
                          names=["fam", "id", "pat", "mat", "sex", "pheno"])
        
        return {'snps': snps, 'fam': fam}
    except Exception as e:
        return None


def prepare_gemma_fam(geno_prefix, phe, output_fam_path, pheno_col=None):
    """Create a new .fam file with phenotype data for GEMMA.

    Reads the original FAM (never modifies it) and writes a new FAM to
    *output_fam_path* with phenotype values merged from *phe*.
    Supports multiple phenotype columns for GEMMA's -n flag.

    Args:
        geno_prefix: PLINK genotype prefix (original, read-only)
        phe: Phenotype DataFrame with sample IDs as index (one or more columns)
        output_fam_path: Where to write the new FAM file
        pheno_col: Column name in phe to use; defaults to all columns

    Returns:
        Path to the newly created FAM file
    """
    fam_file = geno_prefix + '.fam'
    fam = pd.read_csv(fam_file, sep=r"\s+", header=None,
                      names=["fam", "id", "pat", "mat", "sex", "pheno"])

    fam['id'] = fam['id'].astype(str)
    phe = phe.copy()
    phe.index = phe.index.astype(str)

    # Drop the default pheno column — will be replaced
    fam = fam.drop(columns=['pheno'])

    if pheno_col is not None:
        pheno_cols = [pheno_col] if isinstance(pheno_col, str) else list(pheno_col)
    else:
        pheno_cols = list(phe.columns)

    for col in pheno_cols:
        fam[col] = fam['id'].map(phe[col]).fillna(-9)

    matched = (fam[pheno_cols[0]] != -9).sum()
    print(f"[info] prepare_gemma_fam: {matched}/{len(fam)} samples matched, {len(pheno_cols)} phenotype(s)")

    fam.to_csv(output_fam_path, sep='\t', header=False, index=False)
    return output_fam_path


def prepare_gemma_covariates(geno_prefix, cov, output_cov_path):
    """Create a GEMMA covariate file aligned to the PLINK FAM sample order."""
    if isinstance(cov, str):
        cov_df = pd.read_csv(cov, sep=None, engine='python')
    else:
        cov_df = pd.DataFrame(cov)

    sample_col = None
    for col in cov_df.columns:
        col_text = str(col).strip()
        if col_text.lower() in {'id', 'iid', 'sample', 'sample_id', 'sampleid', 'genotype', 'accession'}:
            sample_col = col
            break
    if sample_col is not None:
        cov_df = cov_df.set_index(sample_col)

    cov_df.index = cov_df.index.astype(str)
    cov_df = cov_df.apply(pd.to_numeric, errors='coerce')

    fam_file = geno_prefix + '.fam'
    fam = pd.read_csv(fam_file, sep=r"\s+", header=None,
                      names=["fam", "id", "pat", "mat", "sex", "pheno"])
    sample_ids = fam['id'].astype(str)
    missing = sample_ids[~sample_ids.isin(cov_df.index)]
    if len(missing) > 0:
        raise ValueError(f"covariates missing {len(missing)} samples; first missing: {missing.iloc[0]}")

    aligned = cov_df.loc[sample_ids]
    if aligned.isna().any().any():
        raise ValueError("covariates contain missing or non-numeric values after sample alignment")
    aligned.to_csv(output_cov_path, sep='\t', header=False, index=False)
    return output_cov_path


def gwas_lm(phe, geno_prefix, output_name=None, output_dir=None, num_threads=None):
    """Linear Model GWAS using GEMMA.

    Args:
        phe: Phenotype DataFrame (samples x phenotypes)
        geno_prefix: PLINK genotype prefix
        output_name: Output name for results
        output_dir: Directory for output files (default: ./output)
        num_threads: Number of threads

    Returns:
        List of output files
    """
    num_threads = _resolve_num_threads(num_threads)

    if output_dir is None:
        output_dir = 'output'
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    geno_name = os.path.basename(geno_prefix)
    link_prefix = os.path.join(output_dir, f"{geno_name}.link")

    bed_file = geno_prefix + '.bed'
    bim_file = geno_prefix + '.bim'
    fam_file = geno_prefix + '.fam'

    if not all(os.path.exists(f) for f in [bed_file, bim_file, fam_file]):
        return None

    for ext in ['.bed', '.bim']:
        link = link_prefix + ext
        if os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(geno_prefix + ext), link)

    link_fam = link_prefix + '.fam'
    prepare_gemma_fam(geno_prefix, phe, link_fam)

    cmds = []
    out_names = []
    for i, pheno_name in enumerate(phe.columns):
        pheno_name_safe = str(pheno_name).replace('/', '.').replace(' ', '_')
        if output_name is None:
            out_name = pheno_name_safe
        else:
            out_name = f"{output_name}_{pheno_name_safe}"
        out_names.append(out_name)
        cmds.append(f"{GEMMA_BIN} -bfile {link_prefix} -lm -n {i+1} -outdir {output_dir} -o {out_name}")

    procs_results = _run_shell_commands(cmds, num_threads, "GWAS")

    outputs = []
    for i, out_name in enumerate(out_names):
        if procs_results[i] == 0:
            outputs.append(os.path.join(output_dir, f"{out_name}.assoc.txt"))

    for ext in ['.bed', '.bim', '.fam']:
        if os.path.exists(link_prefix + ext):
            os.remove(link_prefix + ext)

    return outputs


def gwas_lmm(phe, geno_prefix, output_name=None, output_dir=None, num_threads=None, cov=None, kinship_file=None):
    """Linear Mixed Model GWAS using GEMMA.

    Args:
        phe: Phenotype DataFrame (samples x phenotypes)
        geno_prefix: PLINK genotype prefix
        output_name: Output name for results
        output_dir: Directory for output files (default: ./output)
        num_threads: Number of threads
        kinship_file: Optional precomputed GEMMA kinship matrix. When supplied,
            GWAS reuses this file instead of generating ``<output_dir>/<geno>.cXX.txt``.

    Returns:
        List of output files
    """
    num_threads = _resolve_num_threads(num_threads)

    if output_dir is None:
        output_dir = 'output'
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    geno_name = os.path.basename(geno_prefix)
    link_prefix = os.path.join(output_dir, f"{geno_name}.link")

    bed_file = geno_prefix + '.bed'
    bim_file = geno_prefix + '.bim'
    fam_file = geno_prefix + '.fam'

    if not all(os.path.exists(f) for f in [bed_file, bim_file, fam_file]):
        return None

    for ext in ['.bed', '.bim']:
        link = link_prefix + ext
        if os.path.exists(link):
            os.remove(link)
        os.symlink(os.path.abspath(geno_prefix + ext), link)

    link_fam = link_prefix + '.fam'
    prepare_gemma_fam(geno_prefix, phe, link_fam)
    cov_arg = ""
    if cov is not None:
        cov_file = link_prefix + '.covariates.txt'
        prepare_gemma_covariates(geno_prefix, cov, cov_file)
        cov_arg = f" -c {cov_file}"

    if kinship_file is not None:
        kinship_file = os.path.abspath(str(kinship_file))
        if not os.path.isfile(kinship_file):
            raise FileNotFoundError(f"kinship file not found: {kinship_file}")
    else:
        kinship_file = os.path.join(output_dir, f"{geno_name}.cXX.txt")
    if not os.path.exists(kinship_file):
        cmd_kinship = f"{GEMMA_BIN} -bfile {link_prefix} -gk 1 -outdir {output_dir} -o {geno_name}"
        subprocess.run(cmd_kinship, shell=True, capture_output=True)

    cmds = []
    out_names = []
    for i, pheno_name in enumerate(phe.columns):
        pheno_name_safe = str(pheno_name).replace('/', '.').replace(' ', '_')
        if output_name is None:
            out_name = pheno_name_safe
        else:
            out_name = f"{output_name}_{pheno_name_safe}"
        out_names.append(out_name)
        cmds.append(f"{GEMMA_BIN} -bfile {link_prefix} -k {kinship_file} -lmm -n {i+1}{cov_arg} -outdir {output_dir} -o {out_name}")

    procs_results = _run_shell_commands(cmds, num_threads, "GWAS")

    outputs = []
    for i, out_name in enumerate(out_names):
        if procs_results[i] == 0:
            outputs.append(os.path.join(output_dir, f"{out_name}.assoc.txt"))

    for ext in ['.bed', '.bim', '.fam', '.covariates.txt']:
        if os.path.exists(link_prefix + ext):
            os.remove(link_prefix + ext)

    return outputs


def gwas_plink(phe, geno_prefix, pheno_col=None, output_name="gwas"):
    """Run GWAS using PLINK linear regression.
    
    Args:
        phe: Phenotype DataFrame
        geno_prefix: PLINK genotype prefix
        pheno_col: Phenotype column name
        output_name: Output file name
    
    Returns:
        Output file path
    """
    # Create temporary phenotype file
    temp_phe_file = f"{output_name}_temp.phe"
    
    # Prepare phenotype file
    fam = pd.read_csv(geno_prefix + '.fam', sep=r"\s+", header=None,
                      names=["fam", "id", "pat", "mat", "sex", "pheno"])
    
    if pheno_col is not None and pheno_col in phe.columns:
        phe_df = phe[[pheno_col]].reset_index()
        phe_df.columns = ['id', 'pheno']
    else:
        first_phe = phe.columns[0]
        phe_df = phe[[first_phe]].reset_index()
        phe_df.columns = ['id', 'pheno']
    fam['id'] = fam['id'].astype(str)
    fam['fam'] = fam['fam'].astype(str)
    phe_df['id'] = phe_df['id'].astype(str)
    merged = pd.merge(fam[['fam', 'id']], phe_df, on='id', how='left')
    merged['pheno'] = merged['pheno'].fillna(-9)

    merged.to_csv(temp_phe_file, sep=' ', index=False, header=False)

    # Run PLINK GWAS
    cmd = f"{PLINK_BIN} --bfile {geno_prefix} --pheno {temp_phe_file} --linear --out {output_name} --allow-no-sex"

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Cleanup
    if os.path.exists(temp_phe_file):
        os.remove(temp_phe_file)

    if result.returncode == 0:
        return f"{output_name}.assoc.linear"
    return None


def simple_gwas(Y, G, cov=None):
    """Simple linear regression GWAS (pure Python, no external tools).
    
    Args:
        Y: Phenotype vector (n_samples,)
        G: Genotype matrix (n_snps x n_samples)
        cov: Covariate matrix (n_samples x n_cov)
    
    Returns:
        DataFrame with SNP, beta, se, pvalue
    """
    n_snps = G.shape[0]
    results = []
    
    # Add intercept to covariates
    if cov is not None:
        X = np.column_stack([np.ones(len(Y)), cov])
    else:
        X = np.ones((len(Y), 1))
    
    # Remove missing phenotypes
    valid_idx = ~np.isnan(Y)
    Y_valid = Y[valid_idx]
    X_valid = X[valid_idx]
    G_valid = G[:, valid_idx]
    
    # Fit null model
    XtX = np.dot(X_valid.T, X_valid)
    XtX_inv = np.linalg.pinv(XtX)
    beta_null = np.dot(XtX_inv, np.dot(X_valid.T, Y_valid))
    residuals = Y_valid - np.dot(X_valid, beta_null)
    sigma2 = np.var(residuals)
    
    # Test each SNP
    for snp_idx in range(n_snps):
        snp = G_valid[snp_idx]
        snp_valid = ~np.isnan(snp)
        
        if snp_valid.sum() < 10:  # Skip if too few samples
            results.append({'snp': snp_idx, 'beta': np.nan, 'se': np.nan, 'pvalue': np.nan})
            continue
        
        # Combine SNP with covariates
        X_snp = np.column_stack([X_valid, snp])
        
        try:
            XtX_snp = np.dot(X_snp.T, X_snp)
            XtX_snp_inv = np.linalg.pinv(XtX_snp)
            beta = np.dot(XtX_snp_inv, np.dot(X_snp.T, Y_valid))
            
            # Calculate residuals and SE
            residuals_snp = Y_valid - np.dot(X_snp, beta)
            df = len(Y_valid) - X_snp.shape[1]
            sigma2_snp = np.sum(residuals_snp**2) / df
            
            se = np.sqrt(np.diag(XtX_snp_inv) * sigma2_snp)
            beta_snp = beta[-1]
            se_snp = se[-1]
            
            # Wald test
            if se_snp > 0:
                z = beta_snp / se_snp
                pvalue = 2 * stats.norm.cdf(-abs(z))
            else:
                pvalue = np.nan
            
            results.append({'snp': snp_idx, 'beta': beta_snp, 'se': se_snp, 'pvalue': pvalue})
        except:
            results.append({'snp': snp_idx, 'beta': np.nan, 'se': np.nan, 'pvalue': np.nan})
    
    return pd.DataFrame(results)


def generate_clump_input(gwas_dir):
    """Generate clump input files from GWAS results.
    
    Args:
        gwas_dir: Directory containing GWAS association files
    
    Returns:
        Input directory path
    """
    input_dir = os.path.join(os.path.abspath(gwas_dir), 'clump_input')
    if os.path.exists(input_dir):
        shutil.rmtree(input_dir)
    os.makedirs(input_dir)
    
    for fn in glob.glob(os.path.join(gwas_dir.rstrip('/'), '*.assoc.txt')):
        filename = os.path.basename(fn)
        assoc = pd.read_csv(fn, sep='\t')
        
        # Extract SNP and p-value columns
        if 'rs' in assoc.columns and 'p_wald' in assoc.columns:
            assoc = assoc[['rs', 'p_wald']]
            assoc.columns = ['SNP', 'P']
        elif 'snp' in assoc.columns and 'p_score' in assoc.columns:
            assoc = assoc[['snp', 'p_score']]
            assoc.columns = ['SNP', 'P']
        
        out_file = os.path.join(input_dir, filename.replace('.assoc.txt', '.assoc'))
        assoc.to_csv(out_file, index=False, sep='\t')
    
    return input_dir


def gwas_clump(geno_prefix, p1=0.001, p2=0.05, num_threads=None,
               clump_input_dir=None, result_dir=None, clump_kb=500,
               clump_r2=0.1):
    """GWAS result clumping using PLINK.
    
    Args:
        geno_prefix: PLINK genotype prefix
        p1: Clump p1 threshold
        p2: Clump p2 threshold
        num_threads: Number of threads
    
    Returns:
        Results directory path
    """
    num_threads = _resolve_num_threads(num_threads)
    
    if clump_input_dir is None:
        clump_input_dir = './clump_input'
    if result_dir is None:
        result_dir = (
            os.path.join(os.path.dirname(os.path.abspath(clump_input_dir)), 'clump_result')
            if clump_input_dir != './clump_input'
            else './clump_result'
        )
    # MCP long-job metadata is stored in ``result_dir`` before this worker
    # starts. Removing the whole directory deletes the job file, so preserve
    # MCP-owned files and only clean PLINK outputs produced by clumping.
    os.makedirs(result_dir, exist_ok=True)
    # Remove only stale PLINK *clump* products — NEVER a blanket wipe of
    # result_dir. The old blanket glob deleted sibling files (other traits'
    # ``*.assoc.txt``, ``*.qtl.csv``, ...) whenever result_dir was reused.
    _CLUMP_EXTS = ('.clumped', '.clumped.ranges', '.clumped.best')
    for path in glob.glob(os.path.join(result_dir, '*')):
        name = os.path.basename(path)
        if '.mcp_' in name or name.endswith('.mcp_job.json'):
            continue
        if os.path.isfile(path) and name.endswith(_CLUMP_EXTS):
            os.remove(path)

    clump_files = glob.glob(os.path.join(clump_input_dir, '*'))
    if not clump_files:
        print(f"[warn] gwas_clump: no files in {clump_input_dir}")
        return result_dir

    for fn in clump_files:
        phe_name = os.path.basename(fn).split('.')[0]
        out_name = f"{result_dir}/{phe_name}"

        cmd = f"{PLINK_BIN} --bfile {geno_prefix} --clump {fn} " \
              f"--clump-p1 {p1} --clump-p2 {p2} --clump-kb {clump_kb} --clump-r2 {clump_r2} " \
              f"--out {out_name} --clump-allow-overlap --threads {num_threads}"

        subprocess.run(cmd, shell=True, capture_output=True)

    return result_dir


def qq_plot(pvalues, output_file=None):
    """Generate QQ plot from p-values.
    
    Args:
        pvalues: Array of p-values
        output_file: Output file path
    
    Returns:
        Figure data
    """
    # Remove NaN and > 1
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues <= 1]
    pvalues = pvalues[pvalues > 0]
    
    if len(pvalues) == 0:
        return None
    
    # Observed -log10(p)
    observed = -np.log10(sorted(pvalues))
    
    # Expected -log10(p)
    n = len(pvalues)
    expected = -np.log10(np.arange(1, n+1) / n)
    
    return {'observed': observed, 'expected': expected}


def manhattan_plot(gwas_results, output_file=None, significance=5e-8):
    """Generate Manhattan plot from GWAS results.
    
    Args:
        gwas_results: DataFrame with chr, pos, pvalue columns
        output_file: Output file path
        significance: Significance threshold
    
    Returns:
        Plot data
    """
    # Prepare data
    df = gwas_results.copy()
    
    # Calculate -log10(p)
    if 'pvalue' in df.columns:
        df['neg_log_p'] = -np.log10(df['pvalue'])
    elif 'p_wald' in df.columns:
        df['neg_log_p'] = -np.log10(df['p_wald'])
    elif 'p_score' in df.columns:
        df['neg_log_p'] = -np.log10(df['p_score'])
    
    # Get chromosome and position
    if 'chr' not in df.columns and 'chrom' in df.columns:
        df['chr'] = df['chrom']
    if 'pos' not in df.columns and 'ps' in df.columns:
        df['pos'] = df['ps']
    
    return df


def get_top_snps(gwas_file, n=10, pval_col=None, pvalue_cutoff=None):
    n = int(n)
    if isinstance(gwas_file, pd.DataFrame):
        df = gwas_file.copy()
    else:
        df = pd.read_csv(gwas_file, sep='	')
    
    if pval_col is None:
        for col in ['pvalue', 'P', 'p_wald', 'p_score']:
            if col in df.columns:
                pval_col = col
                break
        if pval_col is None:
            return df.head(n)
    
    if pvalue_cutoff is not None and pval_col in df.columns:
        df = df[df[pval_col] <= pvalue_cutoff]
    
    df = df.sort_values(pval_col)
    return df.head(n)


def get_top_snps_from_df(df, n=10, pval_col='pvalue', pvalue_cutoff=None):
    """Get top SNPs from GWAS DataFrame.
    
    Args:
        df: GWAS results DataFrame
        n: Number of top SNPs
        pval_col: P-value column name
        pvalue_cutoff: Optional p-value cutoff filter
    
    Returns:
        DataFrame with top SNPs
    """
    n = int(n)
    df = df.copy()
    
    # Apply p-value cutoff if specified
    if pvalue_cutoff is not None and pval_col in df.columns:
        df = df[df[pval_col] <= pvalue_cutoff]
    
    df = df.sort_values(pval_col)
    return df.head(n)


def annotate_snps(snp_list, gtf_file=None):
    """Annotate SNP list (placeholder - requires reference genome).
    
    Args:
        snp_list: List of SNP IDs
        gtf_file: GTF annotation file
    
    Returns:
        Annotations dictionary
    """
    # This would require a reference genome annotation
    # Placeholder for now
    return {snp: {'gene': 'unknown', 'annotation': 'unknown'} for snp in snp_list}


def calculate_lambda(pvalues):
    """Calculate genomic inflation factor (lambda).
    
    Args:
        pvalues: Array of p-values
    
    Returns:
        Lambda value
    """
    pvalues = pvalues[~np.isnan(pvalues)]
    pvalues = pvalues[pvalues > 0]
    
    if len(pvalues) == 0:
        return 1.0
    
    # Calculate median chi-square
    chi2_obs = stats.chi2.ppf(1 - pvalues, df=1)
    median_chi2 = np.median(chi2_obs)
    
    # Lambda = median(chi2) / median(chi2 under null)
    lambda_gc = median_chi2 / stats.chi2.ppf(0.5, df=1)
    
    return lambda_gc




def add_rs_id_to_vcf(vcf_file, output_file=None):
    """Add rs IDs to VCF file if not present.
    
    If SNP IDs are missing, generates IDs in format: chr:pos
    
    Args:
        vcf_file: Input VCF file
        output_file: Output VCF file (default: add_rs suffix)
    
    Returns:
        Output file path
    """
    if output_file is None:
        output_file = vcf_file.replace('.vcf', '_rs.vcf')
    
    with open(vcf_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                f_out.write(line)
            else:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    # Check if ID column (col 2) is missing or is '.'
                    if parts[2] == '.' or parts[2] == '':
                        # Generate rs ID as chr:pos
                        rs_id = f"{parts[0]}.s_{parts[1]}"
                        parts[2] = rs_id
                    f_out.write('\t'.join(parts) + '\n')
    
    return output_file


def ensure_rs_id(geno_prefix):
    """Ensure genotype files have SNP rs IDs.
    
    Checks and generates rs IDs for PLINK bim file if missing.
    
    Args:
        geno_prefix: PLINK genotype prefix
    
    Returns:
        True if rs IDs were added, False if already present
    """
    bim_file = geno_prefix + '.bim'
    if not os.path.exists(bim_file):
        return False
    
    # Read bim file
    bim = pd.read_csv(bim_file, sep='\t', header=None,
                     names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
    
    # Check if SNP IDs are missing (typically '.' or numeric only)
    missing_rs = (bim['snp'] == '.') | (bim['snp'].str.match(r'^\d+$', na=False))
    
    if missing_rs.any():
        # Generate rs IDs
        bim.loc[missing_rs, 'snp'] = bim.loc[missing_rs].apply(
            lambda x: f"{x['chr']}.s_{x['pos']}", axis=1
        )
        
        # Save updated bim
        bim.to_csv(bim_file, sep='\t', header=False, index=False)
        return True
    
    return False

__all__ = [
    'read_plink_prefix', 'prepare_gemma_fam',
    'gwas_lm', 'gwas_lmm', 'gwas_plink', 'simple_gwas',
    'generate_clump_input', 'gwas_clump',
    'qq_plot', 'manhattan_plot', 'get_top_snps',
    'annotate_snps', 'calculate_lambda', 'add_rs_id_to_vcf', 'ensure_rs_id'
]
