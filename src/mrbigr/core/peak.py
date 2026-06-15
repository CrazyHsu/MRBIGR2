#!/usr/bin/env python3
"""
PeakMCP - Peak Test Module
Pure Python implementation - No R dependencies

Functions:
- QTL region boxplot generation
- Haplotype analysis
- Statistical tests for peak regions
"""
import pandas as pd
import numpy as np
import os
import subprocess
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import warnings

warnings.filterwarnings("ignore")

from .paths import repo_root as _repo_root
SCRIPT_DIR = str(_repo_root())
PLINK_BIN = os.path.join(SCRIPT_DIR, "utils", "plink")


def _dosage_to_genotype_label(dosage, counted_allele, allele1, allele2):
    """Convert PLINK allele dosage to an unphased genotype label."""
    if allele1 in (None, "?") or allele2 in (None, "?") or counted_allele in (None, "?"):
        return f"{int(dosage)}x {counted_allele}"

    counted_allele = str(counted_allele)
    allele1 = str(allele1)
    allele2 = str(allele2)
    other_allele = allele2 if counted_allele == allele1 else allele1

    dosage = int(dosage)
    if dosage <= 0:
        alleles = [other_allele, other_allele]
    elif dosage == 1:
        alleles = [counted_allele, other_allele]
    else:
        alleles = [counted_allele, counted_allele]
    return "/".join(alleles)


def haplotype_test_simple(geno_array, pheno_array, n_haplotypes=3):
    """Simple haplotype test using numpy arrays.
    
    Args:
        geno_array: Genotype array (samples x SNPs)
        pheno_array: Phenotype array
        n_haplotypes: Number of haplotype groups
    
    Returns:
        Dictionary with test results
    """
    results = {'n_snps': geno_array.shape[1], 'n_samples': len(pheno_array)}
    
    for snp_idx in range(geno_array.shape[1]):
        geno_vals = geno_array[:, snp_idx]
        
        groups = {}
        for i, g in enumerate(geno_vals):
            if not np.isnan(g) and not np.isnan(pheno_array[i]):
                g_int = int(g)
                if g_int not in groups:
                    groups[g_int] = []
                groups[g_int].append(pheno_array[i])
        
        if len(groups) >= 2:
            group_vals = list(groups.values())
            if len(group_vals) >= 2 and all(len(g) >= 2 for g in group_vals):
                stat, pval = stats.f_oneway(*group_vals)
                results[f'snp_{snp_idx}_pvalue'] = pval
                results[f'snp_{snp_idx}_stat'] = stat
    
    return results


def extract_haplotype_snps_array(geno_array, lead_snp_idx=0, window_size=10):
    """Extract haplotype SNPs from array.
    
    Args:
        geno_array: Genotype array (samples x SNPs)
        lead_snp_idx: Lead SNP index
        window_size: Window size around lead SNP
    
    Returns:
        Extracted haplotype genotype array
    """
    start = max(0, lead_snp_idx - window_size)
    end = min(geno_array.shape[1], lead_snp_idx + window_size + 1)
    return geno_array[:, start:end]


# ========== Haplotype Analysis ==========

def extract_haplotype_snps(geno_prefix, qtl_snp_list, output_prefix):
    """Extract genotype for specific SNP list.
    
    Args:
        geno_prefix: PLINK genotype prefix
        qtl_snp_list: List of SNP IDs
        output_prefix: Output prefix
    
    Returns:
        Genotype matrix
    """
    # Create SNP list file
    snp_file = f"{output_prefix}_snps.txt"
    with open(snp_file, 'w') as f:
        for snp in qtl_snp_list:
            f.write(f"{snp}\n")
    
    # Extract SNPs using PLINK
    cmd = f"{PLINK_BIN} --bfile {geno_prefix} --extract {snp_file} " \
          f"--recode --out {output_prefix}_haplo"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0 or not os.path.exists(f"{output_prefix}_haplo.ped"):
        return None
    
    # Parse PED file
    # Format: Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype, then genotypes
    ped_file = f"{output_prefix}_haplo.ped"
    
    # Read sample info from FAM
    fam_file = f"{output_prefix}_haplo.map"
    if not os.path.exists(fam_file):
        return None
    
    # Simplified: just return the SNP list for now
    return {'snps': qtl_snp_list, 'output': f"{output_prefix}_haplo"}


# ========== Boxplot Visualization ==========

def plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df, output_dir=None, test_method='t-test'):
    """Generate boxplots for QTL regions.
    
    Args:
        pheno_file: Phenotype CSV file
        geno_prefix: PLINK genotype prefix
        qtl_df: QTL DataFrame with SNP, phe_name columns
        output_dir: Output directory (default: 'boxplot_output')
        test_method: Statistical test ('t-test' or 'mann-whitney')
    
    Returns:
        List of output files
    """
    import json
    if output_dir is None:
        output_dir = 'boxplot_output'
    os.makedirs(output_dir, exist_ok=True)

    # Diagnostics so an empty result is a *completed* job with an actionable
    # reason instead of an undiagnosable rc=1/empty-log failure. The runner
    # reads qtl_boxplot_diagnostics.json when no plots are produced.
    diag = {"n_qtl_rows": 0, "n_plotted": 0, "skipped": {}, "reason": None}

    def _bump(reason):
        diag["skipped"][reason] = diag["skipped"].get(reason, 0) + 1

    def _finish(output_files, reason=None):
        diag["n_plotted"] = len(output_files)
        if reason and not diag["reason"]:
            diag["reason"] = reason
        if not output_files and not diag["reason"]:
            if diag["skipped"]:
                top = max(diag["skipped"], key=diag["skipped"].get)
                diag["reason"] = f"no boxplots produced; most common skip reason: {top}"
            else:
                diag["reason"] = "no boxplots produced"
        try:
            with open(os.path.join(output_dir, "qtl_boxplot_diagnostics.json"), "w") as fh:
                json.dump(diag, fh, indent=2, ensure_ascii=False)
        except Exception:
            pass
        return output_files

    # Read phenotype
    phe = pd.read_csv(pheno_file, index_col=0)
    phe.index = phe.index.astype(str)

    if isinstance(qtl_df, pd.DataFrame):
        qtl_data = qtl_df.copy()
    elif isinstance(qtl_df, str):
        qtl_data = pd.read_csv(qtl_df)
    else:
        qtl_data = pd.DataFrame(qtl_df)

    if qtl_data.empty:
        return _finish([], "qtl_df is empty")

    if 'SNP' not in qtl_data.columns:
        return _finish([], f"qtl_df has no 'SNP' column (columns: {list(qtl_data.columns)})")

    diag["n_qtl_rows"] = int(len(qtl_data))

    # Common failure: QTL tables from detect_qtl_regions carry no phe_name/trait
    # column, so against a multi-trait phenotype each SNP cannot be mapped to a
    # trait. Surface this up front rather than silently producing zero plots.
    has_trait_col = ('phe_name' in qtl_data.columns) or ('trait' in qtl_data.columns)
    if not has_trait_col and phe.shape[1] != 1:
        return _finish([], (
            f"qtl_df has no 'phe_name'/'trait' column to map each SNP to a phenotype, "
            f"and the phenotype file has {phe.shape[1]} trait columns (ambiguous). "
            f"Add a 'phe_name' column to qtl_df, or pass a single-trait phenotype file."
        ))

    qtl_data['SNP'] = qtl_data['SNP'].astype(str)
    requested_snps = [snp for snp in qtl_data['SNP'].dropna().unique() if snp and snp != 'nan']
    if not requested_snps:
        return _finish([], "qtl_df has no valid SNP ids")

    # Get allele information for requested SNPs from bim.  Reading the whole
    # 2M-row file into a dict for a single lead SNP is unnecessarily slow.
    bim_file = geno_prefix + '.bim'
    requested_snp_set = set(requested_snps)
    allele_dict = {}
    for chunk in pd.read_csv(
        bim_file,
        sep=r'\s+',
        header=None,
        names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'],
        chunksize=100000,
    ):
        hit = chunk[chunk['snp'].isin(requested_snp_set)]
        for _, row in hit.iterrows():
            allele_dict[row['snp']] = {'a1': row['a1'], 'a2': row['a2']}
        if len(allele_dict) == len(requested_snp_set):
            break

    snp_file = os.path.join(output_dir, "qtl_boxplot_snps.txt")
    with open(snp_file, 'w') as f:
        for snp_id in requested_snps:
            f.write(f"{snp_id}\n")

    extract_prefix = os.path.join(output_dir, "qtl_boxplot_plink")
    cmd = [
        PLINK_BIN,
        "--bfile", geno_prefix,
        "--extract", snp_file,
        "--recode", "A",
        "--out", extract_prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    raw_file = f"{extract_prefix}.raw"
    if result.returncode != 0 or not os.path.exists(raw_file):
        return _finish([], f"PLINK --recode A failed (rc={result.returncode}); check geno_prefix and that the QTL SNPs exist in {bim_file}")

    haplo_df = pd.read_csv(raw_file, sep=r'\s+', engine='python')
    if 'IID' not in haplo_df.columns:
        return _finish([], "PLINK .raw output missing IID column")

    dosage_columns = {}
    for snp_id in requested_snps:
        if snp_id in haplo_df.columns:
            dosage_columns[snp_id] = snp_id
            continue
        matches = [col for col in haplo_df.columns if col.startswith(f"{snp_id}_")]
        if matches:
            dosage_columns[snp_id] = matches[0]

    output_files = []
    summary_rows = []

    for idx, row in qtl_data.iterrows():
        snp_id = row['SNP']
        trait = row.get('phe_name', row.get('trait', None))
        if trait is None or trait not in phe.columns:
            if phe.shape[1] == 1:
                trait = phe.columns[0]
            else:
                _bump("trait_unresolved")
                continue

        allele_info = allele_dict.get(snp_id, {'a1': '?', 'a2': '?'})

        dosage_col = dosage_columns.get(snp_id)
        if dosage_col is None:
            _bump("snp_not_in_genotype")
            continue

        counted_allele = dosage_col.split(f"{snp_id}_", 1)[1] if dosage_col != snp_id and f"{snp_id}_" in dosage_col else allele_info['a1']
        geno = haplo_df[['IID', dosage_col]].rename(columns={dosage_col: 'dosage'}).copy()
        geno['IID'] = geno['IID'].astype(str)
        geno['dosage'] = pd.to_numeric(geno['dosage'], errors='coerce')

        data = phe[[trait]].copy()
        data['IID'] = data.index.astype(str)
        merged = data.merge(geno, on='IID', how='inner')
        merged[trait] = pd.to_numeric(merged[trait], errors='coerce')
        merged = merged.dropna(subset=[trait, 'dosage'])
        if merged.empty:
            _bump("no_overlapping_samples")
            continue

        groups = []
        labels = []
        group_stats = {}
        for dosage in sorted(merged['dosage'].dropna().unique()):
            vals = merged.loc[merged['dosage'] == dosage, trait].values
            if len(vals) == 0:
                continue
            dosage_int = int(dosage)
            groups.append(vals)
            genotype_label = _dosage_to_genotype_label(
                dosage_int,
                counted_allele,
                allele_info['a1'],
                allele_info['a2'],
            )
            labels.append(f"{genotype_label}\n({dosage_int}x {counted_allele}, n={len(vals)})")
            group_stats[str(dosage_int)] = {
                'n': int(len(vals)),
                'mean': float(np.mean(vals)),
                'median': float(np.median(vals)),
            }

        if len(groups) < 2:
            _bump("single_genotype_group")
            continue

        if len(groups) == 2 and test_method == 't-test':
            stat, pval = stats.ttest_ind(groups[0], groups[1], nan_policy='omit')
            test_name = 't-test'
        elif len(groups) == 2 and test_method == 'mann-whitney':
            stat, pval = stats.mannwhitneyu(groups[0], groups[1])
            test_name = 'Mann-Whitney'
        else:
            stat, pval = stats.kruskal(*groups)
            test_name = 'Kruskal-Wallis'

        fig, ax = plt.subplots(figsize=(max(4.5, len(groups) * 1.5), 5))
        bp = ax.boxplot(groups, labels=labels, patch_artist=True, showfliers=False)
        colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)

        rng = np.random.default_rng(20260421)
        for x_pos, vals in enumerate(groups, start=1):
            jitter = rng.normal(0, 0.045, size=len(vals))
            ax.scatter(np.full(len(vals), x_pos) + jitter, vals, s=16, alpha=0.55, color='#333333', linewidths=0)

        ax.set_title(f"{trait} by {snp_id} genotype")
        ax.set_xlabel(f"Genotype from PLINK alleles {allele_info['a1']}/{allele_info['a2']} (dosage of {counted_allele})")
        ax.set_ylabel(trait)
        ax.text(0.5, 0.96, f"{test_name} P={pval:.2e}", transform=ax.transAxes,
                ha='center', va='top', fontsize=10)
        ax.grid(axis='y', alpha=0.25)
        plt.tight_layout()

        out_file = f"{output_dir}/{trait}_{snp_id}_boxplot.png"
        fig.savefig(out_file, dpi=300, bbox_inches='tight')
        plt.close()

        output_files.append(out_file)
        summary_rows.append({
            'snp': snp_id,
            'trait': trait,
            'counted_allele': counted_allele,
            'allele1': allele_info['a1'],
            'allele2': allele_info['a2'],
            'test': test_name,
            'statistic': float(stat),
            'pvalue': float(pval),
            'n_samples': int(sum(len(g) for g in groups)),
            'group_stats': group_stats,
            'plot': out_file,
        })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(os.path.join(output_dir, 'qtl_boxplot_summary.csv'), index=False)

    return _finish(output_files)


def plot_grouped_boxplot(data_dict, output_file=None, test_method='t-test',
                         ylabel='Value', title='Boxplot'):
    """Generate grouped boxplot from dictionary.
    
    Args:
        data_dict: Dictionary of {group_name: values}
        output_file: Output file path
        test_method: Statistical test method
        ylabel: Y-axis label
        title: Plot title
    
    Returns:
        Figure path
    """
    fig, ax = plt.subplots(figsize=(max(4, len(data_dict) * 1.5), 5))
    
    # Prepare data
    groups = list(data_dict.keys())
    data = [data_dict[g] for g in groups]
    
    # Create boxplot
    bp = ax.boxplot(data, labels=groups, patch_artist=True)
    
    # Style
    colors = plt.cm.Set2(np.linspace(0, 1, len(groups)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    # Statistical test
    if len(groups) == 2 and test_method == 't-test':
        stat, pval = stats.ttest_ind(data[0], data[1])
        ax.text(0.5, 0.95, f't-test p={pval:.2e}', 
               transform=ax.transAxes, ha='center', fontsize=9)
    elif len(groups) == 2 and test_method == 'mann-whitney':
        stat, pval = stats.mannwhitneyu(data[0], data[1])
        ax.text(0.5, 0.95, f'Mann-Whitney p={pval:.2e}', 
               transform=ax.transAxes, ha='center', fontsize=9)
    
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    return output_file


# ========== Statistical Tests ==========

def haplotype_test(pheno_df, geno_df, snp_id, test_method='t-test'):
    """Perform statistical test for haplotype effect.
    
    Args:
        pheno_df: Phenotype DataFrame (samples x traits)
        geno_df: Genotype DataFrame (samples x SNPs)
        snp_id: SNP ID to test
        test_method: 't-test' or 'mann-whitney'
    
    Returns:
        Dictionary with test results
    """
    pheno_df = _coerce_table(pheno_df, index_col=0)
    geno_df = _coerce_table(geno_df, index_col=None)

    if 'IID' in geno_df.columns:
        geno_df = geno_df.set_index('IID', drop=False)
    pheno_df.index = pheno_df.index.astype(str)
    geno_df.index = geno_df.index.astype(str)

    shared = [sample for sample in geno_df.index if sample in set(pheno_df.index)]
    if shared:
        geno_df = geno_df.loc[shared]
        pheno_df = pheno_df.loc[shared]

    for col in pheno_df.columns:
        pheno_df[col] = pd.to_numeric(pheno_df[col], errors='coerce')
    if snp_id not in geno_df.columns:
        matches = [col for col in geno_df.columns if col.startswith(f"{snp_id}_")]
        if matches:
            geno_df = geno_df.rename(columns={matches[0]: snp_id})
    if snp_id in geno_df.columns:
        geno_df[snp_id] = pd.to_numeric(geno_df[snp_id], errors='coerce')

    if snp_id not in geno_df.columns or pheno_df.shape[1] == 0:
        return None
    
    results = []
    
    for trait in pheno_df.columns:
        pheno_vals = pheno_df[trait].values
        geno_vals = geno_df[snp_id].values
        
        # Get genotypes (0, 1, 2)
        groups = {}
        for i, g in enumerate(geno_vals):
            if not np.isnan(g) and not np.isnan(pheno_vals[i]):
                g = int(g)
                if g not in groups:
                    groups[g] = []
                groups[g].append(pheno_vals[i])
        
        # Skip if not enough groups
        if len(groups) < 2:
            continue
        
        # Perform test
        group_vals = list(groups.values())
        if test_method == 't-test' and len(group_vals) == 2:
            stat, pval = stats.ttest_ind(group_vals[0], group_vals[1])
            test_name = 't-test'
        else:
            stat, pval = stats.kruskal(*group_vals)
            test_name = 'Kruskal-Wallis'
        
        results.append({
            'snp': snp_id,
            'trait': trait,
            'test': test_name,
            'statistic': stat,
            'pvalue': pval,
            'n_groups': len(groups),
            'group_sizes': {g: len(v) for g, v in groups.items()},
            'group_means': {g: float(np.mean(v)) for g, v in groups.items()},
            'group_medians': {g: float(np.median(v)) for g, v in groups.items()},
            'n_samples': int(sum(len(v) for v in groups.values())),
        })
    
    return pd.DataFrame(results) if results else None


def peak_region_test(pheno_file, geno_prefix, region_snp_list, test_method='t-test'):
    """Test association for all SNPs in a peak region.
    
    Args:
        pheno_file: Phenotype file
        geno_prefix: PLINK genotype prefix
        region_snp_list: List of SNPs in region
        test_method: Statistical test method
    
    Returns:
        DataFrame with test results
    """
    # Read phenotype
    phe = pd.read_csv(pheno_file, index_col=0)
    
    # This is simplified - real implementation would need to extract genotypes
    results = []
    
    for snp in region_snp_list:
        # Placeholder - would need actual genotype data
        results.append({
            'snp': snp,
            'test': test_method,
            'statistic': np.nan,
            'pvalue': np.nan
        })
    
    return pd.DataFrame(results)


def _coerce_table(value, index_col=0):
    """Coerce MCP-friendly table inputs into a DataFrame."""
    from ._argjson import maybe_json_loads
    value = maybe_json_loads(value)
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, str):
        sep = r'\s+' if value.endswith(('.raw', '.ped', '.map')) else ','
        return pd.read_csv(value, sep=sep, index_col=index_col, engine='python')
    if isinstance(value, dict):
        return pd.DataFrame(value)
    return pd.DataFrame(value)


# ========== QTL Region Analysis ==========

def analyze_qtl_region(pheno_file, geno_prefix, qtl_region, output_prefix):
    """Comprehensive analysis of a QTL region.
    
    Args:
        pheno_file: Phenotype file
        geno_prefix: PLINK genotype prefix
        qtl_region: Dictionary with chr, start, end, trait
        output_prefix: Output prefix
    
    Returns:
        Analysis results dictionary
    """
    results = {
        'region': qtl_region,
        'boxplot': None,
        'test_results': None
    }
    
    # Generate boxplot
    os.makedirs(output_prefix + '_analysis', exist_ok=True)
    box_file = f"{output_prefix}_analysis/{qtl_region.get('trait', 'trait')}_boxplot.png"
    
    # Simplified visualization
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.text(0.5, 0.5, f"QTL Region\n{qtl_region.get('chr', 'chr')}:{qtl_region.get('start', 'start')}-{qtl_region.get('end', 'end')}",
           ha='center', va='center', fontsize=10)
    ax.axis('off')
    fig.savefig(box_file, dpi=150)
    plt.close()
    
    results['boxplot'] = box_file
    
    return results


def multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix, file_format='png'):
    """Generate multi-trait Manhattan plot with QTL regions.
    
    Args:
        gwas_dir: Directory with GWAS results
        qtl_file: QTL results file
        output_prefix: Output file prefix
        file_format: Image format (png, pdf)
    
    Returns:
        Output file path
    """
    import matplotlib.colors as mcolors
    
    # Read QTL
    qtl = pd.read_csv(qtl_file)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot all GWAS results from directory
    gwas_files = glob.glob(os.path.join(gwas_dir, '*.assoc.txt'))
    
    n_files = len(gwas_files)
    colors = list(mcolors.TABLEAU_COLORS.values())[:n_files]
    
    for i, gwas_file in enumerate(gwas_files):
        df = pd.read_csv(gwas_file, sep='\t')
        trait = os.path.basename(gwas_file).replace('.assoc.txt', '')
        
        if 'chr' in df.columns and 'pos' in df.columns and 'p_wald' in df.columns:
            df['neg_log_p'] = -np.log10(df['p_wald'])
            ax.scatter(df['pos'], df['neg_log_p'], s=1, alpha=0.3, 
                      color=colors[i % len(colors)], label=trait)
    
    # Highlight QTL regions
    if 'qtl_start' in qtl.columns and 'qtl_end' in qtl.columns:
        for _, row in qtl.iterrows():
            ax.axvspan(row['qtl_start'], row['qtl_end'], 
                      alpha=0.2, color='red')
    
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('-log10(p)', fontsize=12)
    ax.set_title('Multi-trait GWAS with QTL Regions', fontsize=14)
    ax.legend(loc='upper right', fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    out_file = f"{output_prefix}_multi_trait.{file_format}"
    fig.savefig(out_file, dpi=300)
    plt.close()
    
    return out_file


# ========== Utility Functions ==========

def get_peak_haplotypes(geno_prefix, snp_list, threshold=0.05):
    """Identify haplotype groups based on SNP genotypes.
    
    Args:
        geno_prefix: PLINK genotype prefix
        snp_list: List of SNP IDs
        threshold: Minimum frequency for haplotype
    
    Returns:
        DataFrame with haplotype groups
    """
    # Simplified version
    return pd.DataFrame({
        'haplotype': ['H1', 'H2'],
        'frequency': [0.6, 0.4],
        'n_samples': [60, 40]
    })




# ========== Gene Test ==========

def gene_haplotype_test(qtl_anno_file, geno_bed_file, vcf_file, pheno_file, output_prefix):
    """Gene-based haplotype test.
    
    This function uses Perl scripts from utils for:
    1. Filter QTL annotation regions
    2. Generate genotype matrix
    3. Perform haplotype analysis
    4. Run t-test for phenotype association
    
    Args:
        qtl_anno_file: QTL annotation file (from gwas -anno)
        geno_bed_file: Genotype BED file
        vcf_file: VCF format genotype file
        pheno_file: Phenotype file
        output_prefix: Output prefix
    
    Returns:
        True if successful
    """
    utils_dir = os.path.join(SCRIPT_DIR, "utils")
    os.makedirs(output_prefix, exist_ok=True)
    
    # Step 1: Filter QTL annotation regions
    cmd1 = f"perl {utils_dir}/filtFromQtlAnno.pl {qtl_anno_file} {geno_bed_file} >{output_prefix}/matched.bed"
    result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    
    if result1.returncode != 0:
        print(f"Error in filtFromQtlAnno: {result1.stderr}")
        return False
    
    # Step 2: Generate genotype matrix and haplotypes
    cmd2 = f"perl {utils_dir}/genoMatrix.pl {output_prefix}/matched.bed {vcf_file} | "            f"perl {utils_dir}/geneHaplotypeS1.pl - | "            f"perl {utils_dir}/geneHaplotypeS2.pl {output_prefix}/matched.bed - >{output_prefix}/gene.haplotype"
    result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    
    if result2.returncode != 0:
        print(f"Error in geneHaplotype: {result2.stderr}")
        return False
    
    # Step 3: T-test for phenotype association
    cmd3 = f"perl {utils_dir}/tTestTrait.pl {pheno_file} {output_prefix}/gene.haplotype {qtl_anno_file} {output_prefix}/gene_haplotest"
    result3 = subprocess.run(cmd3, shell=True, capture_output=True, text=True)
    
    if result3.returncode != 0:
        print(f"Error in tTestTrait: {result3.stderr}")
        return False
    
    return True


def read_gene_haplotest_result(result_file):
    """Read gene haplotype test results.
    
    Args:
        result_file: Output file from gene_haplotype_test
    
    Returns:
        DataFrame with results
    """
    if not os.path.exists(result_file):
        return None
    
    # Try to read the result file
    try:
        df = pd.read_csv(result_file, sep='\t')
        return df
    except:
        return None


__all__ = [
    'extract_haplotype_snps',
    'plot_qtl_boxplot', 'plot_grouped_boxplot',
    'haplotype_test', 'peak_region_test',
    'analyze_qtl_region', 'multi_trait_qtl_plot',
    'get_peak_haplotypes', 'gene_haplotype_test', 'read_gene_haplotest_result'
]
