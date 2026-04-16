#!/usr/bin/env python3
"""
QTLMCP - QTL Analysis Module
Pure Python implementation - No R dependencies

Functions:
- QTL interval detection from GWAS results
- Peak SNP identification
- QTL region visualization
- Basic haplotype analysis
"""
import pandas as pd
import numpy as np
import os
import subprocess
import glob
import shutil
import warnings

warnings.filterwarnings("ignore")

# External tools
from .paths import repo_root as _repo_root
SCRIPT_DIR = str(_repo_root())
PLINK_BIN = os.path.join(SCRIPT_DIR, "utils", "plink")


def _coerce_qtl_df(qtl_df):
    """Coerce structured inputs into a QTL DataFrame."""
    if isinstance(qtl_df, pd.DataFrame):
        return qtl_df.copy()
    if isinstance(qtl_df, dict):
        return pd.DataFrame(qtl_df)
    return pd.DataFrame(qtl_df)


def _ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


# ========== QTL Detection ==========

def detect_qtl(gwas_file, p1=1e-7, p2=1e-5, p2n=5, window=500000):
    """Detect QTL regions from GWAS results.
    
    Args:
        gwas_file: GWAS results file path or DataFrame
        p1: Strong significance threshold (index SNP threshold)
        p2: Suggestive threshold (secondary SNP threshold)
        p2n: Minimum number of secondary significant SNPs in a QTL
        window: Window for merging nearby SNPs
    
    Returns:
        DataFrame with QTL intervals
    """
    if isinstance(gwas_file, pd.DataFrame):
        df = gwas_file.copy()
    else:
        df = pd.read_csv(gwas_file, sep='\t')
    
    # Determine column names
    if 'p_wald' in df.columns:
        p_col = 'p_wald'
    elif 'p_score' in df.columns:
        p_col = 'p_score'
    else:
        return None
    
    if 'chr' in df.columns:
        chr_col = 'chr'
    elif 'chrom' in df.columns:
        chr_col = 'chrom'
    else:
        return None
    
    if 'pos' in df.columns:
        pos_col = 'pos'
    elif 'ps' in df.columns:
        pos_col = 'ps'
    else:
        return None
    
    # Candidate SNPs include both strong and suggestive hits.
    candidate_snps = df[df[p_col] < p2].copy()
    if len(candidate_snps) == 0:
        return pd.DataFrame()
    
    # Group by chromosome and merge nearby regions
    qtls = []
    for chr_val in candidate_snps[chr_col].unique():
        chr_snps = candidate_snps[candidate_snps[chr_col] == chr_val].sort_values(pos_col)
        
        # Build candidate clusters from suggestive SNPs.
        current_cluster = None
        for _, snp in chr_snps.iterrows():
            if current_cluster is None:
                current_cluster = [snp]
            elif snp[pos_col] - current_cluster[-1][pos_col] < window:
                current_cluster.append(snp)
            else:
                cluster_df = pd.DataFrame(current_cluster)
                n_secondary = int((cluster_df[p_col] < p2).sum())
                has_strong = bool((cluster_df[p_col] < p1).any())
                if has_strong or n_secondary >= p2n:
                    lead = cluster_df.loc[cluster_df[p_col].idxmin()]
                    qtls.append({
                        'CHR': chr_val,
                        'qtl_start': int(cluster_df[pos_col].min()),
                        'qtl_end': int(cluster_df[pos_col].max()),
                        'SNP': lead.get('rs', lead.get('snp', f"{chr_val}:{lead[pos_col]}")),
                        'P': lead[p_col],
                        'qtl_length': int(cluster_df[pos_col].max() - cluster_df[pos_col].min())
                    })
                current_cluster = [snp]
        
        if current_cluster is not None:
            cluster_df = pd.DataFrame(current_cluster)
            n_secondary = int((cluster_df[p_col] < p2).sum())
            has_strong = bool((cluster_df[p_col] < p1).any())
            if has_strong or n_secondary >= p2n:
                lead = cluster_df.loc[cluster_df[p_col].idxmin()]
                qtls.append({
                    'CHR': chr_val,
                    'qtl_start': int(cluster_df[pos_col].min()),
                    'qtl_end': int(cluster_df[pos_col].max()),
                    'SNP': lead.get('rs', lead.get('snp', f"{chr_val}:{lead[pos_col]}")),
                    'P': lead[p_col],
                    'qtl_length': int(cluster_df[pos_col].max() - cluster_df[pos_col].min())
                })
    
    result = pd.DataFrame(qtls)
    if len(result) > 0:
        result['qtl_length'] = result['qtl_end'] - result['qtl_start']
    
    return result


def get_lead_snp(gwas_file, region_chr, region_start, region_end):
    """Get lead SNP (most significant) in a region.
    
    Args:
        gwas_file: GWAS results
        region_chr: Chromosome
        region_start: Region start
        region_end: Region end
    
    Returns:
        DataFrame with lead SNP info
    """
    df = pd.read_csv(gwas_file, sep='\t')
    
    # Column mapping
    chr_col = 'chr' if 'chr' in df.columns else 'chrom'
    pos_col = 'pos' if 'pos' in df.columns else 'ps'
    p_col = 'p_wald' if 'p_wald' in df.columns else 'p_score'
    
    # Filter to region
    region_snps = df[
        (df[chr_col].astype(str) == str(region_chr)) &
        (df[pos_col] >= region_start) &
        (df[pos_col] <= region_end)
    ]
    
    region_snps = region_snps.dropna(subset=[p_col])
    if len(region_snps) == 0:
        return None

    lead = region_snps.loc[region_snps[p_col].idxmin()]
    
    return pd.DataFrame([{
        'chr': lead[chr_col],
        'pos': lead[pos_col],
        'snp': lead.get('rs', lead.get('snp', 'unknown')),
        'ref': lead.get('a1', 'N/A'),
        'alt': lead.get('a0', 'N/A'),
        'beta': lead.get('beta', lead.get('eff', 0)),
        'pvalue': lead[p_col],
        'af': lead.get('af', 0)
    }])


def identify_peak_snps(gwas_dir, p_threshold=1e-5, max_peaks=50):
    """Identify peak SNPs from multiple GWAS results.
    
    Args:
        gwas_dir: Directory with GWAS results
        p_threshold: P-value threshold
        max_peaks: Maximum number of peaks to return
    
    Returns:
        Combined DataFrame of peak SNPs
    """
    peak_list = []
    
    for gwas_file in glob.glob(os.path.join(gwas_dir, '*.assoc.txt')):
        df = pd.read_csv(gwas_file, sep='\t')
        
        # Determine columns
        p_col = 'p_wald' if 'p_wald' in df.columns else 'p_score'
        chr_col = 'chr' if 'chr' in df.columns else 'chrom'
        pos_col = 'pos' if 'pos' in df.columns else 'ps'
        
        # Get significant SNPs
        sig = df[df[p_col] < p_threshold]
        
        # Add trait name from filename
        trait = os.path.basename(gwas_file).replace('.assoc.txt', '')
        sig['trait'] = trait
        
        peaks = sig[[chr_col, pos_col, 'rs', p_col, 'beta']].copy()
        peaks.columns = ['chr', 'pos', 'snp', 'pvalue', 'beta']
        peak_list.append(peaks)
    
    if len(peak_list) == 0:
        return pd.DataFrame()
    
    combined = pd.concat(peak_list, ignore_index=True)
    combined = combined.sort_values('pvalue').head(max_peaks)
    
    return combined


# ========== QTL Region Analysis ==========

def extract_qtl_genotypes(geno_prefix, qtl_regions, output_prefix):
    """Extract genotypes for QTL regions.
    
    Args:
        geno_prefix: PLINK genotype prefix
        qtl_regions: DataFrame with qtl_start, qtl_end, CHR columns
        output_prefix: Output file prefix
    
    Returns:
        List of output file paths
    """
    qtl_regions = _coerce_qtl_df(qtl_regions)
    _ensure_parent_dir(output_prefix)
    output_files = []
    
    for idx, qtl in qtl_regions.iterrows():
        chr_val = qtl['CHR']
        start = int(qtl['qtl_start'])
        end = int(qtl['qtl_end'])
        
        # Use PLINK to extract region
        output_file = f"{output_prefix}_qtl_{idx}"
        cmd = f"{PLINK_BIN} --bfile {geno_prefix} --chr {chr_val} " \
              f"--from-bp {start} --to-bp {end} " \
              f"--make-bed --out {output_file}"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0 and os.path.exists(f"{output_file}.bed"):
            output_files.append(output_file)
    
    return output_files


def calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix):
    """Calculate haplotypes for given SNP list.
    
    Args:
        geno_prefix: PLINK genotype prefix
        qtl_snp_list: List of SNP IDs
        output_prefix: Output prefix
    
    Returns:
        Haplotype DataFrame
    """
    _ensure_parent_dir(output_prefix)
    # Create SNP list file
    snp_file = f"{output_prefix}_snps.txt"
    with open(snp_file, 'w') as f:
        for snp in qtl_snp_list:
            f.write(f"{snp}\n")
    
    # Extract SNPs
    cmd = f"{PLINK_BIN} --bfile {geno_prefix} --extract {snp_file} " \
          f"--recode --out {output_prefix}_haplo"
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    # Read PED file and calculate haplotypes
    # This is a simplified version
    ped_file = f"{output_prefix}_haplo.ped"
    if not os.path.exists(ped_file):
        return None
    
    # Parse PED file (simplified)
    haplo_df = pd.read_csv(ped_file, sep=r'\s+', header=None, engine='python')
    
    # Calculate haplotype (each row has 2 alleles per SNP)
    # This is a basic implementation
    return haplo_df


# ========== QTL Visualization ==========

def plot_qtl_region(gwas_file, chr_val, start, end, output_file=None, 
                   highlight_snps=None):
    """Plot GWAS results for a specific QTL region.
    
    Args:
        gwas_file: GWAS results file
        chr_val: Chromosome
        start: Region start
        end: Region end
        output_file: Output file path
        highlight_snps: List of SNP IDs to highlight
    
    Returns:
        Figure path
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    df = pd.read_csv(gwas_file, sep='\t')
    
    # Filter to region
    chr_col = 'chr' if 'chr' in df.columns else 'chrom'
    pos_col = 'pos' if 'pos' in df.columns else 'ps'
    p_col = 'p_wald' if 'p_wald' in df.columns else 'p_score'
    
    region = df[
        (df[chr_col].astype(str) == str(chr_val)) &
        (df[pos_col] >= start) &
        (df[pos_col] <= end)
    ]
    
    if len(region) == 0:
        return None
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    region = region.sort_values(pos_col)
    ax.scatter(region[pos_col], -np.log10(region[p_col]), s=10, c='steelblue')
    
    # Highlight specific SNPs
    if highlight_snps:
        highlight = region[region['rs'].isin(highlight_snps)]
        ax.scatter(highlight[pos_col], -np.log10(highlight[p_col]), 
                  s=50, c='red', zorder=5)
    
    # Significance line
    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.5)
    ax.axhline(y=-np.log10(1e-5), color='blue', linestyle='--', alpha=0.5)
    
    ax.set_xlabel(f"Position (bp)", fontsize=12)
    ax.set_ylabel("-log10(p)", fontsize=12)
    ax.set_title(f"QTL Region: Chr{chr_val}:{start}-{end}", fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_file is None:
        output_file = f"qtl_chr{chr_val}_{start}_{end}.png"
    
    _ensure_parent_dir(output_file)
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_file


def qtl_summary(qtl_dir, output_prefix='qtl_summary'):
    """Generate summary of all QTL regions.
    
    Args:
        qtl_dir: Directory with GWAS results
        output_prefix: Output file prefix
    
    Returns:
        DataFrame with all QTL
    """
    all_qtls = []
    
    for gwas_file in glob.glob(os.path.join(qtl_dir, '*.assoc.txt')):
        trait = os.path.basename(gwas_file).replace('.assoc.txt', '')
        
        qtls = detect_qtl(gwas_file)
        if qtls is not None and len(qtls) > 0:
            qtls['trait'] = trait
            all_qtls.append(qtls)
    
    if len(all_qtls) == 0:
        return pd.DataFrame()
    
    result = pd.concat(all_qtls, ignore_index=True)
    
    _ensure_parent_dir(output_prefix)
    # Save to file
    result.to_csv(f"{output_prefix}.csv", index=False)
    
    return result


# ========== QTL to Gene Mapping ==========

def map_qtl_to_genes(qtl_df, annotation_file):
    """Map QTL regions to genes using annotation.
    
    Args:
        qtl_df: QTL DataFrame with CHR, qtl_start, qtl_end columns
        annotation_file: GTF/GFF or BED annotation file
    
    Returns:
        QTL with gene annotations
    """
    from anno import get_genes_in_region
    
    qtl_df = _coerce_qtl_df(qtl_df)
    results = []
    for _, qtl in qtl_df.iterrows():
        genes = get_genes_in_region(
            qtl['CHR'], 
            qtl['qtl_start'], 
            qtl['qtl_end'], 
            annotation_file
        )
        
        result = qtl.to_dict()
        if genes is not None and not genes.empty:
            result['n_genes'] = len(genes)
            gene_ids = genes['gene_id'].dropna().tolist() if 'gene_id' in genes.columns else []
            result['genes'] = ';'.join(gene_ids) if gene_ids else 'none'
        else:
            result['n_genes'] = 0
            result['genes'] = 'none'
        results.append(result)
    
    return pd.DataFrame(results)


# ========== QTL Enrichment (Basic) ==========

def qtl_enrichment_test(qtl_genes, background_genes, gene_sets):
    """Test QTL gene enrichment (simplified).
    
    Args:
        qtl_genes: List of genes in QTL
        background_genes: Background gene list
        gene_sets: Dictionary of gene set names to gene lists
    
    Returns:
        DataFrame with enrichment results
    """
    from scipy.stats import fisher_exact
    
    results = []
    
    qtl_set = set(qtl_genes)
    bg_set = set(background_genes)
    
    for gs_name, gs_genes in gene_sets.items():
        gs_set = set(gs_genes)
        
        # Contingency table
        #            In QTL   Not in QTL
        # In GS        a         b
        # Not in GS    c         d
        
        a = len(qtl_set & gs_set)
        b = len(gs_set - qtl_set)
        c = len(qtl_set - gs_set)
        d = len(bg_set - qtl_set - gs_set)
        
        if a > 0 and d > 0:
            odds, pval = fisher_exact([[a, b], [c, d]])
            
            results.append({
                'gene_set': gs_name,
                'overlap': a,
                'total_in_gs': len(gs_set),
                'odds_ratio': odds,
                'pvalue': pval
            })
    
    if len(results) == 0:
        return pd.DataFrame(columns=['gene_set', 'overlap', 'total_in_gs', 'odds_ratio', 'pvalue'])

    return pd.DataFrame(results).sort_values('pvalue')


# ========== Utility Functions ==========

def get_qtl_overlap(qtl1, qtl2):
    """Find overlapping QTL between two QTL sets.
    
    Args:
        qtl1: First QTL DataFrame
        qtl2: Second QTL DataFrame
    
    Returns:
        DataFrame with overlapping QTL
    """
    overlaps = []
    
    for _, q1 in qtl1.iterrows():
        for _, q2 in qtl2.iterrows():
            if q1['CHR'] == q2['CHR']:
                # Check overlap
                start = max(q1['qtl_start'], q2['qtl_start'])
                end = min(q1['qtl_end'], q2['qtl_end'])
                
                if end > start:
                    overlaps.append({
                        'CHR': q1['CHR'],
                        'start': start,
                        'end': end,
                        'qtl1_trait': q1.get('phe_name', q1.get('trait', 'unknown')),
                        'qtl2_trait': q2.get('phe_name', q2.get('trait', 'unknown')),
                        'length': end - start
                    })
    
    if len(overlaps) == 0:
        return pd.DataFrame(columns=['CHR', 'start', 'end', 'qtl1_trait', 'qtl2_trait', 'length'])

    return pd.DataFrame(overlaps)


def export_qtl_bed(qtl_df, output_file):
    """Export QTL as BED file for visualization.
    
    Args:
        qtl_df: QTL DataFrame
        output_file: Output BED file
    
    Returns:
        File path
    """
    qtl_df = _coerce_qtl_df(qtl_df)
    _ensure_parent_dir(output_file)
    bed = qtl_df[['CHR', 'qtl_start', 'qtl_end']].copy()
    if 'phe_name' in qtl_df.columns:
        bed['name'] = qtl_df['phe_name']
    elif 'trait' in qtl_df.columns:
        bed['name'] = qtl_df['trait']
    else:
        bed['name'] = 'QTL'
    
    bed['score'] = -np.log10(qtl_df['P']) if 'P' in qtl_df.columns else 0
    bed['strand'] = '.'
    
    bed.to_csv(output_file, sep='\t', header=False, index=False)
    
    return output_file


__all__ = [
    'detect_qtl', 'get_lead_snp', 'identify_peak_snps',
    'extract_qtl_genotypes', 'calculate_qtl_haplo',
    'plot_qtl_region', 'qtl_summary',
    'map_qtl_to_genes', 'qtl_enrichment_test',
    'get_qtl_overlap', 'export_qtl_bed'
]
