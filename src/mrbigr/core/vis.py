#!/usr/bin/env python3
"""
VisMCP - Visualization Module
Pure Python implementation - No R dependencies

Functions:
- Manhattan plot (replacing CMplot.R)
- QQ plot
- PCA/t-SNE scatter plots
- LD heatmap
- Phenotype distribution plots
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats
from scipy.stats import norm
import os
import warnings

warnings.filterwarnings("ignore")


# ========== GWAS Visualization ==========

def manhattan_plot(gwas_file, output_file=None, significance=5e-8, suggest=1e-5, 
                   chrom_col='chr', pos_col='pos', pval_col='p_wald',
                   title='Manhattan Plot', dpi=300):
    """Generate Manhattan plot from GWAS results.
    
    Args:
        gwas_file: GWAS results file (TSV)
        output_file: Output file path
        significance: Genome-wide significance threshold (default 5e-8)
        suggest: Suggestive threshold (default 1e-5)
        chrom_col: Chromosome column name
        pos_col: Position column name
        pval_col: P-value column name
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    # Read GWAS results
    df = pd.read_csv(gwas_file, sep='\t')
    
    # Check columns
    if chrom_col not in df.columns:
        # Try alternative names
        if 'chrom' in df.columns:
            chrom_col = 'chrom'
        elif '#chr' in df.columns:
            chrom_col = '#chr'
    
    if pos_col not in df.columns:
        if 'ps' in df.columns:
            pos_col = 'ps'
        elif 'bp' in df.columns:
            pos_col = 'bp'
    
    if pval_col not in df.columns:
        if 'p_score' in df.columns:
            pval_col = 'p_score'
        elif 'pvalue' in df.columns:
            pval_col = 'pvalue'
    
    # Calculate -log10(p)
    df['neg_log_p'] = -np.log10(df[pval_col])
    
    # Handle chromosome sorting
    chr_order = []
    for c in df[chrom_col].unique():
        try:
            chr_order.append((int(c), c))
        except:
            chr_order.append((999, c))  # Put non-numeric chromosomes at end
    
    chr_order.sort(key=lambda x: x[0])
    chr_map = {c: i for i, (_, c) in enumerate(chr_order)}
    
    df['chr_num'] = df[chrom_col].map(chr_map)
    df = df.sort_values(['chr_num', pos_col])
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Color scheme
    colors = ['#1f77b4', '#aec7e8']  # Blue, light blue
    
    # Plot by chromosome with cumulative x offset
    offset = 0
    chr_centers = []
    chr_labels = []

    for i, (chr_num, chr_data) in enumerate(df.groupby('chr_num')):
        chr_name = chr_data[chrom_col].iloc[0]
        color = colors[i % len(colors)]
        n = len(chr_data)

        x_pos = np.arange(n) + offset
        y_vals = chr_data['neg_log_p'].values

        ax.scatter(x_pos, y_vals, c=color, s=3, alpha=0.7)

        chr_labels.append(str(chr_name))
        chr_centers.append(offset + n // 2)
        offset += n

    # Add significance lines
    ax.axhline(y=-np.log10(significance), color='red', linestyle='--',
               linewidth=1, label=f'GWAS ({significance})')
    ax.axhline(y=-np.log10(suggest), color='blue', linestyle='--',
               linewidth=1, label=f'Suggestive ({suggest})')

    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('-log10(p)', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(chr_centers)
    ax.set_xticklabels(chr_labels, rotation=45, fontsize=8)
    ax.set_xlim(0, offset)
    
    plt.tight_layout()
    
    # Save
    if output_file is None:
        output_file = 'Manhattan.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


def qq_plot(gwas_file, output_file=None, pval_col='p_wald', title='Q-Q Plot', dpi=300):
    """Generate Q-Q plot from GWAS results.
    
    Args:
        gwas_file: GWAS results file (TSV)
        output_file: Output file path
        pval_col: P-value column name
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    # Read GWAS results
    df = pd.read_csv(gwas_file, sep='\t')
    
    if pval_col not in df.columns:
        if 'p_score' in df.columns:
            pval_col = 'p_score'
        elif 'pvalue' in df.columns:
            pval_col = 'pvalue'
    
    # Get p-values
    pvals = df[pval_col].dropna()
    pvals = pvals[pvals > 0]
    pvals = pvals[pvals <= 1]
    
    if len(pvals) == 0:
        return None
    
    # Observed -log10(p)
    observed = -np.log10(sorted(pvals))
    
    # Expected -log10(p)
    n = len(pvals)
    expected = -np.log10(np.arange(1, n+1) / n)
    
    # Calculate lambda
    chi2_obs = stats.chi2.ppf(1 - pvals, df=1)
    median_chi2 = np.median(chi2_obs)
    lambda_gc = median_chi2 / stats.chi2.ppf(0.5, df=1)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot
    ax.scatter(expected, observed, c='blue', s=10, alpha=0.5)
    
    # Add diagonal line
    max_val = max(max(expected), max(observed))
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1, label='Expected')
    
    # Labels
    ax.set_xlabel('Expected -log10(p)', fontsize=12)
    ax.set_ylabel('Observed -log10(p)', fontsize=12)
    ax.set_title(f'{title}\nλ = {lambda_gc:.3f}', fontsize=14)
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save
    if output_file is None:
        output_file = 'QQplot.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


# ========== Population Structure Visualization ==========

def pca_plot(pca_file, output_file=None, pc_x='PC1', pc_y='PC2', 
             color_by=None, title='PCA Plot', dpi=300):
    """Generate PCA scatter plot.
    
    Args:
        pca_file: PCA results CSV file
        output_file: Output file path
        pc_x: X-axis PC (default PC1)
        pc_y: Y-axis PC (default PC2)
        color_by: Column to color by
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    df = pd.read_csv(pca_file, index_col=0)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if color_by and color_by in df.columns:
        # Color by category
        categories = df[color_by].unique()
        cmap = plt.cm.get_cmap('tab10', len(categories))
        for i, cat in enumerate(categories):
            mask = df[color_by] == cat
            ax.scatter(df.loc[mask, pc_x], df.loc[mask, pc_y], 
                      c=[cmap(i)], s=30, alpha=0.7, label=str(cat))
        ax.legend()
    else:
        ax.scatter(df[pc_x], df[pc_y], c='steelblue', s=30, alpha=0.7)
    
    ax.set_xlabel(pc_x, fontsize=12)
    ax.set_ylabel(pc_y, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'PCA_plot.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


def tsne_plot(tsne_file, output_file=None, dim_x='dim1', dim_y='dim2',
              color_by=None, title='t-SNE Plot', dpi=300):
    """Generate t-SNE scatter plot.
    
    Args:
        tsne_file: t-SNE results CSV file
        output_file: Output file path
        dim_x: X-axis dimension
        dim_y: Y-axis dimension
        color_by: Column to color by
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    df = pd.read_csv(tsne_file, index_col=0)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if color_by and color_by in df.columns:
        categories = df[color_by].unique()
        cmap = plt.cm.get_cmap('tab10', len(categories))
        for i, cat in enumerate(categories):
            mask = df[color_by] == cat
            ax.scatter(df.loc[mask, dim_x], df.loc[mask, dim_y], 
                      c=[cmap(i)], s=30, alpha=0.7, label=str(cat))
        ax.legend()
    else:
        ax.scatter(df[dim_x], df[dim_y], c='steelblue', s=30, alpha=0.7)
    
    ax.set_xlabel(dim_x, fontsize=12)
    ax.set_ylabel(dim_y, fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'tSNE_plot.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


# ========== Genotype Visualization ==========

def ld_heatmap(geno_prefix, output_file=None, max_snps=500, dpi=300):
    """Generate LD heatmap from genotype data.
    
    Args:
        geno_prefix: PLINK genotype prefix
        output_file: Output file path
        max_snps: Maximum SNPs to include
        dpi: Resolution
    
    Returns:
        Figure path
    """
    from geno import read_plink_bed
    
    # Read genotype data
    G, snps, fam = read_plink_bed(geno_prefix, max_snps=max_snps)
    
    G_clean = G.astype(float)
    G_clean[G_clean == -1] = np.nan
    
    # Calculate LD matrix (correlation)
    # Handle missing values
    n_snps = G_clean.shape[0]
    ld_matrix = np.zeros((n_snps, n_snps))
    
    for i in range(n_snps):
        for j in range(i, n_snps):
            valid = ~np.isnan(G_clean[i]) & ~np.isnan(G_clean[j])
            if valid.sum() > 10:
                corr = np.corrcoef(G_clean[i, valid], G_clean[j, valid])[0, 1]
                ld_matrix[i, j] = corr
                ld_matrix[j, i] = corr
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot heatmap
    im = ax.imshow(ld_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('r²', fontsize=12)
    
    # Labels
    ax.set_xlabel('SNP', fontsize=12)
    ax.set_ylabel('SNP', fontsize=12)
    ax.set_title('Linkage Disequilibrium Heatmap', fontsize=14)
    
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'LD_heatmap.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


# ========== Phenotype Visualization ==========

def phenotype_hist(pheno_file, output_file=None, title='Phenotype Distribution', dpi=300):
    """Generate histogram for each phenotype.
    
    Args:
        pheno_file: Phenotype CSV file
        output_file: Output file path
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    df = pd.read_csv(pheno_file, index_col=0)
    
    n_traits = df.shape[1]
    n_cols = min(4, n_traits)
    n_rows = (n_traits + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
    if n_traits == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for i, col in enumerate(df.columns):
        ax = axes[i]
        data = df[col].dropna()
        ax.hist(data, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_xlabel(col, fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'{col}: μ={data.mean():.2f}, σ={data.std():.2f}', fontsize=10)
    
    # Hide empty subplots
    for i in range(n_traits, len(axes)):
        axes[i].set_visible(False)
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'phenotype_hist.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


def phenotype_boxplot(pheno_file, output_file=None, group_col=None, title='Phenotype Boxplot', dpi=300):
    """Generate boxplot for phenotypes.
    
    Args:
        pheno_file: Phenotype CSV file
        output_file: Output file path
        group_col: Column to group by
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    df = pd.read_csv(pheno_file, index_col=0)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create boxplot data
    data_to_plot = []
    labels = []
    
    for col in df.columns:
        data_to_plot.append(df[col].dropna().values)
        labels.append(col)
    
    bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
    
    # Style
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    
    ax.set_ylabel('Value', fontsize=12)
    ax.set_title(title, fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'phenotype_boxplot.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


def phenotype_correlation(pheno_file, output_file=None, method='pearson', title='Phenotype Correlation', dpi=300):
    """Generate phenotype correlation heatmap.
    
    Args:
        pheno_file: Phenotype CSV file
        output_file: Output file path
        method: Correlation method ('pearson' or 'spearman')
        title: Plot title
        dpi: Resolution
    
    Returns:
        Figure path
    """
    df = pd.read_csv(pheno_file, index_col=0)
    
    # Calculate correlation matrix
    if method == 'spearman':
        corr = df.corr(method='spearman')
    else:
        corr = df.corr(method='pearson')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot heatmap
    im = ax.imshow(corr.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Correlation', fontsize=12)
    
    # Labels
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(corr.columns, fontsize=8)
    
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    
    if output_file is None:
        output_file = 'phenotype_correlation.png'
    fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return output_file


# ========== Multi-function Plotting ==========

def gwas_summary_plot(gwas_dir, output_prefix='gwas_summary', dpi=300):
    """Generate comprehensive GWAS summary plots (Manhattan + QQ).
    
    Args:
        gwas_dir: Directory containing GWAS results
        output_prefix: Output file prefix
        dpi: Resolution
    
    Returns:
        Dictionary with output paths
    """
    results = {}
    
    # Find GWAS files
    gwas_files = [f for f in os.listdir(gwas_dir) if f.endswith('.assoc.txt')]
    
    for gwas_file in gwas_files:
        gwas_path = os.path.join(gwas_dir, gwas_file)
        trait_name = gwas_file.replace('.assoc.txt', '')
        
        # Manhattan plot
        manhattan_out = f"{output_prefix}_{trait_name}_manhattan.png"
        manhattan_plot(gwas_path, manhattan_out, title=f'Manhattan Plot - {trait_name}', dpi=dpi)
        results[f'{trait_name}_manhattan'] = manhattan_out
        
        # QQ plot
        qq_out = f"{output_prefix}_{trait_name}_qq.png"
        qq_plot(gwas_path, qq_out, title=f'Q-Q Plot - {trait_name}', dpi=dpi)
        results[f'{trait_name}_qq'] = qq_out
    
    return results


__all__ = [
    'manhattan_plot', 'qq_plot',
    'pca_plot', 'tsne_plot',
    'ld_heatmap',
    'phenotype_hist', 'phenotype_boxplot', 'phenotype_correlation',
    'gwas_summary_plot'
]