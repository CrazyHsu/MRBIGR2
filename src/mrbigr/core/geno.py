#!/usr/bin/env python3
"""
GenoMCP - Genotype Processing Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/geno.py
"""
import pandas as pd
import numpy as np
import os
import struct
import subprocess
import shutil
from pathlib import Path
from collections import Counter
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import warnings

warnings.filterwarnings("ignore")


# Directories - use relative paths from script location
from .paths import repo_root as _repo_root
SCRIPT_DIR = _repo_root()
OUTPUT_DIR = SCRIPT_DIR / "output"
DATA_DIR = SCRIPT_DIR / "data"
PLINK_BIN = str(SCRIPT_DIR / "utils" / "plink")
if not os.path.isfile(PLINK_BIN):
    raise FileNotFoundError(f"Bundled PLINK binary not found: {PLINK_BIN}")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)




def read_plink_bed(input_prefix, max_snps=None):
    """Read PLINK bed file efficiently.
    
    Args:
        input_prefix: PLINK file prefix
        max_snps: Max SNPs to read (for memory efficiency)
    
    Returns:
        genotype_matrix: (n_snps x n_samples) with values 0,1,2 or -1 (missing)
        snps_df: DataFrame with SNP info
        fam_df: DataFrame with sample info
    """
    if os.path.isabs(input_prefix) or os.path.exists(input_prefix + '.bed'):
        bed_file = input_prefix + '.bed'
        bim_file = input_prefix + '.bim'
        fam_file = input_prefix + '.fam'
    else:
        bed_file = str(DATA_DIR / f"{input_prefix}.bed")
        bim_file = str(DATA_DIR / f"{input_prefix}.bim")
        fam_file = str(DATA_DIR / f"{input_prefix}.fam")
    
    # Read bim file
    snps = pd.read_csv(bim_file, sep="\t", header=None, 
                       names=["chr", "snp", "cm", "pos", "a1", "a2"])
    n_snps_total = len(snps)
    
    # Read fam file
    fam = pd.read_csv(fam_file, sep=r"\s+", header=None,
                      names=["fam", "id", "pat", "mat", "sex", "pheno"])
    n_samples = len(fam)
    
    snps_to_read = max_snps if max_snps and max_snps < n_snps_total else n_snps_total
    
    bytes_per_snp = (n_samples + 3) // 4
    # PLINK BED 2-bit lookup: 00->0(homA1), 01->-1(missing), 10->1(het), 11->2(homA2)
    _GENO_LUT = np.array([0, -1, 1, 2], dtype=np.int8)

    with open(bed_file, 'rb') as f:
        magic = f.read(3)
        if magic != b'\x6c\x1b\x01':
            raise ValueError("Not a valid PLINK bed file")
        raw = np.frombuffer(f.read(bytes_per_snp * snps_to_read), dtype=np.uint8)

    raw = raw.reshape(snps_to_read, bytes_per_snp)
    unpacked = np.empty((snps_to_read, bytes_per_snp * 4), dtype=np.int8)
    for shift in range(4):
        unpacked[:, shift::4] = _GENO_LUT[(raw >> (shift * 2)) & 0x03]
    genotype_data = unpacked[:, :n_samples]

    return genotype_data, snps.iloc[:snps_to_read], fam

# ========== Genotype Format Conversion ==========

def convert_hapmap_genotype(g, n):
    """Convert hapmap genotype to PLINK encoding."""
    c = [1, 4, 16, 64]
    s = 0
    for i in range(len(g)):
        if g[i][0] != g[i][1]:
            s += 2 * c[i]
        elif g[i] == 'NN':
            s += 1 * c[i]
        elif g[i][0] == n[0]:
            s += 0 * c[i]
        else:
            s += 3 * c[i]
    return s


def vcf_to_plink(vcf_file, output_prefix):
    """Convert VCF to PLINK format using PLINK.
    
    Args:
        vcf_file: Input VCF file
        output_prefix: Output PLINK prefix
    
    Returns:
        True if successful
    """
    try:
        cmd = f"{PLINK_BIN} --vcf {vcf_file} --out {output_prefix} --double-id --biallelic-only strict"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False


def plink_to_vcf(bed_prefix, output_prefix):
    """Convert PLINK to VCF format.
    
    Args:
        bed_prefix: Input PLINK prefix
        output_prefix: Output VCF prefix
    
    Returns:
        True if successful
    """
    try:
        cmd = f"{PLINK_BIN} --bfile {bed_prefix} --recode vcf-iid --out {output_prefix}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False


def hapmap_to_plink(hapmap_file, output_prefix):
    """Convert HapMap format to PLINK.
    
    Args:
        hapmap_file: Input HapMap file
        output_prefix: Output PLINK prefix
    
    Returns:
        True if successful
    """
    try:
        with open(hapmap_file) as h, open(output_prefix + '.bed', 'wb') as b, \
             open(output_prefix + '.bim', 'w') as bim, open(output_prefix + '.fam', 'w') as fam:
            
            # Write magic bytes for PLINK bed file (SNP-major mode)
            b.write(struct.pack('B', 108))
            b.write(struct.pack('B', 27))
            b.write(struct.pack('B', 1))
            
            out = list()
            for idx, l in enumerate(h):
                l = l.strip().split('\t')
                if idx == 0:
                    # Header line - sample IDs
                    samples = l[11:]
                    for s in samples:
                        fam.write(' '.join([s, s, '0', '0', '0', '-9']) + '\n')
                else:
                    # Data lines - SNP genotypes
                    g_seq = ''.join(l[11:])
                    nlu_num = list(set(g_seq))
                    
                    # Skip multi-allelic SNPs
                    if len(nlu_num) > 2 and 'N' not in nlu_num:
                        print(f"Warning: More than two nucleotide letters at SNP {l[2]}, skipping")
                        continue
                    if len(nlu_num) == 1:
                        print(f"Warning: Only one nucleotide letter at SNP {l[2]}, skipping")
                        continue
                    
                    # Determine alleles
                    if 'N' in nlu_num:
                        c = Counter(g_seq)
                        del c['N']
                        n = sorted(c, key=lambda x: c[x])
                    else:
                        g_seq_sort = ''.join(sorted(g_seq))
                        major = g_seq_sort[len(g_seq) // 2]
                        if g_seq_sort[0] == major:
                            n = [g_seq_sort[-1], major]
                        else:
                            n = [g_seq_sort[0], major]
                    
                    # Write BIM
                    bim.write('\t'.join([l[2], l[0], '0', l[3]] + n) + '\n')
                    
                    # Convert genotypes
                    for num in range(11, len(l), 4):
                        if num + 4 < len(l):
                            out.append(convert_hapmap_genotype(l[num:num + 4], n))
                        else:
                            out.append(convert_hapmap_genotype(l[num:len(l)], n))
            
            # Write genotype data
            b.write(struct.pack('B' * len(out), *out))
        
        return True
    except Exception as e:
        print(f"Error in hapmap_to_plink: {e}")
        return False


# ========== Genotype QC ==========

def subset_plink(input_prefix, output_prefix, chromosomes=None, proportion=None, seed=42):
    """Subset a PLINK dataset by chromosomes or random variant proportion.

    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix
        chromosomes: Chromosome string/list, e.g. "1,2,3" or ["1", "2"]
        proportion: Random SNP retention proportion in (0, 1]
        seed: Random seed used with proportion-based subset

    Returns:
        True if successful
    """
    has_chromosomes = chromosomes is not None and str(chromosomes).strip() != ""
    has_proportion = proportion is not None

    if has_chromosomes == has_proportion:
        raise ValueError("Specify exactly one of 'chromosomes' or 'proportion'.")

    cmd = [PLINK_BIN, "--bfile", input_prefix, "--make-bed", "--out", output_prefix]

    if has_chromosomes:
        if isinstance(chromosomes, (list, tuple, set)):
            chrom_tokens = [str(c).strip() for c in chromosomes if str(c).strip()]
        else:
            chrom_tokens = [c for c in str(chromosomes).replace(",", " ").split() if c]
        if not chrom_tokens:
            raise ValueError("No valid chromosome values were provided.")
        cmd.extend(["--chr", *chrom_tokens])
    else:
        proportion = float(proportion)
        if not (0 < proportion <= 1):
            raise ValueError("'proportion' must be within (0, 1].")
        cmd.extend(["--thin", str(proportion), "--seed", str(int(seed))])

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

def snp_qc(input_prefix, output_prefix, maf=0.05, missing_rate=0.2, mind=0.2):
    """SNP quality control using PLINK.
    
    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix
        maf: Minor allele frequency threshold
        missing_rate: Genotype missing rate threshold
        mind: Sample missing rate threshold
    
    Returns:
        True if successful
    """
    try:
        cmd = f"{PLINK_BIN} --bfile {input_prefix} --out {output_prefix} " \
              f"--maf {maf} --geno {missing_rate} --mind {mind} --make-bed"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False


def snp_impute(input_prefix, output_prefix, method='mean'):
    """Impute missing genotypes.
    
    Note: Full imputation requires reference panels. This is a simple mean imputation.
    
    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix  
        method: Imputation method ('mean' or 'median')
    
    Returns:
        True if successful
    """
    # Read genotype data
    try:
        
        G, snps, fam = read_plink_bed(input_prefix)
        
        G_imputed = G.astype(np.float32)
        for i in range(G_imputed.shape[0]):
            missing = G_imputed[i, :] == -1
            if missing.any() and (~missing).any():
                valid_vals = G_imputed[i, ~missing]
                fill = valid_vals.mean() if method == 'mean' else np.median(valid_vals)
                G_imputed[i, missing] = round(fill)
        
        # Save as PLINK bed file
        n_samples = G_imputed.shape[1]
        n_snps = G_imputed.shape[0]
        
        # Write .bed file — vectorized
        _ENC_LUT = np.zeros(256, dtype=np.uint8)
        _ENC_LUT[0] = 0    # hom A1 -> 00
        _ENC_LUT[1] = 2    # het    -> 10
        _ENC_LUT[2] = 3    # hom A2 -> 11
        # -1 (255 as uint8) -> missing -> 01
        _ENC_LUT[255] = 1

        G_write = G_imputed.astype(np.int8).view(np.uint8)
        encoded = _ENC_LUT[G_write]

        bytes_per_snp = (n_samples + 3) // 4
        pad = bytes_per_snp * 4 - n_samples
        if pad > 0:
            encoded = np.hstack([encoded, np.zeros((n_snps, pad), dtype=np.uint8)])

        packed = (encoded[:, 0::4]
                  | (encoded[:, 1::4] << 2)
                  | (encoded[:, 2::4] << 4)
                  | (encoded[:, 3::4] << 6)).astype(np.uint8)

        with open(output_prefix + '.bed', 'wb') as f:
            f.write(b'\x6c\x1b\x01')
            f.write(packed.tobytes())
        
        # Copy .bim and .fam
        shutil.copy(input_prefix + '.bim', output_prefix + '.bim')
        shutil.copy(input_prefix + '.fam', output_prefix + '.fam')
        
        return True
    except Exception as e:
        print(f"Error in snp_impute: {e}")
        return False


# ========== PCA and Relatedness ==========

def calculate_pca(input_prefix, n_components=10, max_snps=20000):
    """Calculate PCA from genotype data.
    
    Args:
        input_prefix: PLINK genotype prefix
        n_components: Number of principal components
        max_snps: Maximum SNPs to use for PCA
    
    Returns:
        DataFrame with PC scores
    """
    
    
    # Read genotype data
    G, snps, fam = read_plink_bed(input_prefix, max_snps=max_snps)
    n_samples = G.shape[1]
    
    G_clean = G.astype(np.float32)
    for i in range(G_clean.shape[0]):
        missing = G_clean[i, :] == -1
        if missing.any() and (~missing).any():
            G_clean[i, missing] = G_clean[i, ~missing].mean()
    
    # PCA
    X = G_clean.T  # (n_samples x n_snps)
    pca = PCA(n_components=min(n_components, n_samples - 1))
    pcs = pca.fit_transform(X)
    
    # Create output DataFrame
    pc_df = pd.DataFrame(pcs[:, :n_components], columns=[f"PC{i+1}" for i in range(n_components)])
    pc_df.insert(0, "sample", fam["id"].values)
    
    return pc_df, pca.explained_variance_ratio_[:n_components]


def calculate_ibd(input_prefix, output_prefix):
    """Calculate Identity-By-Descent (IBD) using PLINK.
    
    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output prefix
    
    Returns:
        True if successful
    """
    try:
        cmd = f"{PLINK_BIN} --bfile {input_prefix} --genome --out {output_prefix}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False


def calculate_kinship(input_prefix, output_prefix):
    """Calculate kinship/relatedness matrix.
    
    Args:
        input_prefix: PLINK genotype prefix
        output_prefix: Output file prefix
    
    Returns:
        Kinship matrix file path
    """
    output_prefix = os.path.abspath(output_prefix)
    output_dir = os.path.dirname(output_prefix) or os.getcwd()
    output_stem = os.path.basename(output_prefix)
    os.makedirs(output_dir, exist_ok=True)

    gemma_bin = str(SCRIPT_DIR / "utils" / "gemma.linux")
    cmd = [gemma_bin, "-bfile", input_prefix, "-gk", "1", "-outdir", output_dir, "-o", output_stem]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode == 0:
        return os.path.join(output_dir, f"{output_stem}.cXX.txt")
    return None


# ========== Utility Functions ==========

def read_plink(input_prefix):
    """Read PLINK files into DataFrames.
    
    Args:
        input_prefix: PLINK file prefix
    
    Returns:
        Dictionary with 'bim', 'fam', 'geno' (if loaded)
    """
    if os.path.isabs(input_prefix) or os.path.exists(input_prefix + '.bim'):
        bim_file = input_prefix + '.bim'
        fam_file = input_prefix + '.fam'
    else:
        bim_file = str(DATA_DIR / f"{input_prefix}.bim")
        fam_file = str(DATA_DIR / f"{input_prefix}.fam")
    bim = pd.read_csv(bim_file, sep='\t', header=None,
                     names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                     names=['fam', 'id', 'pat', 'mat', 'sex', 'pheno'])
    
    return {'bim': bim, 'fam': fam}


def get_snp_stats(input_prefix):
    """Get basic SNP statistics.
    
    Args:
        input_prefix: PLINK genotype prefix
    
    Returns:
        Dictionary with statistics
    """
    bim, fam = read_plink(input_prefix).values()
    
    return {
        'n_snps': len(bim),
        'n_samples': len(fam),
        'chromosomes': sorted(bim['chr'].unique().tolist())
    }




# ========== SNP Pruning ==========

def snp_pruning(input_prefix, output_prefix, window=50, shift=5, r2=0.5,maf=0.05):
    """LD-based SNP pruning using PLINK.
    
    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix
        window: Window size (kb)
        shift: Shift step
        r2: R2 threshold for pruning
        maf: Minor allele frequency threshold
    
    Returns:
        True if successful
    """
    try:
        cmd = f"{PLINK_BIN} --bfile {input_prefix} --out {output_prefix} "               f"--indep-pairwise {window} {shift} {r2} --maf {maf}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # Create pruned dataset
        prune_file = f"{output_prefix}.prune.in"
        if os.path.exists(prune_file):
            cmd2 = f"{PLINK_BIN} --bfile {input_prefix} --extract {prune_file} "                    f"--make-bed --out {output_prefix}_pruned"
            result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
            return result2.returncode == 0
        
        return result.returncode == 0
    except Exception:
        return False


def snp_clumping(input_prefix, output_prefix, r2=0.5, maf=0.05, window_kb=250):
    """Genotype-based LD clumping (pure Python, equivalent to bigsnpr::snp_clumping).

    Greedy algorithm per chromosome: iterate through SNPs sorted by position,
    keep a SNP only if it is not in high LD (r2 > threshold) with any
    already-selected SNP within a distance window.

    Args:
        input_prefix: Input PLINK prefix
        output_prefix: Output PLINK prefix
        r2: r-squared threshold for LD
        maf: Minor allele frequency filter
        window_kb: LD window in kilobases

    Returns:
        True if successful
    """
    try:
        snps = pd.read_csv(input_prefix + '.bim', sep='\t', header=None,
                           names=['chr', 'snp', 'cm', 'pos', 'a1', 'a2'])
        n_snps_total = len(snps)
        window_bp = window_kb * 1000

        G, _, _ = read_plink_bed(input_prefix)
        # G is (n_snps, n_samples) — transpose to (n_samples, n_snps)
        G = G.T.astype(np.float32)
        n_samples = G.shape[0]
        G[G == -1] = np.nan

        af = np.nanmean(G, axis=0) / 2.0
        maf_ok = np.minimum(af, 1.0 - af) >= maf

        # Standardize each column (mean 0, unit norm); NaN → 0
        col_mean = np.nanmean(G, axis=0)
        col_std = np.nanstd(G, axis=0)
        col_std[col_std == 0] = 1.0
        nan_mask = np.isnan(G)
        G = (G - col_mean) / col_std
        G[nan_mask] = 0.0

        norms = np.linalg.norm(G, axis=0)
        norms[norms == 0] = 1.0

        keep_global = []

        for chrom in sorted(snps['chr'].unique()):
            chr_idx = np.where((snps['chr'].values == chrom) & maf_ok)[0]
            if len(chr_idx) == 0:
                continue
            pos = snps['pos'].values[chr_idx]
            order = np.argsort(pos)
            chr_idx = chr_idx[order]
            pos = pos[order]

            sel_idx = []
            sel_pos = []

            for k, gi in enumerate(chr_idx):
                if not sel_idx:
                    sel_idx.append(gi)
                    sel_pos.append(pos[k])
                    continue

                p = pos[k]
                in_ld = False
                for si, sp in zip(reversed(sel_idx), reversed(sel_pos)):
                    if p - sp > window_bp:
                        break
                    dot = np.dot(G[:, gi], G[:, si])
                    r2_val = (dot / (norms[gi] * norms[si])) ** 2
                    if r2_val > r2:
                        in_ld = True
                        break
                if not in_ld:
                    sel_idx.append(gi)
                    sel_pos.append(p)

            keep_global.extend(sel_idx)

        print(f"[info] snp_clumping: {n_snps_total} -> {len(keep_global)} SNPs "
              f"(r2={r2}, maf={maf}, window={window_kb}kb)")

        keep_file = f"{output_prefix}.clump.in"
        keep_names = snps['snp'].values[keep_global]
        with open(keep_file, 'w') as f:
            for s in keep_names:
                f.write(s + '\n')

        cmd = (f"{PLINK_BIN} --bfile {input_prefix} --extract {keep_file} "
               f"--make-bed --out {output_prefix}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"[error] snp_clumping plink: {result.stderr[:500] if result.stderr else 'no stderr'}")
        return result.returncode == 0
    except Exception as e:
        print(f"[error] snp_clumping: {e}")
        return False




# ========== Tree Building ==========

def ped_to_fasta(ped_file, output_prefix):
    """Convert PED file to FASTA format.
    
    Args:
        ped_file: Input PED file
        output_prefix: Output FASTA prefix
    
    Returns:
        FASTA file path
    """
    fasta_file = f"{output_prefix}.fasta"
    
    with open(ped_file, 'r') as fi, open(fasta_file, 'w') as fo:
        for line in fi:
            tmp = line.strip().split()
            if len(tmp) >= 7:
                name = tmp[1]
                seq = ''.join(tmp[6:])
                seq = seq.replace('0', 'N')
                fo.write(f">{name}\n{seq}\n")
    
    return fasta_file


def build_phylogenetic_tree(bed_prefix, output_prefix):
    """Build phylogenetic tree using FastTree.
    
    Args:
        bed_prefix: PLINK genotype prefix
        output_prefix: Output prefix
    
    Returns:
        Tree file path (.nwk)
    """
    import re
    
    # Convert to PED
    ped_prefix = f"{output_prefix}_tree"
    cmd1 = f"{PLINK_BIN} --bfile {bed_prefix} --recode --out {ped_prefix}"
    result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    
    if result1.returncode != 0:
        return None
    
    # Convert PED to FASTA
    fasta_file = ped_to_fasta(f"{ped_prefix}.ped", f"{output_prefix}")
    
    # Run FastTree
    tree_file = f"{output_prefix}.tree.nwk"
    fasttree_bin = str(SCRIPT_DIR / "utils" / "FastTree")
    cmd2 = f"{fasttree_bin} -nt -gtr -quiet {fasta_file} > {tree_file}"
    result2 = subprocess.run(cmd2, shell=True, stderr=subprocess.PIPE, text=True)
    
    if result2.returncode != 0:
        return None
    
    # Cleanup intermediate files
    for ext in ['.ped', '.map', '.nof']:
        f = f"{ped_prefix}{ext}"
        if os.path.exists(f):
            os.remove(f)
    if os.path.exists(fasta_file):
        os.remove(fasta_file)
    
    return tree_file


__all__ = [
    'vcf_to_plink', 'plink_to_vcf', 'hapmap_to_plink',
    'subset_plink', 'snp_qc', 'snp_impute',
    'calculate_pca', 'calculate_ibd', 'calculate_kinship',
    'read_plink', 'get_snp_stats', 'snp_pruning', 'snp_clumping', 'ped_to_fasta', 'build_phylogenetic_tree'
]
