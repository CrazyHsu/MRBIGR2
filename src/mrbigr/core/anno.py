#!/usr/bin/env python3
"""
AnnoMCP - Annotation Module
Supports both GTF and ANNOVAR annotation
"""
import pandas as pd
import numpy as np
import os
import subprocess
import warnings

warnings.filterwarnings("ignore")

try:
    import pyranges as pr
    HAS_PYRANGES = True
except ImportError:
    HAS_PYRANGES = False


def _normalize_chr_value(value):
    if pd.isna(value):
        return None
    chrom = str(value).strip()
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    return chrom


def _extract_attr(attr_str, key):
    if pd.isna(attr_str):
        return None
    for item in str(attr_str).split(';'):
        item = item.strip()
        if not item:
            continue
        if item.startswith(f'{key} '):
            parts = item.split('"')
            if len(parts) >= 2:
                return parts[1]
            return item.split(None, 1)[1].strip('"')
        if item.startswith(f'{key}='):
            return item.split('=', 1)[1].strip('"')
    return None


def _load_gene_annotation(annotation_file):
    if annotation_file is None:
        return None
    if isinstance(annotation_file, pd.DataFrame):
        genes = annotation_file.copy()
    else:
        annotation_path = str(annotation_file)
        lower_path = annotation_path.lower()
        if lower_path.endswith(('.gtf', '.gff', '.gtf.gz', '.gff.gz')):
            genes = get_genes_from_gtf(annotation_path)
        elif lower_path.endswith('.pkl'):
            genes = pd.read_pickle(annotation_path)
        else:
            genes = pd.read_csv(annotation_path)
    if genes is None:
        return None
    genes = genes.copy()
    if 'feature' in genes.columns:
        genes = genes[genes['feature'].astype(str).str.lower() == 'gene'].copy()
    has_gene_id = 'gene_id' in genes.columns
    genes['chr'] = genes['chr'].map(_normalize_chr_value)
    genes['start'] = pd.to_numeric(genes['start'], errors='raise').astype(int)
    genes['end'] = pd.to_numeric(genes['end'], errors='raise').astype(int)
    if 'gene_id' not in genes.columns:
        genes['gene_id'] = None
    if 'gene_name' not in genes.columns:
        genes['gene_name'] = None
    if 'strand' not in genes.columns:
        genes['strand'] = None
    if has_gene_id:
        genes = genes.dropna(subset=['gene_id']).drop_duplicates(subset=['gene_id'])
    return genes


def _coerce_snp_df(snps):
    from ._argjson import maybe_json_loads
    snps = maybe_json_loads(snps)
    if isinstance(snps, pd.DataFrame):
        snp_df = snps.copy()
    elif isinstance(snps, dict):
        if any(isinstance(v, (list, tuple, pd.Series, np.ndarray)) for v in snps.values()):
            snp_df = pd.DataFrame(snps)
        else:
            snp_df = pd.DataFrame([snps])
    elif isinstance(snps, list):
        if not snps:
            snp_df = pd.DataFrame(columns=['chr', 'pos'])
        elif all(isinstance(s, str) for s in snps):
            snp_df = pd.DataFrame([s.split(':', 1) for s in snps], columns=['chr', 'pos'])
        else:
            snp_df = pd.DataFrame(snps)
    else:
        raise TypeError("snps must be a DataFrame, dict, list[dict], or list['chr:pos']")
    if 'pos' not in snp_df.columns and 'ps' in snp_df.columns:
        snp_df = snp_df.rename(columns={'ps': 'pos'})
    if 'chr' not in snp_df.columns or 'pos' not in snp_df.columns:
        raise ValueError("SNP input must contain 'chr' and 'pos' columns")
    snp_df['chr'] = snp_df['chr'].map(_normalize_chr_value)
    snp_df['pos'] = pd.to_numeric(snp_df['pos'], errors='raise').astype(int)
    return snp_df


def parse_gtf(gtf_file):
    """Parse GTF annotation file."""
    if not os.path.exists(gtf_file):
        return None
    try:
        df = pd.read_csv(
            gtf_file,
            sep='\t',
            header=None,
            names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'],
            dtype={'chr': str, 'start': int, 'end': int},
            comment='#',
            compression='infer',
        )
        df['chr'] = df['chr'].map(_normalize_chr_value)
        df['gene_id'] = df['attributes'].apply(lambda x: _extract_attr(x, 'gene_id'))
        df['gene_name'] = df['attributes'].apply(lambda x: _extract_attr(x, 'gene_name'))
        df['transcript_id'] = df['attributes'].apply(lambda x: _extract_attr(x, 'transcript_id'))
        df['gene_type'] = df['attributes'].apply(
            lambda x: _extract_attr(x, 'gene_type') or _extract_attr(x, 'gene_biotype')
        )
        return df
    except Exception as e:
        print(f"Error parsing GTF: {e}")
        return None


def get_genes_from_gtf(gtf_file):
    """Extract gene list from GTF."""
    df = parse_gtf(gtf_file)
    if df is None:
        return None
    genes = df[df['feature'] == 'gene'][['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name', 'gene_type']].copy()
    genes = genes.dropna(subset=['gene_id']).drop_duplicates(subset=['gene_id'])
    return genes


def annotate_snps_simple(snps, gtf_file, window=5000):
    """Annotate SNPs using GTF."""
    genes = _load_gene_annotation(gtf_file)
    if genes is None:
        return None
    genes = genes[['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name']].copy()
    snp_df = _coerce_snp_df(snps)
    results = []
    for _, snp in snp_df.iterrows():
        chr_val, pos = snp['chr'], int(snp['pos'])
        overlapping = genes[(genes['chr'] == chr_val) & (genes['start'] <= pos + window) & (genes['end'] >= pos - window)]
        if len(overlapping) > 0:
            for _, gene in overlapping.iterrows():
                if gene['start'] <= pos <= gene['end']:
                    region, dist = "exon", 0
                else:
                    dist = min(abs(pos - gene['start']), abs(pos - gene['end']))
                    region = "upstream" if pos < gene['start'] else "downstream"
                gene_label = gene['gene_name'] if pd.notna(gene['gene_name']) else gene['gene_id']
                results.append({
                    'chr': chr_val,
                    'pos': pos,
                    'gene_id': gene['gene_id'],
                    'gene_name': gene['gene_name'],
                    'distance': dist,
                    'region': region,
                    'annotation': f"{gene_label}:{region}" if pd.notna(gene_label) else region
                })
        else:
            results.append({'chr': chr_val, 'pos': pos, 'gene_id': None, 'gene_name': None, 'distance': None, 'region': 'intergenic', 'annotation': 'intergenic'})
    return pd.DataFrame(results)


def annotate_vcf_simple(vcf_file, gtf_file, output_file=None):
    """Annotate VCF using GTF."""
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variants.append({'chr': parts[0], 'pos': int(parts[1]), 'id': parts[2] if parts[2] != '.' else None, 'ref': parts[3], 'alt': parts[4]})
    vcf_df = pd.DataFrame(variants)
    annotated = annotate_snps_simple(vcf_df, gtf_file)
    if annotated is not None and output_file:
        annotated.to_csv(output_file, sep='\t', index=False)
    return annotated


def get_genes_in_region(chr_name, start, end, annotation_file, mode='overlap'):
    """Get genes in genomic region."""
    genes = _load_gene_annotation(annotation_file)
    if genes is None:
        return None
    chr_name = _normalize_chr_value(chr_name)
    genes = genes[genes['chr'] == chr_name]
    if mode == 'contain':
        return genes[(genes['start'] >= start) & (genes['end'] <= end)]
    if mode != 'overlap':
        raise ValueError("mode must be 'overlap' or 'contain'")
    return genes[(genes['start'] <= end) & (genes['end'] >= start)]


def query_gene_position(gene_name, gtf_file):
    """Get gene position from GTF."""
    genes = get_genes_from_gtf(gtf_file)
    if genes is None:
        return None
    return genes[(genes['gene_name'] == gene_name) | (genes['gene_id'] == gene_name)]


def calculate_tss_distance(snp_chr, snp_pos, gtf_file):
    """Calculate distance to nearest TSS."""
    genes = _load_gene_annotation(gtf_file)
    if genes is None:
        return None
    genes = genes[genes['chr'] == _normalize_chr_value(snp_chr)]
    distances = []
    for _, gene in genes.iterrows():
        tss = gene['start'] if gene['strand'] == '+' else gene['end']
        distances.append({'gene_id': gene['gene_id'], 'gene_name': gene['gene_name'], 'strand': gene['strand'], 'tss': tss, 'distance': abs(snp_pos - tss)})
    return pd.DataFrame(distances).sort_values('distance')


def predict_variant_effect(snp_chr=None, snp_pos=None, ref_allele=None, alt_allele=None, gtf_file=None, 
                          simple_ref=None, simple_alt=None, region_type=None):
    """Predict variant effect.
    
    Full version:
        predict_variant_effect(snp_chr, snp_pos, ref_allele, alt_allele, gtf_file)
    Simple version (for quick tests):
        predict_variant_effect(simple_ref='A', simple_alt='T', region_type='exon')
    """
    # Simple version for quick testing
    if simple_ref is not None and simple_alt is not None:
        ref = str(simple_ref).upper()
        alt = str(simple_alt).upper()
        
        # Simple effect prediction based on allele change
        if region_type == 'exon':
            if len(ref) == len(alt):
                if ref != alt:
                    effect = 'missense' if len(ref) == 1 else 'indel'
                else:
                    effect = 'synonymous'
            else:
                effect = 'frameshift' if abs(len(ref) - len(alt)) % 3 != 0 else 'inframe_indel'
        elif region_type == 'intron':
            effect = 'intronic'
        elif region_type == 'promoter':
            effect = 'regulatory'
        elif region_type == 'utr':
            effect = 'UTR'
        else:
            effect = 'intergenic'
        
        return {'effect': effect, 'ref': ref, 'alt': alt, 'region': region_type}
    
    # Full version requiring GTF
    if gtf_file is None:
        return {'effect': 'unknown', 'gene': None, 'error': 'gtf_file required for full annotation'}

    chr_name = _normalize_chr_value(snp_chr)
    lower_path = str(gtf_file).lower()
    if lower_path.endswith(('.gtf', '.gff', '.gtf.gz', '.gff.gz')):
        gtf_df = parse_gtf(gtf_file)
        if gtf_df is None:
            return {'effect': 'unknown', 'gene': None}
        gene_rows = gtf_df[(gtf_df['feature'] == 'gene') & (gtf_df['chr'] == chr_name)]
        gene_hits = gene_rows[(gene_rows['start'] <= snp_pos) & (gene_rows['end'] >= snp_pos)]
        if len(gene_hits) == 0:
            return {'effect': 'intergenic', 'gene': None}
        gene = gene_hits.iloc[0]
        feature_hits = gtf_df[
            (gtf_df['chr'] == chr_name) &
            (gtf_df['gene_id'] == gene['gene_id']) &
            (gtf_df['start'] <= snp_pos) &
            (gtf_df['end'] >= snp_pos)
        ]
        features = set(feature_hits['feature'].dropna().tolist())
        if 'CDS' in features:
            effect = 'coding'
        elif 'exon' in features:
            effect = 'exonic'
        elif {'five_prime_utr', 'three_prime_utr', 'UTR'} & features:
            effect = 'UTR'
        else:
            effect = 'intronic'
    else:
        genes = _load_gene_annotation(gtf_file)
        if genes is None:
            return {'effect': 'unknown', 'gene': None}
        gene_hits = genes[(genes['chr'] == chr_name) & (genes['start'] <= snp_pos) & (genes['end'] >= snp_pos)]
        if len(gene_hits) == 0:
            return {'effect': 'intergenic', 'gene': None}
        gene = gene_hits.iloc[0]
        effect = 'genic'
    return {
        'effect': effect,
        'gene_id': gene['gene_id'],
        'gene_name': gene['gene_name'],
        'chr': chr_name,
        'pos': snp_pos,
        'ref': ref_allele,
        'alt': alt_allele
    }


def create_annotation_db(gtf_file, output_prefix):
    """Create annotation database from GTF."""
    genes = get_genes_from_gtf(gtf_file)
    if genes is None:
        return None
    db_file = f"{output_prefix}_db.pkl"
    genes.to_pickle(db_file)
    return db_file


def qtl_annotation(qtl_file, annotation_file):
    """Annotate QTL regions with genes."""
    qtl_df = pd.read_csv(qtl_file) if not isinstance(qtl_file, pd.DataFrame) else qtl_file.copy()
    genes = _load_gene_annotation(annotation_file)
    if genes is None:
        return None
    chr_col = 'chr' if 'chr' in qtl_df.columns else 'CHR'
    has_range = 'qtl_start' in qtl_df.columns and 'qtl_end' in qtl_df.columns
    if not has_range:
        pos_col = 'pos' if 'pos' in qtl_df.columns else ('ps' if 'ps' in qtl_df.columns else 'position')
    qtl_df = qtl_df.copy()
    qtl_df[chr_col] = qtl_df[chr_col].map(_normalize_chr_value)
    results = []
    for _, qtl in qtl_df.iterrows():
        chr_val = qtl[chr_col]
        if has_range:
            q_start, q_end = int(qtl['qtl_start']), int(qtl['qtl_end'])
            overlap = genes[(genes['chr'] == chr_val) & (genes['end'] >= q_start) & (genes['start'] <= q_end)]
        else:
            pos = int(qtl[pos_col])
            overlap = genes[(genes['chr'] == chr_val) & (genes['start'] <= pos) & (genes['end'] >= pos)]
        if len(overlap) > 0:
            for _, gene in overlap.iterrows():
                results.append({**qtl.to_dict(), 'gene_id': gene['gene_id'], 'gene_name': gene['gene_name']})
        else:
            results.append({**qtl.to_dict(), 'gene_id': None, 'gene_name': None})
    return pd.DataFrame(results)


def find_annovar():
    """Find ANNOVAR installation path."""
    bundled_utils = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "utils"
    )
    possible_paths = [
        bundled_utils,
        '/opt/annovar',
        os.path.expanduser('~/annovar'),
    ]
    for path in possible_paths:
        if os.path.exists(path):
            annotate_variation = os.path.join(path, 'annotate_variation.pl')
            if os.path.exists(annotate_variation):
                return path
    annovar_dir = os.environ.get('ANNOVAR_DIR')
    if annovar_dir and os.path.exists(annovar_dir):
        return annovar_dir
    return None


def annotate_with_annovar(vcf_file, database_dir, output_prefix, buildver='hg38'):
    """Annotate VCF file using ANNOVAR."""
    annovar_path = find_annovar()
    if annovar_path is None:
        print("Warning: ANNOVAR not found. Please install ANNOVAR or set ANNOVAR_DIR.")
        return False
    table_annovar = os.path.join(annovar_path, 'table_annovar.pl')
    if not os.path.exists(table_annovar):
        print(f"Error: table_annovar.pl not found at {table_annovar}")
        return False
    os.makedirs(output_prefix, exist_ok=True)
    cmd = f"perl {table_annovar} {vcf_file} {database_dir} -buildver {buildver} -out {output_prefix}/annovar_out -remove -protocol refGene,cytoBand,exac03,gnomad30,dbnsfp42 -operation g,r,f,f -nastring . -vcfinput"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in ANNOVAR annotation: {result.stderr}")
        return False
    return True


def convert_to_annovar(input_file, output_file, format='avinput'):
    """Convert formats to ANNOVAR input."""
    annovar_path = find_annovar()
    if annovar_path is None:
        print("Warning: ANNOVAR not found.")
        return False
    convert2annovar = os.path.join(annovar_path, 'convert2annovar.pl')
    if not os.path.exists(convert2annovar):
        print(f"Error: convert2annovar.pl not found")
        return False
    if format == 'vcf':
        cmd = f"perl {convert2annovar} {input_file} -format vcf > {output_file}"
    else:
        cmd = f"perl {convert2annovar} {input_file} > {output_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in format conversion: {result.stderr}")
        return False
    return True


def gene_based_annotation(avinput_file, database_dir, output_prefix, buildver='hg38'):
    """Gene-based annotation using ANNOVAR."""
    annovar_path = find_annovar()
    if annovar_path is None:
        print("Warning: ANNOVAR not found.")
        return False
    annotate_variation = os.path.join(annovar_path, 'annotate_variation.pl')
    if not os.path.exists(annotate_variation):
        print(f"Error: annotate_variation.pl not found")
        return False
    os.makedirs(output_prefix, exist_ok=True)
    cmd = f"perl {annotate_variation} -gene {avinput_file} {database_dir} -buildver {buildver} -out {output_prefix}/gene_anno"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in gene-based annotation: {result.stderr}")
        return False
    return True


def filter_annovar_variants(input_file, database_dir, output_file, buildver='hg38', filter_type='gene'):
    """Filter variants using ANNOVAR."""
    annovar_path = find_annovar()
    if annovar_path is None:
        print("Warning: ANNOVAR not found.")
        return False
    annotate_variation = os.path.join(annovar_path, 'annotate_variation.pl')
    if filter_type == 'gene':
        cmd = f"perl {annotate_variation} -filter {input_file} {database_dir} -buildver {buildver} -out {output_file} -gene"
    else:
        cmd = f"perl {annotate_variation} -filter {input_file} {database_dir} -buildver {buildver} -out {output_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error in filtering: {result.stderr}")
        return False
    return True


def read_annovar_output(annovar_out_file):
    """Read ANNOVAR table output."""
    if not os.path.exists(annovar_out_file):
        return None
    try:
        df = pd.read_csv(annovar_out_file, sep='\t')
        return df
    except Exception as e:
        print(f"Error reading ANNOVAR output: {e}")
        return None


__all__ = ['parse_gtf', 'get_genes_from_gtf', 'annotate_snps_simple', 
           'annotate_vcf_simple', 'get_genes_in_region', 'query_gene_position',
           'calculate_tss_distance', 'predict_variant_effect', 
           'create_annotation_db', 'qtl_annotation', 'HAS_PYRANGES',
           'find_annovar', 'annotate_with_annovar', 'convert_to_annovar',
           'gene_based_annotation', 'filter_annovar_variants', 'read_annovar_output']
