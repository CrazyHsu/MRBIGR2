#!/usr/bin/env python3
"""
GenotypeMCP Server - FastMCP server for genomic selection and GWAS
Pure Python implementation - no R dependencies

Uses modular structure:
- pheno.py: Phenotype analysis (filtering, scaling, imputation, QC)
- geno.py: Genotype processing (QC, PCA, IBD, format conversion)  
- gwas.py: GWAS analysis (LM, LMM, QTL detection)
- vis.py: Visualization (Manhattan, QQ, PCA, heatmaps)
- anno.py: Annotation (GTF parsing, VCF annotation, QTL annotation)
- qtl.py: QTL analysis (QTL detection, peak SNPs, gene mapping)
- multi.py: Parallel processing utilities
- peak.py: Peak test (boxplot, haplotype analysis)
"""
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

try:
    from fastmcp import FastMCP
    HAS_FASTMCP = True
except ImportError:
    HAS_FASTMCP = False
    print("Warning: fastmcp not installed")

# Import new modular components
try:
    import pheno
    import gwas
    import geno
    import vis
    import anno
    import qtl
    import multi
    import peak
    import mr
    import go
    import net
    HAS_MODULES = True
except ImportError as e:
    HAS_MODULES = False
    print(f"Warning: modules not available: {e}")

def run_mcp_server():
    if not HAS_FASTMCP:
        print("Error: fastmcp required")
        sys.exit(1)
    
    mcp = FastMCP("GenotypeMCP")

    # ========== Genotype Tools (from geno.py) ==========
    
    @mcp.tool()
    def run_snp_qc(input_prefix, output_prefix, maf=0.05, missing_rate=0.2, mind=0.2):
        """SNP quality control using PLINK."""
        return geno.snp_qc(input_prefix, output_prefix, maf=maf, missing_rate=missing_rate, mind=mind)

    @mcp.tool()
    def subset_genotype(input_prefix, output_prefix, chromosomes=None, proportion=None, seed=42):
        """Subset PLINK genotype data by chromosome or random proportion."""
        return geno.subset_plink(input_prefix, output_prefix, chromosomes=chromosomes, proportion=proportion, seed=seed)

    @mcp.tool()
    def run_genotype_pca(input_prefix, n_components=10, max_snps=20000):
        """PCA analysis for genotype data."""
        pc_df, var_ratio = geno.calculate_pca(input_prefix, n_components=n_components, max_snps=max_snps)
        return {"pc_data": pc_df.to_dict(), "variance_explained": var_ratio.tolist()}

    @mcp.tool()
    def run_calculate_ibd(input_prefix, output_prefix):
        """Calculate Identity by Descent (IBD) matrix."""
        return geno.calculate_ibd(input_prefix, output_prefix)

    @mcp.tool()
    def convert_vcf(vcf_file, output_prefix):
        """Convert VCF to PLINK format."""
        return geno.vcf_to_plink(vcf_file, output_prefix)

    @mcp.tool()
    def convert_hapmap(hapmap_file, output_prefix):
        """Convert HapMap format to PLINK."""
        return geno.hapmap_to_plink(hapmap_file, output_prefix)

    @mcp.tool()
    def run_plink_to_vcf(bed_prefix, output_prefix):
        """Convert PLINK to VCF format."""
        return geno.plink_to_vcf(bed_prefix, output_prefix)

    @mcp.tool()
    def calculate_kinship(input_prefix, output_prefix):
        """Calculate kinship/relatedness matrix using GEMMA."""
        return geno.calculate_kinship(input_prefix, output_prefix)

    @mcp.tool()
    def impute_genotype(input_prefix, output_prefix, method='mean'):
        """Impute missing genotype values."""
        return geno.snp_impute(input_prefix, output_prefix, method=method)

    @mcp.tool()
    def run_snp_pruning(input_prefix, output_prefix, window=50, shift=5, r2=0.5, maf=0.05):
        """LD-based SNP pruning using PLINK."""
        return geno.snp_pruning(input_prefix, output_prefix, window=window, shift=shift, r2=r2, maf=maf)

    @mcp.tool()
    def run_snp_clumping(input_prefix, output_prefix, r2=0.5, maf=0.05, window_kb=250):
        """LD-based SNP clumping."""
        return geno.snp_clumping(input_prefix, output_prefix, r2=r2, maf=maf, window_kb=window_kb)

    @mcp.tool()
    def get_snp_statistics(input_prefix):
        """Get basic SNP statistics."""
        return geno.get_snp_stats(input_prefix)

    # ========== Phenotype Tools (from pheno.py) ==========
    
    @mcp.tool()
    def filter_abundance(d, abundance):
        """Filter features by minimum abundance threshold."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.abundance_filter(df, abundance).to_dict()

    @mcp.tool()
    def filter_missing(d, missing_ratio):
        """Filter features by missing ratio threshold."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.missing_filter(df, missing_ratio).to_dict()

    @mcp.tool()
    def scale_phenotype(d, method='zscore'):
        """Scale phenotype data (log2, log10, zscore, minmax, robust, boxcox)."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.scale_wrapper(df, method).to_dict()

    @mcp.tool()
    def impute_phenotype(d, method='mean'):
        """Impute missing phenotype values."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.pheno_imputer(df, method=method).to_dict()

    @mcp.tool()
    def remove_outliers(d, method='zscore'):
        """Remove outliers from phenotype data."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.outlier(df, method=method).to_dict()

    @mcp.tool()
    def run_blup(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUP via REML mixed model."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.blup(df, method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def run_blue(d, method='matrix', y_col='y', geno_col='genotype', env_col='env', rep_col=None):
        """BLUE via OLS fixed-effects model."""
        import pandas as pd
        df = pd.DataFrame(d)
        return pheno.blue(df, method=method, y_col=y_col, geno_col=geno_col, env_col=env_col, rep_col=rep_col).to_dict()

    @mcp.tool()
    def correct_trait(pc, y):
        """Correct phenotype using principal components."""
        import pandas as pd
        import numpy as np
        pc_df = pd.DataFrame(pc)
        y_arr = np.array(y)
        return pheno.trait_correct(pc_df, y_arr).tolist()

    # ========== GWAS Tools (from gwas.py) ==========
    
    def _auto_plot(result_files):
        """Generate Manhattan + QQ plots for GWAS result files."""
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
        phe_df = pd.DataFrame(phe)
        result_files = gwas.gwas_lm(phe_df, geno_prefix, output_name=output_name, output_dir=output_dir)
        if auto_plot and result_files:
            return {"gwas_results": result_files, "plots": _auto_plot(result_files)}
        return result_files

    @mcp.tool()
    def run_gwas_lmm(phe, geno_prefix, output_name=None, output_dir=None, auto_plot=False):
        """Linear Mixed Model GWAS using GEMMA. Set auto_plot=True to generate Manhattan and QQ plots."""
        import pandas as pd
        phe_df = pd.DataFrame(phe)
        result_files = gwas.gwas_lmm(phe_df, geno_prefix, output_name=output_name, output_dir=output_dir)
        if auto_plot and result_files:
            return {"gwas_results": result_files, "plots": _auto_plot(result_files)}
        return result_files

    @mcp.tool()
    def run_gwas_plink(phe, geno_prefix, pheno_col=None, output_name="gwas", auto_plot=False):
        """GWAS using PLINK linear regression. Set auto_plot=True to generate Manhattan and QQ plots."""
        import pandas as pd
        phe_df = pd.DataFrame(phe)
        result_file = gwas.gwas_plink(phe_df, geno_prefix, pheno_col=pheno_col, output_name=output_name)
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
        y_arr = np.array(Y)
        g_arr = np.array(G)
        cov_arr = np.array(cov) if cov is not None else None
        result = gwas.simple_gwas(y_arr, g_arr, cov=cov_arr)
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
        arr = np.array(pvalues)
        return gwas.calculate_lambda(arr)

    @mcp.tool()
    def qq_plot_data(pvalues):
        """Generate QQ plot data from p-values."""
        import numpy as np
        arr = np.array(pvalues)
        return gwas.qq_plot(arr)

    @mcp.tool()
    def manhattan_plot_data(gwas_results):
        """Generate Manhattan plot data from GWAS results."""
        import pandas as pd
        df = pd.DataFrame(gwas_results)
        return gwas.manhattan_plot(df).to_dict()


    # ========== Visualization Tools (from vis.py) ==========
    
    @mcp.tool()
    def plot_manhattan(gwas_file, output_file=None, significance=5e-8, suggest=1e-5):
        """Generate Manhattan plot from GWAS results."""
        return vis.manhattan_plot(gwas_file, output_file=output_file, 
                                   significance=significance, suggest=suggest)

    @mcp.tool()
    def plot_qq(gwas_file, output_file=None):
        """Generate Q-Q plot from GWAS results."""
        return vis.qq_plot(gwas_file, output_file=output_file)

    @mcp.tool()
    def plot_pca(pca_file, output_file=None, pc_x='PC1', pc_y='PC2'):
        """Generate PCA scatter plot."""
        return vis.pca_plot(pca_file, output_file=output_file, pc_x=pc_x, pc_y=pc_y)

    @mcp.tool()
    def plot_tsne(tsne_file, output_file=None):
        """Generate t-SNE scatter plot."""
        return vis.tsne_plot(tsne_file, output_file=output_file)

    @mcp.tool()
    def plot_ld_heatmap(geno_prefix, output_file=None, max_snps=500):
        """Generate LD heatmap from genotype data."""
        return vis.ld_heatmap(geno_prefix, output_file=output_file, max_snps=max_snps)

    @mcp.tool()
    def plot_phenotype_hist(pheno_file, output_file=None):
        """Generate histogram for phenotype distribution."""
        return vis.phenotype_hist(pheno_file, output_file=output_file)

    @mcp.tool()
    def plot_phenotype_boxplot(pheno_file, output_file=None):
        """Generate boxplot for phenotypes."""
        return vis.phenotype_boxplot(pheno_file, output_file=output_file)

    @mcp.tool()
    def plot_phenotype_correlation(pheno_file, output_file=None, method='pearson'):
        """Generate phenotype correlation heatmap."""
        return vis.phenotype_correlation(pheno_file, output_file=output_file, method=method)

    @mcp.tool()
    def gwas_summary(gwas_dir, output_prefix='gwas_summary'):
        """Generate comprehensive GWAS summary plots (Manhattan + QQ)."""
        return vis.gwas_summary_plot(gwas_dir, output_prefix=output_prefix)


    # ========== Annotation Tools (from anno.py) ==========
    
    @mcp.tool()
    def parse_gtf(gtf_file):
        """Parse GTF annotation file."""
        result = anno.parse_gtf(gtf_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def annotate_snps(snp_df, gtf_file, window=5000):
        """Annotate SNPs using GTF file."""
        result = anno.annotate_snps_simple(snp_df, gtf_file, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def qtl_annotation(qtl_file, annotation_file):
        """Annotate QTL regions with genes."""
        result = anno.qtl_annotation(qtl_file, annotation_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def create_annotation_db(gtf_file, output_prefix):
        """Create local annotation database from GTF."""
        return anno.create_annotation_db(gtf_file, output_prefix)

    @mcp.tool()
    def predict_variant_effect(snp_chr=None, snp_pos=None, ref_allele=None, alt_allele=None,
                               gtf_file=None, simple_ref=None, simple_alt=None, region_type=None):
        """Predict variant functional effects for a single variant."""
        return anno.predict_variant_effect(
            snp_chr=snp_chr,
            snp_pos=snp_pos,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            gtf_file=gtf_file,
            simple_ref=simple_ref,
            simple_alt=simple_alt,
            region_type=region_type,
        )

    @mcp.tool()
    def get_genes_in_region(chr_name, start, end, annotation_file, mode='overlap'):
        """Get all genes in a genomic region."""
        result = anno.get_genes_in_region(chr_name, start, end, annotation_file, mode=mode)
        return result.to_dict() if result is not None else None

    # ========== QTL Tools (from qtl.py) ==========
    
    @mcp.tool()
    def detect_qtl_regions(gwas_file, p1=1e-7, p2=1e-5, p2n=5, window=500000):
        """Detect QTL regions from GWAS results."""
        result = qtl.detect_qtl(gwas_file, p1=p1, p2=p2, p2n=p2n, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_lead_snp(gwas_file, region_chr, region_start, region_end):
        """Get lead SNP in a region."""
        result = qtl.get_lead_snp(gwas_file, region_chr, region_start, region_end)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def identify_peak_snps(gwas_dir, p_threshold=1e-5, max_peaks=50):
        """Identify peak SNPs from GWAS results."""
        result = qtl.identify_peak_snps(gwas_dir, p_threshold=p_threshold, max_peaks=max_peaks)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def extract_qtl_genotypes(geno_prefix, qtl_regions, output_prefix):
        """Extract genotypes for QTL regions."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_regions) if isinstance(qtl_regions, dict) else qtl_regions
        return qtl.extract_qtl_genotypes(geno_prefix, qtl_data, output_prefix)

    @mcp.tool()
    def calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix):
        """Calculate haplotypes for SNPs in a QTL region."""
        result = qtl.calculate_qtl_haplo(geno_prefix, qtl_snp_list, output_prefix)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def plot_qtl_region(gwas_file, chr_val, start, end, output_file=None):
        """Plot GWAS results for QTL region."""
        return qtl.plot_qtl_region(gwas_file, chr_val, start, end, output_file=output_file)

    @mcp.tool()
    def qtl_summary(gwas_dir, output_prefix='qtl_summary'):
        """Generate summary of all QTL regions."""
        result = qtl.qtl_summary(gwas_dir, output_prefix=output_prefix)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def map_qtl_to_genes(qtl_df, annotation_file):
        """Map QTL regions to genes."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        result = qtl.map_qtl_to_genes(qtl_data, annotation_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def qtl_enrichment_test(qtl_genes, background_genes, gene_sets):
        """Test enrichment of QTL genes against gene sets."""
        result = qtl.qtl_enrichment_test(qtl_genes, background_genes, gene_sets)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_qtl_overlap(qtl1, qtl2):
        """Find overlapping QTL intervals between two QTL sets."""
        import pandas as pd
        qtl1_data = pd.DataFrame(qtl1) if isinstance(qtl1, dict) else qtl1
        qtl2_data = pd.DataFrame(qtl2) if isinstance(qtl2, dict) else qtl2
        result = qtl.get_qtl_overlap(qtl1_data, qtl2_data)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_qtl_bed(qtl_df, output_file):
        """Export QTL as BED file."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        return qtl.export_qtl_bed(qtl_data, output_file)

    # ========== Network Tools (from net.py) ==========

    @mcp.tool()
    def module_identify(edge_weight, module_size=3):
        """Identify network modules using ClusterONE or a NetworkX fallback."""
        import pandas as pd
        edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
        result = net.module_identify(edge_df, module_size=module_size)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def hub_identify(edge_weight, cluster_one_res):
        """Identify hub genes in each network module."""
        import pandas as pd
        edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
        cluster_df = pd.DataFrame(cluster_one_res) if isinstance(cluster_one_res, dict) else cluster_one_res
        cluster_res, hub_res = net.hub_identify(edge_df, cluster_df)
        return {
            "cluster_result": cluster_res.to_dict() if cluster_res is not None else None,
            "hub_result": hub_res.to_dict() if hub_res is not None else None,
        }


    # ========== Peak Test Tools (from peak.py) ==========
    
    @mcp.tool()
    def plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df):
        """Generate boxplots for QTL regions."""
        return peak.plot_qtl_boxplot(pheno_file, geno_prefix, qtl_df)

    @mcp.tool()
    def plot_grouped_boxplot(data_dict, output_file=None, test_method='t-test'):
        """Generate grouped boxplot."""
        return peak.plot_grouped_boxplot(data_dict, output_file=output_file, test_method=test_method)

    @mcp.tool()
    def haplotype_test(pheno_df, geno_df, snp_id, test_method='t-test'):
        """Perform statistical test for haplotype effect."""
        result = peak.haplotype_test(pheno_df, geno_df, snp_id, test_method=test_method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix):
        """Generate multi-trait Manhattan plot with QTL regions."""
        return peak.multi_trait_qtl_plot(gwas_dir, qtl_file, output_prefix)

    # ========== MultiProcess Tools (from multi.py) ==========
    
    @mcp.tool()
    def parallel_run_commands(cmds, num_threads=4):
        """Run commands in parallel."""
        return multi.parallel_run(multi.run_cmd, [(cmd,) for cmd in cmds], num_threads=num_threads)

    @mcp.tool()
    def get_optimal_threads(max_threads=None):
        """Get optimal number of threads."""
        return multi.get_optimal_threads(max_threads)

    # ========== Mendelian Randomization Tools (from mr.py) ==========

    @mcp.tool()
    def run_mr_analysis(exposure_gwas, outcome_gwas, method='ivw'):
        """Mendelian Randomization analysis (IVW, MR-Egger)."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        result = mr.mr_analysis(exp_df, out_df, method=method)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def calculate_mr_causal_estimate(exposure_gwas, outcome_gwas, snp_col='rs', beta_col='beta', se_col='se', method='ivw'):
        """Calculate MR causal estimate between two traits."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        result = mr.mr_causal_estimate(
            exp_df, out_df, snp_col=snp_col, beta_col=beta_col, se_col=se_col, method=method
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def test_mr_pleiotropy(exposure_gwas, outcome_gwas):
        """Test for horizontal pleiotropy (MR-Egger intercept)."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        return mr.test_pleiotropy(exp_df, out_df)

    @mcp.tool()
    def test_mr_heterogeneity(exposure_gwas, outcome_gwas):
        """Test for heterogeneity using IVW method."""
        import pandas as pd
        exp_df = pd.DataFrame(exposure_gwas) if isinstance(exposure_gwas, dict) else exposure_gwas
        out_df = pd.DataFrame(outcome_gwas) if isinstance(outcome_gwas, dict) else outcome_gwas
        return mr.heterogeneity_test(exp_df, out_df)

    @mcp.tool()
    def run_qtl_target_analysis(qtl_df, tf_genes, target_genes, window=500000):
        """QTL targeting analysis - find which TFs target which genes."""
        import pandas as pd
        qtl_data = pd.DataFrame(qtl_df) if isinstance(qtl_df, dict) else qtl_df
        result = mr.qtl_target_analysis(qtl_data, tf_genes, target_genes, window=window)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def format_qtl_for_mr(qtl_file, output_file=None):
        """Format QTL file for MR analysis."""
        result = mr.format_qtl_for_mr(qtl_file, output_file=output_file)
        return result.to_dict() if result is not None else None

    # ========== GO Enrichment Tools (from go.py) ==========

    @mcp.tool()
    def run_go_enrichment(gene_list, organism='Mouse', pvalue_cutoff=0.05, qvalue_cutoff=0.05,
                          mode='local', gene_sets_file=None, background_genes=None,
                          go_obo_file=None, obo_cache_dir=None, auto_download_obo=True,
                          drop_unmapped_terms=True):
        """GO enrichment analysis (BP, MF, CC)."""
        import pandas as pd
        genes = pd.DataFrame({'genes': gene_list}) if isinstance(gene_list, list) else gene_list
        result = go.go_enrich(
            genes,
            organism=organism,
            pvalue_cutoff=pvalue_cutoff,
            qvalue_cutoff=qvalue_cutoff,
            mode=mode,
            gene_sets_file=gene_sets_file,
            background_genes=background_genes,
            go_obo_file=go_obo_file,
            obo_cache_dir=obo_cache_dir,
            auto_download_obo=auto_download_obo,
            drop_unmapped_terms=drop_unmapped_terms,
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def run_gsea_analysis(gene_rank, gene_symbol_col='gene', score_col='score', organism='Mouse',
                          mode='local', gene_sets_file=None):
        """Gene Set Enrichment Analysis (GSEA)."""
        import pandas as pd
        df = pd.DataFrame(gene_rank) if isinstance(gene_rank, dict) else gene_rank
        result = go.gsea_enrich(
            None,
            gene_rank=df,
            gene_symbol_col=gene_symbol_col,
            score_col=score_col,
            organism=organism,
            mode=mode,
            gene_sets_file=gene_sets_file,
        )
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def plot_go_enrichment(enrich_result, output_prefix, plot_types=None, top_n=20, figsize=None,
                           split_ontology=True, label_col='Term',
                           bar_x='GeneRatio', bar_color='Adjusted P-value',
                           dot_x='GeneRatio', dot_color='Adjusted P-value',
                           dot_size='Count', sort_by=None, ascending=None):
        """Generate GO enrichment plots."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        if plot_types is None:
            plot_types = ['barplot', 'dotplot']
        if figsize is None:
            figsize = (8, 6)
        return go.go_plot(
            df,
            output_prefix,
            plot_types=plot_types,
            top_n=top_n,
            figsize=tuple(figsize),
            split_ontology=split_ontology,
            label_col=label_col,
            bar_x=bar_x,
            bar_color=bar_color,
            dot_x=dot_x,
            dot_color=dot_color,
            dot_size=dot_size,
            sort_by=sort_by,
            ascending=ascending,
        )

    @mcp.tool()
    def plot_gsea_results(gsea_result, output_prefix, top_n=20, figsize=None, format='png',
                          plot_mode='curve', terms=None, curve_terms=1, trace_terms=3, rank_metric=None):
        """Generate GSEA plots (curve, trace, NES summary)."""
        import pandas as pd
        df = pd.DataFrame(gsea_result) if isinstance(gsea_result, dict) else gsea_result
        if figsize is None:
            figsize = (8, 6)
        return go.gsea_plot(
            df,
            output_prefix,
            top_n=top_n,
            figsize=tuple(figsize),
            format=format,
            plot_mode=plot_mode,
            terms=terms,
            curve_terms=curve_terms,
            trace_terms=trace_terms,
            rank_metric=rank_metric,
        )

    @mcp.tool()
    def run_kegg_enrichment(gene_list, organism='mmu', pvalue_cutoff=0.05):
        """KEGG pathway enrichment analysis."""
        import pandas as pd
        genes = pd.DataFrame({'genes': gene_list}) if isinstance(gene_list, list) else gene_list
        result = go.kegg_enrich(genes, organism=organism, pvalue_cutoff=pvalue_cutoff)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def get_enrichr_libraries(organism='Mouse'):
        """Get available gene set libraries."""
        return go.enrichr_library_list(organism=organism)

    @mcp.tool()
    def extract_go_from_gtf(gtf_file, output_file=None):
        """Extract GO annotations from GTF file."""
        result = go.go_annotation_from_gtf(gtf_file, output_file=output_file)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def simplify_go_results(enrich_result, similarity_cutoff=0.7):
        """Simplify GO terms by removing redundancy."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        result = go.simplify_go_terms(df, similarity_cutoff=similarity_cutoff)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def export_go_report(enrich_result, output_file, format='excel'):
        """Export GO enrichment results."""
        import pandas as pd
        df = pd.DataFrame(enrich_result) if isinstance(enrich_result, dict) else enrich_result
        return go.export_go_report(df, output_file, format=format)

    mcp.run()


if __name__ == "__main__":
    if HAS_FASTMCP and HAS_MODULES:
        run_mcp_server()
    else:
        print("GenotypeMCP Server")
        print("=" * 40)
        print("Requirements: pip install fastmcp")
        print("\nAvailable modules:")
        print("  pheno: abundance_filter, missing_filter, scale_wrapper,")
        print("         pheno_imputer, outlier, blup, blue, trait_correct")
        print("  geno:  snp_qc, calculate_pca, calculate_ibd,")
        print("         vcf_to_plink, hapmap_to_plink, calculate_kinship,")
        print("         snp_pruning, snp_clumping")
        print("  gwas:  gwas_lm, gwas_lmm, gwas_plink, simple_gwas,")
        print("         add_rs_id_to_vcf, ensure_rs_id")
        print("         generate_clump, plink_clump, get_top_snps, calculate_lambda")
        print("  vis:   manhattan_plot, qq_plot, pca_plot, tsne_plot,")
        print("         ld_heatmap, phenotype_hist, phenotype_boxplot, gwas_summary")
        print("  anno:  parse_gtf, annotate_snps, qtl_annotation,")
        print("         create_annotation_db, predict_variant_effect, get_genes_in_region")
        print("  qtl:   detect_qtl_regions, get_lead_snp, identify_peak_snps,")
        print("         extract_qtl_genotypes, calculate_qtl_haplo,")
        print("         plot_qtl_region, qtl_summary, map_qtl_to_genes,")
        print("         qtl_enrichment_test, get_qtl_overlap, export_qtl_bed")
        print("  mr:    run_mr_analysis, calculate_mr_causal_estimate,")
        print("         test_mr_pleiotropy, test_mr_heterogeneity,")
        print("         run_qtl_target_analysis, format_qtl_for_mr")
        print("  go:    run_go_enrichment, run_gsea_analysis, plot_go_enrichment,")
        print("         plot_gsea_results, run_kegg_enrichment, get_enrichr_libraries,")
        print("         extract_go_from_gtf, simplify_go_results, export_go_report")
        print("  net:   module_identify, hub_identify")
        print("  multi: parallel_run_commands, get_optimal_threads")
        print("  peak:  plot_qtl_boxplot, plot_grouped_boxplot, haplotype_test, multi_trait_qtl_plot")
