#!/bin/bash
# Quick setup script for MRBIGR2

echo "Setting up MRBIGR2..."

# Create directories
mkdir -p data output

# Install Python dependencies
pip install pandas numpy scipy scikit-learn gseapy networkx matplotlib seaborn statsmodels -q

# Test all modules
cd "$(dirname "$0")"
python3 -c "
import sys
sys.path.insert(0, 'src')

# Test all 11 modules
print('Testing modules...')
from pheno import blup, scale_wrapper
from geno import get_snp_stats, calculate_pca
from gwas import gwas_lm, qq_plot
from vis import pca_plot, phenotype_boxplot
from anno import parse_gtf, annotate_snps_simple
from qtl import detect_qtl, get_lead_snp
from peak import plot_qtl_boxplot, haplotype_test
from mr import mr_analysis, heterogeneity_test
from go import go_enrich, go_annotation_from_gtf
from net import check_java, check_cluster_one
from multi import batch_process, chunk_process

print('All 11 modules imported successfully!')
"

echo ""
echo "MRBIGR2 setup complete!"
echo ""
echo "Usage:"
echo "  # Import in Python"
echo "  import sys; sys.path.insert(0, 'src')"
echo "  from geno import *"
echo "  from gwas import *"
echo "  from pheno import *"
echo "  # ... other modules"
echo ""
echo "  # Run MCP server"
echo "  python src/server.py"
echo ""
echo "Test data: ./test_run/real_data/"
