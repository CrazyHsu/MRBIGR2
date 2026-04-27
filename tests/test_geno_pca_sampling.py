"""Regression tests for PCA SNP sampling."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mrbigr.core import geno


def test_calculate_pca_uses_even_snp_sampling() -> None:
    captured = {}
    original_read = geno.read_plink_bed
    original_pca = geno.PCA

    class FakePCA:
        def __init__(self, n_components):
            self.n_components = n_components
            self.explained_variance_ratio_ = np.array([0.7, 0.3], dtype=float)

        def fit_transform(self, X):
            captured["X_shape"] = X.shape
            return np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)

    def fake_read_plink_bed(input_prefix, max_snps=None, snp_selection="head"):
        captured["input_prefix"] = input_prefix
        captured["max_snps"] = max_snps
        captured["snp_selection"] = snp_selection
        G = np.array([[0, 1], [1, 0]], dtype=np.int8)
        snps = pd.DataFrame({"snp": ["rs1", "rs2"]})
        fam = pd.DataFrame({"id": ["s1", "s2"]})
        return G, snps, fam

    try:
        geno.read_plink_bed = fake_read_plink_bed
        geno.PCA = FakePCA
        pc_df, var_ratio = geno.calculate_pca("geno_prefix", n_components=2, max_snps=20000)
    finally:
        geno.read_plink_bed = original_read
        geno.PCA = original_pca

    assert captured["input_prefix"] == "geno_prefix"
    assert captured["max_snps"] == 20000
    assert captured["snp_selection"] == "even"
    assert captured["X_shape"] == (2, 2)
    assert list(pc_df.columns) == ["sample", "PC1", "PC2"]
    assert list(pc_df["sample"]) == ["s1", "s2"]
    assert np.allclose(var_ratio, np.array([0.7, 0.3]))
