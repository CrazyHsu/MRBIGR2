"""Input coercion tests for QTL haplotype SNP lists."""
from __future__ import annotations

from pathlib import Path

import pytest

from mrbigr.core import qtl


def test_coerce_snp_list_accepts_list() -> None:
    assert qtl._coerce_snp_list(["chr1.s_1", "chr2.s_2"]) == ["chr1.s_1", "chr2.s_2"]


def test_coerce_snp_list_accepts_json_string() -> None:
    assert qtl._coerce_snp_list('["chr1.s_1", "chr2.s_2"]') == ["chr1.s_1", "chr2.s_2"]


def test_coerce_snp_list_accepts_delimited_string() -> None:
    assert qtl._coerce_snp_list("chr1.s_1, chr2.s_2\nchr3.s_3") == [
        "chr1.s_1",
        "chr2.s_2",
        "chr3.s_3",
    ]


def test_coerce_snp_list_accepts_single_snp_string() -> None:
    assert qtl._coerce_snp_list("chr1.s_1") == ["chr1.s_1"]


def test_coerce_snp_list_accepts_file_path(tmp_path: Path) -> None:
    snp_file = tmp_path / "snps.txt"
    snp_file.write_text("chr1.s_1\nchr2.s_2\n", encoding="utf-8")

    assert qtl._coerce_snp_list(str(snp_file)) == ["chr1.s_1", "chr2.s_2"]


def test_coerce_snp_list_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="qtl_snp_list"):
        qtl._coerce_snp_list("")
