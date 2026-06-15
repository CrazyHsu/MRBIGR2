"""Regression: plot_qtl_boxplot / peak_boxplot_runner must produce a
*diagnosable* result and never return None.

The longjob worker codes a ``None`` return as ``returncode=1`` (failed) with
empty stdout/stderr. Previously, when no boxplots could be produced (e.g. the
QTL table has no ``phe_name``/``trait`` column to map each SNP to a multi-trait
phenotype), the runner returned ``None`` -> an undiagnosable rc=1 failure.
It must instead complete with an actionable ``reason``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
REPO_SRC = REPO_ROOT / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from mrbigr.mcp import runners


def test_peak_boxplot_runner_returns_reason_not_none_when_unmappable(tmp_path: Path) -> None:
    # Multi-trait phenotype + QTL table without a phe_name/trait column => the
    # SNP->trait mapping is ambiguous. This returns at the upfront check, before
    # any PLINK call, so no genotype is needed.
    pheno = tmp_path / "pheno.csv"
    pheno.write_text("ID,traitA,traitB\ns1,1.0,2.0\ns2,3.0,4.0\n", encoding="utf-8")
    qtl = tmp_path / "qtl.csv"
    qtl.write_text("CHR,qtl_start,qtl_end,SNP,P\n1,100,200,snp1,1e-8\n", encoding="utf-8")
    out_dir = tmp_path / "box"
    result_json = out_dir / "result.json"

    result = runners.peak_boxplot_runner(
        pheno_file=str(pheno),
        geno_prefix=str(tmp_path / "nogeno"),  # never read for this path
        qtl_csv=str(qtl),
        output_dir=str(out_dir),
        test_method="t-test",
        result_json=str(result_json),
    )

    # Never None: the worker would otherwise mark this rc=1/failed with empty logs.
    assert result is not None
    assert result_json.is_file()
    payload = json.loads(result_json.read_text(encoding="utf-8"))
    assert payload["n_plots"] == 0
    assert payload["output_files"] == []
    assert "phe_name" in (payload["reason"] or "")


def test_peak_boxplot_core_writes_diagnostics_json(tmp_path: Path) -> None:
    from mrbigr.core import peak

    pheno = tmp_path / "pheno.csv"
    pheno.write_text("ID,traitA,traitB\ns1,1.0,2.0\ns2,3.0,4.0\n", encoding="utf-8")
    qtl = tmp_path / "qtl.csv"
    qtl.write_text("CHR,SNP\n1,snp1\n", encoding="utf-8")
    out_dir = tmp_path / "box"

    outputs = peak.plot_qtl_boxplot(str(pheno), str(tmp_path / "nogeno"), str(qtl), output_dir=str(out_dir))

    assert outputs == []
    diag_file = out_dir / "qtl_boxplot_diagnostics.json"
    assert diag_file.is_file()
    diag = json.loads(diag_file.read_text(encoding="utf-8"))
    assert diag["n_plotted"] == 0
    assert diag["reason"]
