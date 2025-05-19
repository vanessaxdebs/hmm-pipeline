import pytest
from pathlib import Path
from pipeline.modules.build import build_hmm

def test_hmm_creation(tmp_path):
    test_aln = tmp_path / "test.sto"
    test_aln.write_text("# STOCKHOLM 1.0\nseq1 ACDEF\n//")
    build_hmm(test_aln, tmp_path / "test.hmm")
    assert (tmp_path / "test.hmm").exists()
