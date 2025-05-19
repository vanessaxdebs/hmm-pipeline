import subprocess
from pathlib import Path

def build_hmm(alignment: Path, hmm_output: Path):
    """Build HMM from Stockholm alignment."""
    subprocess.run(
        f"hmmbuild {hmm_output} {alignment}",
        shell=True,
        check=True
    )
