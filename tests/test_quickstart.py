import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "quickstart.py"


def test_quickstart_script_runs_end_to_end():
    out = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                         timeout=120, check=True).stdout
    assert "[FLAG] panel_consistency" in out
    assert "[PASS] control_check" in out
    assert "[FLAG] size_confounding" in out and "[PASS] size_confounding" in out
