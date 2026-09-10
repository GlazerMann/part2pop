from pathlib import Path
import subprocess
import sys

def test_analysis_factory_compiles_without_syntax_warnings():
    repository_root = Path(__file__).resolve().parents[2]
    factory_dir = (
        repository_root
        / "src"
        / "part2pop"
        / "analysis"
        / "population"
        / "factory"
    )

    result = subprocess.run(
        [
            sys.executable,
            "-W"
            "error::SyntaxWarning",
            "-m",
            "compileall",
            "-q",
            "-f",
            str(factory_dir),
        ],
        cwd=repository_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
