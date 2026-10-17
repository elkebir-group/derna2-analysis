"""The pipeline's wiring: every documented target resolves, run.sh works with no extra
arguments, and the driver configs render without host paths."""

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake") or str(Path(sys.executable).parent / "snakemake")
needs_snakemake = pytest.mark.skipif(
    not Path(SNAKEMAKE).exists() and not shutil.which("snakemake"),
    reason="snakemake not installed (pipeline extra)",
)


def dry_run(*targets):
    """run.sh with no options of its own: `--dry-run` goes after `--`, straight to snakemake,
    so run.sh's own argument array stays empty, the case that breaks bash before 4.4."""
    env = {"SNAKEMAKE": SNAKEMAKE, "PATH": "/usr/bin:/bin"}
    cmd = ["bash", str(ROOT / "pipelines" / "run.sh"), "--", "--dry-run", *targets]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, env=env)


def readme_targets():
    """Rule names the README tells a reader to run: every `name` in its target tables."""
    text = (ROOT / "README.md").read_text()
    return sorted(set(re.findall(r"^\| `([a-z0-9_]+)` \|", text, flags=re.M)))


@needs_snakemake
def test_run_sh_with_no_arguments_reaches_snakemake():
    # The empty-array expansion under `set -u` fails on bash before 4.4; the default
    # invocation is the one a reader tries first.
    got = dry_run()
    assert got.returncode == 0, got.stderr[-2000:]
    assert "Building DAG of jobs" in got.stdout


@needs_snakemake
@pytest.mark.parametrize("target", readme_targets())
def test_every_readme_target_is_a_rule(target):
    got = dry_run(target)
    assert "MissingRuleException" not in got.stderr and "Unknown target" not in got.stderr, (
        got.stderr[-2000:]
    )


def readme_path_targets():
    """Output paths the README's examples pass to run.sh."""
    text = (ROOT / "README.md").read_text()
    return sorted(set(re.findall(r"run\.sh -- (data/\S+)", text)))


@needs_snakemake
@pytest.mark.parametrize("target", readme_path_targets())
def test_every_readme_path_target_has_a_rule(target):
    got = dry_run(target)
    assert got.returncode == 0, got.stderr[-2000:]


def test_the_readme_has_path_examples():
    assert readme_path_targets()


def test_the_readme_lists_targets():
    assert len(readme_targets()) >= 10


@pytest.mark.parametrize("name", ["e1", "e2"])
def test_driver_configs_render_without_host_paths(tmp_path, name):
    template = ROOT / "drivers" / "benchmark" / f"{name}.config.json"
    out = tmp_path / "config.json"
    values = {
        k: str(tmp_path / k) for k in ["bin", "bin_e1", "inputs", "out", "root", "derna1_params"]
    }
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "pipelines" / "scripts" / "render_driver_config.py"),
            str(template),
            str(out),
            *[f"{k}={v}" for k, v in values.items() if "{" + k + "}" in template.read_text()],
        ],
        check=True,
    )
    text = out.read_text()
    json.loads(text)
    assert "/scratch/" not in text and "/home/" not in text
    # The beam width placeholder is the driver's own and must survive rendering.
    if name == "e1":
        assert "b{w}" in text or "{w}" in text


def test_render_refuses_a_missing_value(tmp_path):
    template = ROOT / "drivers" / "benchmark" / "e1.config.json"
    got = subprocess.run(
        [
            sys.executable,
            str(ROOT / "pipelines" / "scripts" / "render_driver_config.py"),
            str(template),
            str(tmp_path / "c.json"),
            f"bin={tmp_path}",
        ],
        capture_output=True,
        text=True,
    )
    assert got.returncode != 0 and "no value" in got.stderr
