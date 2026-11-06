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


# A dry run resolves input *files*, so every assertion about its exit status needs the
# deposits, which the repository does not carry. The assertions that do not need them --
# that run.sh reaches snakemake at all, and that a target names a rule -- are kept
# unconditional, so a checkout without deposits still checks them.
HAS_DEPOSITS = (ROOT / "data" / "results" / "e1" / "sweep.jsonl").exists()
needs_deposits = pytest.mark.skipif(
    not HAS_DEPOSITS, reason="no deposits here: a dry run cannot resolve its input files"
)


@needs_snakemake
def test_run_sh_with_no_arguments_reaches_snakemake():
    # The empty-array expansion under `set -u` fails on bash before 4.4; the default
    # invocation is the one a reader tries first. Reaching the DAG is what that tests, and
    # it holds with or without deposits; the exit status does not.
    got = dry_run()
    assert "Building DAG of jobs" in got.stdout, got.stderr[-2000:]
    if HAS_DEPOSITS:
        assert got.returncode == 0, got.stderr[-2000:]


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
    # Without deposits snakemake still reports an unknown path as a missing *rule*, which is
    # the half of this that does not depend on the data being present.
    assert "MissingRuleException" not in got.stderr, got.stderr[-2000:]
    if HAS_DEPOSITS:
        assert got.returncode == 0, got.stderr[-2000:]


def test_the_readme_has_path_examples():
    assert readme_path_targets()


def test_the_readme_lists_targets():
    assert len(readme_targets()) >= 10


@pytest.mark.parametrize("name", ["e1", "e3"])
def test_driver_configs_render_without_host_paths(tmp_path, name):
    template = ROOT / "drivers" / name / f"{name}.config.json"
    out = tmp_path / "config.json"
    values = {
        k: str(tmp_path / k)
        for k in ["bin", "derna2_bin", "inputs", "out", "root", "lineardesign_dir", "lcds_dir"]
    }
    values["pin"] = "0123abcd"
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
    # The pin is a commit, not a path, and reaches the rows as written.
    if name in ("e1", "e3"):
        assert json.loads(text)["derna2_pin"] == "0123abcd"


@pytest.mark.parametrize("name", ["e1", "e3"])
def test_every_input_a_driver_config_names_is_in_the_repo(name):
    template = (ROOT / "drivers" / name / f"{name}.config.json").read_text()
    named = set(re.findall(r'"\{inputs\}/([^"]+)"', template))
    assert named
    assert sorted(p for p in named if not (ROOT / "inputs" / p).exists()) == []


def test_render_refuses_a_missing_value(tmp_path):
    template = ROOT / "drivers" / "e1" / "e1.config.json"
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
