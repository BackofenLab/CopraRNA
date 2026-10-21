import json
import subprocess
import sys

import pytest

from coprarna.cli import main


def test_single_helper_uses_workdir_and_clean_stdout(tmp_path, capsys):
    (tmp_path / "input.fa").write_text(">NC_000001 description\nAC\nGT\n")
    assert main(["helper", "--workdir", str(tmp_path), "format-fasta", "input.fa"]) == 0
    assert capsys.readouterr().out == ">NC_000001\nACGT\n"


def _manifest(tmp_path, second="input.fa"):
    for name in ("job1", "job2"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "input.fa").write_text(f">{name}\nACGT\n")
    items = [
        {
            "id": name,
            "helper": "format-fasta",
            "args": [second if name == "job2" else "input.fa"],
            "workdir": name,
            "stdout": "formatted.fa",
        }
        for name in ("job1", "job2")
    ]
    manifest = tmp_path / "jobs.json"
    manifest.write_text(json.dumps(items))
    return manifest, items


def test_parallel_batch_preserves_job_order_and_files(tmp_path, capsys):
    manifest, _ = _manifest(tmp_path)
    assert main(["batch", str(manifest), "--jobs", "2"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert [r["id"] for r in result] == ["job1", "job2"]
    assert all(r["status"] == "ok" for r in result)
    assert (tmp_path / "job1/formatted.fa").read_text() == ">job1\nACGT\n"
    assert (tmp_path / "job2/formatted.fa").read_text() == ">job2\nACGT\n"


def test_batch_reports_failures_without_losing_successful_jobs(tmp_path, capsys):
    manifest, _ = _manifest(tmp_path, "missing.fa")
    assert main(["batch", str(manifest), "--jobs", "2"]) == 1
    result = json.loads(capsys.readouterr().out)
    assert [r["status"] for r in result] == ["ok", "error"]
    assert (tmp_path / "job1/formatted.fa").exists()
    assert not (tmp_path / "job2/formatted.fa").exists()


@pytest.mark.parametrize(
    "problem",
    [
        "same-directory",
        "same-id",
        "outside-stdout",
        "existing-stdout",
        "bad-args",
        "unknown-helper",
    ],
)
def test_batch_preflight_prevents_partial_execution(tmp_path, capsys, problem):
    manifest, items = _manifest(tmp_path)
    if problem == "same-directory":
        items[1]["workdir"] = "job1"
    elif problem == "same-id":
        items[1]["id"] = "job1"
    elif problem == "outside-stdout":
        items[1]["stdout"] = "../outside.fa"
    elif problem == "existing-stdout":
        (tmp_path / "job2/formatted.fa").write_text("keep this")
    elif problem == "bad-args":
        items[1]["args"] = []
    else:
        items[1]["helper"] = "nonexistent-helper"
    manifest.write_text(json.dumps(items))
    assert main(["batch", str(manifest), "--jobs", "2"]) == 2
    assert not (tmp_path / "job1/formatted.fa").exists()
    if problem == "existing-stdout":
        assert (tmp_path / "job2/formatted.fa").read_text() == "keep this"


def test_module_cli_help_and_failure_exit_codes(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "coprarna", "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0
    assert "comparative RNA target prediction" in result.stdout
    assert "Traceback" not in result.stderr
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "coprarna",
            "helper",
            "--workdir",
            str(tmp_path),
            "format-fasta",
            "missing.fa",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert result.stdout == ""
    assert "Traceback" not in result.stderr


def test_prediction_batch_options_are_passed_per_job(tmp_path, capsys, monkeypatch):
    from coprarna import pipeline

    manifest, items = _manifest(tmp_path)
    for i, item in enumerate(items):
        item.update(helper="predict", options={"cores": i + 1, "cons": 2})
    manifest.write_text(json.dumps(items))
    seen = []

    def fake_prediction(srna, directory, config):
        seen.append((directory.name, config.cores, config.cons))
        return ""

    monkeypatch.setattr(pipeline, "run_prediction", fake_prediction)
    assert main(["batch", str(manifest), "--jobs", "2"]) == 0
    assert sorted(seen) == [("job1", 1, 2), ("job2", 2, 2)]
