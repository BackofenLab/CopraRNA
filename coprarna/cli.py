"""Command-line access to validated migration helpers and isolated batch jobs."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from importlib.resources import files
import json
import logging
from pathlib import Path
import sys

from . import __version__
from .helpers import HELPERS, resolve, run_helper

LOG = logging.getLogger("coprarna")


def _parser():
    parser = argparse.ArgumentParser(
        description="CopraRNA comparative RNA target prediction with Python 3."
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    helper = commands.add_parser(
        "helper", help="Execute a validated port in one working directory"
    )
    helper.add_argument(
        "--legacy", action="store_true", help="Reproduce the frozen legacy defects"
    )
    helper.add_argument("--workdir", type=Path, default=Path.cwd())
    helper.add_argument("name", help="Helper name or original script basename")
    helper.add_argument("args", nargs=argparse.REMAINDER)
    batch = commands.add_parser(
        "batch", help="Run independent prediction or helper jobs from a JSON manifest"
    )
    batch.add_argument("manifest", type=Path)
    batch.add_argument(
        "--jobs", type=int, default=1, help="Maximum concurrent jobs (default: 1)"
    )
    run = commands.add_parser(
        "run", help="Run a complete comparative target prediction"
    )
    run.add_argument("--workdir", type=Path, default=Path.cwd())
    run.add_argument("--srnaseq", "-srnaseq", type=Path, default=Path("input_sRNA.fa"))
    from .pipeline import Config
    from dataclasses import fields

    for field in fields(Config):
        if field.name == "seed":
            run.add_argument("--seed", type=int, default=None)
        elif isinstance(field.default, bool):
            run.add_argument("--" + field.name, "-" + field.name, action="store_true")
        else:
            run.add_argument(
                "--" + field.name,
                "-" + field.name,
                type=type(field.default),
                default=field.default,
            )
    commands.add_parser("status", help="Print the complete migration inventory as JSON")
    commands.add_parser("helpers", help="List helper names and positional arguments")
    return parser


def _load_jobs(manifest: Path) -> list[dict]:
    raw = json.loads(manifest.read_text())
    if not isinstance(raw, list) or not raw:
        raise ValueError("Batch manifest must be a nonempty JSON list")
    jobs, identifiers, directories = [], set(), []
    for item in raw:
        if not isinstance(item, dict) or set(item) - {
            "id",
            "helper",
            "args",
            "workdir",
            "stdout",
            "options",
        }:
            raise ValueError(
                "Each job must contain only id, helper, args, workdir and optional stdout"
            )
        if not all(key in item for key in ("id", "helper", "args", "workdir")):
            raise ValueError("Each job requires id, helper, args and workdir")
        identifier = item["id"]
        if (
            not isinstance(identifier, str)
            or not identifier
            or identifier in identifiers
        ):
            raise ValueError("Job IDs must be unique nonempty strings")
        if (
            not isinstance(item["helper"], str)
            or not isinstance(item["args"], list)
            or not all(isinstance(a, str) for a in item["args"])
        ):
            raise ValueError(f"Invalid helper or argument list for job {identifier}")
        settings = item.get("options", {})
        if not isinstance(settings, dict) or (settings and item["helper"] != "predict"):
            raise ValueError("options is a configuration object for predict jobs only")
        if item["helper"] == "predict":
            from .pipeline import Config

            try:
                Config(**settings)
            except TypeError as exc:
                raise ValueError(str(exc)) from exc
        name = resolve(item["helper"])
        spec = HELPERS[name]
        if len(item["args"]) < len(spec.arguments) or (
            not spec.variadic and len(item["args"]) != len(spec.arguments)
        ):
            raise ValueError(
                f"Wrong argument count for job {identifier}: {name} expects {spec.arguments}"
            )
        if not isinstance(item["workdir"], str):
            raise ValueError(f"Invalid working directory for job {identifier}")
        directory = (manifest.resolve().parent / item["workdir"]).resolve()
        if not directory.is_dir():
            raise ValueError(
                f"Working directory does not exist for job {identifier}: {directory}"
            )
        if any(
            directory == other
            or directory in other.parents
            or other in directory.parents
            for other in directories
        ):
            raise ValueError(
                "Batch jobs require separate, non-nested working directories"
            )
        output = None
        if "stdout" in item:
            if not isinstance(item["stdout"], str) or not item["stdout"]:
                raise ValueError(f"Invalid stdout path for job {identifier}")
            relative = Path(item["stdout"])
            output = (directory / relative).resolve()
            if relative.is_absolute() or directory not in output.parents:
                raise ValueError(
                    "A job stdout path must stay inside its working directory"
                )
            if output.exists():
                raise ValueError(f"Refusing to overwrite stdout file: {output}")
            if not output.parent.is_dir():
                raise ValueError(
                    f"Stdout parent directory does not exist: {output.parent}"
                )
        identifiers.add(identifier)
        directories.append(directory)
        jobs.append(
            dict(
                id=identifier,
                helper=name,
                args=item["args"],
                workdir=directory,
                output=output,
                options=settings,
            )
        )
    return jobs


def _run_job(job: dict) -> dict:
    record = {"id": job["id"], "helper": job["helper"]}
    stream = None
    try:
        if job["output"] is not None:
            stream = job["output"].open("x", encoding="utf-8", newline="")
        if job["helper"] == "predict":
            from .pipeline import Config, run_prediction

            result = run_prediction(
                Path(job["args"][0]), job["workdir"], Config(**job["options"])
            )
        else:
            result = run_helper(job["helper"], job["args"], job["workdir"])
        if stream is None:
            record["stdout"] = result
        else:
            stream.write(result)
            record["stdout_file"] = str(job["output"])
        record["status"] = "ok"
    except Exception as exc:
        record.update(status="error", error=str(exc))
        LOG.error("Job %s: %s", job["id"], exc)
        if stream is not None:
            stream.close()
            job["output"].unlink()
    finally:
        if stream is not None:
            stream.close()
    return record


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr
    )
    try:
        if args.command == "run":
            from .pipeline import Config, run_prediction
            from dataclasses import fields

            config = Config(
                **{field.name: getattr(args, field.name) for field in fields(Config)}
            )
            run_prediction(args.srnaseq, args.workdir, config)
        elif args.command == "helpers":
            for name, spec in HELPERS.items():
                print(
                    name
                    + (" " + " ".join(spec.arguments) if spec.arguments else "")
                    + (" ..." if spec.variadic else "")
                )
        elif args.command == "status":
            print(files("coprarna").joinpath("migration.json").read_text(), end="")
        elif args.command == "helper":
            sys.stdout.write(
                run_helper(args.name, args.args, args.workdir, legacy=args.legacy)
            )
        elif args.command == "batch":
            if args.jobs < 1:
                raise ValueError("--jobs must be at least 1")
            jobs = _load_jobs(args.manifest)
            with ThreadPoolExecutor(max_workers=args.jobs) as pool:
                records = list(pool.map(_run_job, jobs))
            json.dump(records, sys.stdout, indent=2)
            print()
            return int(any(record["status"] == "error" for record in records))
    except (
        OSError,
        ValueError,
        IndexError,
        KeyError,
        ArithmeticError,
        RuntimeError,
    ) as exc:
        LOG.error("%s", exc)
        return 2
    except KeyboardInterrupt:
        LOG.error("Interrupted")
        return 130
    return 0
