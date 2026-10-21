#!/usr/bin/env python3
"""Explicitly refresh the golden files after reviewing a legacy fixture change."""
import json
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.support.legacy import FIXTURES, capture, interpreter


def main():
    for specfile in sorted(FIXTURES.glob("*/case.json")):
        case = specfile.parent
        spec = json.loads(specfile.read_text())
        with tempfile.TemporaryDirectory(prefix="coprarna-oracle-") as directory:
            outputs = capture(case, Path(directory))
        expected = case / "expected"
        expected.mkdir(exist_ok=True)
        for name, data in outputs.items():
            (expected / name).parent.mkdir(parents=True, exist_ok=True)
            (expected / name).write_bytes(data)
        print(f"Captured {case.name} using {interpreter(spec['language'])}")


if __name__ == "__main__":
    main()
