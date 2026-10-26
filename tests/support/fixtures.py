"""Read independently captured scientific regression fixtures."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests" / "fixtures"


def canonicalize(name: str, value: bytes, spec: dict) -> bytes:
    if name in spec.get("unordered", []):
        header, *rows = value.splitlines(keepends=True)
        return header + b"".join(sorted(rows))
    return value
