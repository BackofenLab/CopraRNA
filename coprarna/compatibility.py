"""Thread-local choice between corrected behavior and the frozen legacy oracle."""

from contextlib import contextmanager
from contextvars import ContextVar

_LEGACY = ContextVar("coprarna_legacy", default=False)


def legacy_enabled():
    return _LEGACY.get()


@contextmanager
def legacy_behavior(enabled=True):
    token = _LEGACY.set(bool(enabled))
    try:
        yield
    finally:
        _LEGACY.reset(token)
