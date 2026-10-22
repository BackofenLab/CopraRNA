"""Thread-local choice between corrected behavior and the frozen legacy oracle."""

from contextlib import contextmanager
from contextvars import ContextVar

_LEGACY = ContextVar("coprarna_legacy", default=False)


def legacy_enabled():
    """Report whether retained historical defects are enabled in this context."""
    return _LEGACY.get()


@contextmanager
def legacy_behavior(enabled=True):
    """Scope historical helper behavior and restore the prior setting on exit.

    ContextVar isolation keeps concurrent batch jobs independent and restoration
    works on exceptions. Production prediction explicitly selects corrected behavior.
    """
    token = _LEGACY.set(bool(enabled))
    try:
        yield
    finally:
        _LEGACY.reset(token)
