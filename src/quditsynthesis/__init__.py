"""QuditSynthesis — quantum circuit synthesis over cyclotomic unitary groups.

Two backends share one API: "rust" (default, compiled PyO3 extension) and "python"
(pure-Python reference, fallback when the extension isn't built). Select with the
QUDITSYNTHESIS_BACKEND environment variable or set_backend() — before importing
names from this package, since `from quditsynthesis import operator` binds the
then-active backend.
"""

import os as _os
import warnings as _warnings

__all__ = [
    "cyclotomic_ring",
    "cyclotomic_element",
    "operator",
    "state",
    "gauss_sequence",
    "circulant",
    "multiply_many",
    "multiply_selected",
    "set_backend",
    "get_backend",
]

# Public names delegated to the active backend, and their rust-sidecar attribute names.
_BACKENDS = ("rust", "python")
_DELEGATED = __all__[: __all__.index("set_backend")]
_RUST_NAMES = {
    "cyclotomic_ring": "cyclotomic_ring_rust",
    "cyclotomic_element": "cyclotomic_element_rust",
    "operator": "operator_rust",
    "state": "state_rust",
    "gauss_sequence": "gauss_sequence_rust",
    "circulant": "circulant_rust",
    "multiply_many": "multiply_many_rust",
    "multiply_selected": "multiply_selected_rust",
}

_active = None  # (name, module) of the selected backend


def set_backend(name):
    """Select the compute backend: ``"rust"`` (default) or ``"python"``.

    Falls back to ``"python"`` with a warning when the rust extension is not
    built. Existing objects keep their original backend — backends do not mix.
    """
    global _active
    name = str(name).strip().lower()
    if name not in _BACKENDS:
        raise ValueError(f"unknown backend {name!r} (expected 'rust' or 'python')")
    if name == "rust":
        try:
            from quditsynthesis import datastructures_rust as module
        except ImportError:
            _warnings.warn(
                "quditsynthesis._rust extension not built; falling back to the "
                "python backend. Build with: maturin develop --release",
                RuntimeWarning,
                stacklevel=2,
            )
            name = "python"
    if name == "python":
        from quditsynthesis import datastructures as module
    _active = (name, module)


def get_backend():
    """Name of the active backend: ``"rust"`` or ``"python"``."""
    return _active[0]


set_backend(_os.environ.get("QUDITSYNTHESIS_BACKEND", "rust"))


def __getattr__(name):
    # PEP 562: resolve names against the active backend at access time, so set_backend() affects names not yet imported.
    if name in _DELEGATED:
        backend_name, module = _active
        if backend_name == "rust":
            name = _RUST_NAMES[name]
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
