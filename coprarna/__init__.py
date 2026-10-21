"""Python 3 ports of CopraRNA, with explicit legacy compatibility."""

__version__ = "3.0.0.dev0"

# R's finite-difference fit is sensitive to BLAS reduction order. Match the
# validated Haswell/Zen AVX2 kernels before NumPy/SciPy initialize OpenBLAS.
# Preserve an explicit caller setting and avoid unsupported CPU instructions.
import os
import platform
from pathlib import Path

if platform.system() == "Linux" and platform.machine() == "x86_64":
    try:
        if "avx2" in Path("/proc/cpuinfo").read_text():
            os.environ.setdefault("OPENBLAS_CORETYPE", "Haswell")
    except OSError:
        pass
