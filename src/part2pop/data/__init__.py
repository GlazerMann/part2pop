import io
from functools import lru_cache
from typing import TextIO
import importlib_resources

@lru_cache(maxsize=32)
def _read_dataset(filename: str, encoding: str) -> str:
    """Read and cache a data file from package resources."""
    resource = importlib_resources.files("part2pop.data").joinpath(filename)
    with resource.open("r", encoding=encoding) as f:
        return f.read()

def open_dataset(filename: str, encoding: str = "utf-8") -> TextIO:
    """Return a fresh text stream over cached package-data contents.

    The public contract is a readable text stream with an independent cursor
    and lifetime on every call. Callers must not rely on OS-backed file
    attributes such as ``name`` or ``fileno()``.
    """
    return io.StringIO(_read_dataset(filename, encoding))
