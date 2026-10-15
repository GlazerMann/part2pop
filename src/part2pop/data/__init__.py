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
    """Return a new text stream over cached packaged resource contents.

    Cache entries are keyed by (filename, encoding), with up to 32 entries;
    resources are assumed immutable for the process lifetime. The stream has
    its own cursor and is not a real OS file handle. Runtime disk changes may
    not be seen until the underlying cache is cleared or the process restarts.
    """
    return io.StringIO(_read_dataset(filename, encoding))
