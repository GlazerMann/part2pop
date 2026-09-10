import io
from functools import lru_cache
import importlib_resources

@lru_cache(maxsize=32)
def _read_dataset(filename: str, encoding: str) -> str:
    """Read and cache a data file from package resources."""
    resource = importlib_resources.files("part2pop.data").joinpath(filename)
    with resource.open("r", encoding=encoding) as f:
        return f.read()

def open_dataset(filename: str, encoding: str = "utf-8") -> io.StringIO:
    """Return a fresh readable stream over cached package-data text.
    Every call has an independent cursor and lifetime. The returned object is
    intentionally an in-memory text stream; OS-backed file attributes such as
    ``name`` and ``fileno()`` are not part of this helper's contract.
    """
    return io.StringIO(_read_dataset(filename, encoding))
