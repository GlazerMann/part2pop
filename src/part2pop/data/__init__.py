import io
from functools import lru_cache
import importlib_resources

@lru_cache(maxsize=64)
def _read_dataset(filename: str, encoding: str) -> str:
    """Read and cache a data file from package resources."""
    resource = importlib_resources.files("part2pop.data").joinpath(filename)
    with resource.open("r", encoding=encoding) as f:
        return f.read()

def open_dataset(filename: str, encoding: str = "utf-8") -> io.StringIO:
    """Return a readable stream for a cached package data file."""
    return io.StringIO(_read_dataset(filename, encoding))
