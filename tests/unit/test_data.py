import io

from part2pop import data


class _FakeResource:
    def __init__(self):
        self.opens = 0

    def joinpath(self, filename):
        assert filename == "sample.txt"
        return self

    def open(self, mode, encoding):
        assert mode == "r"
        assert encoding == "utf-8"
        self.opens += 1
        return io.StringIO("alpha\nbeta\n")


def test_open_dataset_caches_contents_but_returns_independent_streams(monkeypatch):
    resource = _FakeResource()
    monkeypatch.setattr(data.importlib_resources, "files", lambda _: resource)
    data._read_dataset.cache_clear()

    with data.open_dataset("sample.txt") as first:
        assert first.readline() == "alpha\n"

    with data.open_dataset("sample.txt") as second:
        assert second.read() == "alpha\nbeta\n"

    assert resource.opens == 1
    data._read_dataset.cache_clear()
