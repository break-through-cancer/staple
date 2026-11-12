"""Offline tests for the staple_sheets preprocess. No Cirro session or network is used.

The Cirro SDK import is stubbed so the tests run anywhere; the flag logic is pure.
"""

import sys
import types

# Stub the SDK module the preprocess imports, so this test needs no Cirro install.
_helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
_helpers.PreprocessDataset = object
sys.modules.setdefault("cirro", types.ModuleType("cirro"))
sys.modules.setdefault("cirro.helpers", types.ModuleType("cirro.helpers"))
sys.modules["cirro.helpers.preprocess_dataset"] = _helpers

import preprocess  # noqa: E402


class FakeDataset:
    def __init__(self, params):
        self.params = dict(params)
        self.logger = types.SimpleNamespace(info=lambda *_a, **_k: None)

    def remove_param(self, name):
        self.params.pop(name)

    def add_param(self, name, value):
        self.params[name] = value


def test_visium_format_sets_only_visium():
    assert preprocess.format_flags("visium", "") == {"visium": True}


def test_anndata_format_sets_only_anndata():
    assert preprocess.format_flags("anndata", None) == {"anndata": True}


def test_xenium_format_sets_only_xenium():
    assert preprocess.format_flags("xenium", None) == {"xenium": True}


def test_visium_hd_keeps_table_name():
    assert preprocess.format_flags("visium_hd", "square_016um") == {"visium_hd": "square_016um"}


def test_visium_hd_without_table_name_sets_nothing():
    assert preprocess.format_flags("visium_hd", "") == {}


def test_table_name_ignored_for_other_formats():
    assert preprocess.format_flags("visium", "square_016um") == {"visium": True}


def test_prepare_replaces_stale_flags():
    ds = FakeDataset({"format": "anndata", "visium": True, "visium_hd": "square_002um"})
    preprocess.prepare_params(ds)
    assert ds.params == {"format": "anndata", "anndata": True}


def test_prepare_does_not_write_a_samplesheet():
    ds = FakeDataset({"format": "xenium", "input": "from-the-form.csv"})
    preprocess.prepare_params(ds)
    assert ds.params["input"] == "from-the-form.csv"


if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main([__file__, "-q"]))
