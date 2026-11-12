#!/usr/bin/env python3
"""Preprocess for the staple_sheets pipeline.

The samplesheet is a structured input: Cirro writes it from the table on the run form
and passes its path to the workflow as `input`. This script does not build it.

What it does: the run form has one `format` choice, but the workflow branches on one
flag per format (visium, anndata, xenium, visium_hd). Those flags are derived here, the
same way the legacy staple preprocess derived them.
"""

from cirro.helpers.preprocess_dataset import PreprocessDataset

FORMAT_FLAGS = ("visium", "anndata", "xenium", "visium_hd")


def format_flags(format_name, visium_hd_value):
    """The per-format flags for one run form choice.

    Only truthy flags are returned: the workflow treats an absent flag and a false one
    the same way, so absent is what a run without that format carries. A Visium HD table
    name only counts when the format is visium_hd.
    """
    flags = {
        "visium": format_name == "visium",
        "anndata": format_name == "anndata",
        "xenium": format_name == "xenium",
        "visium_hd": visium_hd_value if format_name == "visium_hd" and visium_hd_value else None,
    }
    return {name: value for name, value in flags.items() if value}


def prepare_params(ds: PreprocessDataset) -> None:
    flags = format_flags(ds.params.get("format"), ds.params.get("visium_hd"))

    # Clear any flag already in the params, so only the one this run's format sets remains.
    for name in FORMAT_FLAGS:
        if name in ds.params:
            ds.remove_param(name)

    for name, value in flags.items():
        ds.add_param(name, value)

    ds.logger.info(f"Format flags for format {ds.params.get('format')!r}: {flags}")


def main():
    ds = PreprocessDataset.from_running()
    prepare_params(ds)


if __name__ == "__main__":
    main()
