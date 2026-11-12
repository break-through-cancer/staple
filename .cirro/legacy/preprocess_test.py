import os.path
import pandas as pd

from cirro.helpers.preprocess_dataset import PreprocessDataset

import preprocess


def test_preprocess(mocker):
    mock_dataset = mocker.Mock(spec=PreprocessDataset)

    logger = mocker.Mock()

    mock_dataset.logger = logger

    mock_dataset.params = {
        "sample": "sample",
        "outdir": "s3://bucket/output",
        "data_directory": "s3://bucket/data",
        "expression_profile": "expression_profile"
    }

    try:
        with mocker.patch("cirro.helpers.preprocess_dataset.PreprocessDataset.from_running", return_value=mock_dataset):
            preprocess.main()

        assert os.path.exists("samplesheet.csv")

        df = pd.read_csv("samplesheet.csv")
        for colname in preprocess.SAMPLESHEET_REQUIRED_COLUMNS:
            assert colname in df.columns
    finally:
        os.remove("samplesheet.csv")
