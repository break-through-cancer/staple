#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
from urllib.parse import urlparse
from pathlib import Path

SAMPLESHEET_REQUIRED_COLUMNS = ("sample",
                                "data_directory",
                                "expression_profile"
)

# Helper function to check if a string is a URL
def is_url(string):
    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except ValueError:
        return False


def prepare_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    ds.logger.info([ds.params])

    format_name = ds.params.get('format')
    visium_hd_value = ds.params.get('visium_hd')

    for param_name in ('visium', 'anndata', 'xenium', 'visium_hd'):
        if param_name in ds.params:
            ds.remove_param(param_name)

    if format_name == 'visium':
        ds.add_param('visium', True)
    elif format_name == 'anndata':
        ds.add_param('anndata', True)
    elif format_name == 'xenium':
        ds.add_param('xenium', True)
    elif format_name == 'visium_hd' and visium_hd_value:
        ds.add_param('visium_hd', visium_hd_value)
    
    # If the reference_scrna is not a URL, we assume it is a file mask string
    # to look for in the data directory downstream
    if 'reference_scrna' in ds.params and not is_url(ds.params['reference_scrna']):
        ds.params['expression_profile'] = ds.params['reference_scrna']

    samplesheet = samplesheet_from_files(ds)
    ds.logger.info(f'Prepared samplesheet from files: {samplesheet.to_dict()}')
    
    #check is pipeline uses Cirro samplesheet, and if not prepare it from params
    if samplesheet.empty:
        ds.logger.warning("No files found in dataset. Preparing samplesheet from params.")
        samplesheet = samplesheet_from_params(ds)
        if samplesheet.empty:
            raise ValueError("No files found in dataset and unable to prepare samplesheet from params.")
        ds.logger.info("Prepared samplesheet from params.")
    
    # Ensure all required columns are present (populate missing)
    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in samplesheet.columns:
            ds.logger.warning(f"Samplesheet is missing required column '{colname}'. Populating with NaN.")
            samplesheet[colname] = np.nan

    # Save to a file
    samplesheet.to_csv("cirro-samplesheet.csv", index=None)

    # Clear params that we wrote to the samplesheet
    # cleared params will not overload the nextflow.params
    to_remove = []
    for k in ds.params:
        if k in SAMPLESHEET_REQUIRED_COLUMNS:
            ds.logger.info(f"Removing param '{k}' from dataset params as it is now in the samplesheet.")
            to_remove.append(k)

    for k in to_remove:
        ds.remove_param(k)

    ds.add_param("input", "cirro-samplesheet.csv")

    # Log the samplesheet
    ds.logger.info(samplesheet.to_dict())


def samplesheet_from_files(ds):
    pipeline_param_names = [c for c in SAMPLESHEET_REQUIRED_COLUMNS]
    files = ds.files
    
    ds.logger.info(f'found files in ds.files: {files.to_dict()}')

    # Assumes samplesheet associates sample with a file in the sample's root directory
    # Convert s3 link to PosixPath and derive parent; convert back into string
    # Path converts s3:// to s3:/, so revert proper s3 prefix afterwards
    files['data_directory'] = files['file'].apply(lambda x: str(Path(x).parent).replace('s3:/', 's3://'))
    files = files[['sample','data_directory']]

    samplesheet = pd.merge(ds.samplesheet, files, on='sample', how='left')
    

    return samplesheet

def samplesheet_from_params(ds):

    data_params = pd.DataFrame({
        'sample':[x['name'] for x in ds.metadata['inputs']],
        'data_directory': [x['dataPath'] for x in ds.metadata['inputs']]
    })

    return data_params

def main():
    ds = PreprocessDataset.from_running()

    prepare_samplesheet(ds)

    # log
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
