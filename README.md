# PatentChem
[//]: # (Badges)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/612404368.svg)](https://zenodo.org/badge/latestdoi/612404368)


This code downloads the weekly bulk files from the [USPTO](https://developer.uspto.gov/product/patent-grant-dataxml), selects the patents relevant to chemistry, and queries the chemistry patent claims sections for given keywords to find molecule SMILES related to those keywords.

## Installation
```
conda env create -f environment.yml
```
Alternatively, the following may be faster:
```
conda install -c conda-forge mamba
mamba env create -f environment.yml
```

## Example Usage
The scripts should be run in the following order. Steps 1 and 2 only need to be run once, then steps 3 and 4 can be run as many times as desired using different query terms.


### 1. Download
Download all USPTO patents for the given years: 

`python download.py --years 2021 --data_dir /data/patents_data/`

Additional options:
* `--force_redownload`: download files even if they already exist
* `--no_uncompress`: do not uncompress the `*.zip` and `*.tar` files after selecting the chemistry-related ones (to save space if downloading all years at once)
* `--remove_compressed`: remove the original .tar and .zip files for the weekly releases after decompressing them

### 2. Select
Select the patents related to chemistry based on the presence of CDX files:

`python select_chem.py --years 2021 --data_dir /data/patents_data/`

Additional options:
* `--remove_compressed`: remove all original `*.zip` and `*.tar` files after unzipping those relevant to chemistry

### 3. Search
Provide a list of query terms to search for in the claims sections of the chemistry-related patents:

`python search.py --years 2021 --data_dir /data/patents_data/ --naming opd --subject_list "organic electronic" "photodiode" "organic photovoltaic"`

### 4. Clean
Output a yaml file with the SMILES strings found related to the query terms and all of the patents related to each SMILES. Also output a text file with all unique SMILES strings from the query.

`python clean.py --years 2021 --naming opd`

Additional options:
* `--output_dir`: parent directory of `--naming` if different from `.`
* `--charged_only`: include only charged molecules in output
* `--neutral_only`: include only neutral molecules in output
* `--mw_min`: include only molecules with molecular weight greater than this
* `--mw_max`: include only molecules with molecular weight less than this

### Other Useful Commands
* If you ran `download.py` with `--no_uncompress` and now want to uncompress the files: `python utils/unzip.py --years 2021 --data_dir /data/patents_data/`
* If you've run `select_chem.py` but want to reset the data directory to its state after downloading (with unzipped files): `python utils/reset.py --years 2021 --data_dir /data/patents_data/`
  * This assumes you didn't use `--remove_compressed` with `download.py`. If you used `--remove_compressed`, you'll have to redo the `download.py` step.
  * Add `--remove_uncompressed_only` if you want to reset to a state after downloading with `--no_uncompress` (all files other than the original tar/zip files will be deleted)

## Download File Sizes

The following file sizes are taken from the [USPTO Bulk Data Storage System](https://bulkdata.uspto.gov) using the URLs `https://bulkdata.uspto.gov/data/patent/grant/redbook/<YEAR>/` and converting from bytes to GB. The 2023 file size is as of 23 May 2023. The download speed seems to be restricted to ~5-10MB/s, which means downloading the full set could require > 4 days if done in series. Alternatively, you can run multiple downloads in parallel. Either way, we recommend starting the downloads in a [tmux](https://github.com/tmux/tmux/wiki) session to run in the background. Note that these file sizes are for the *compressed* (zip or tar) files, so the total space required to store this data is larger than what is reported below. Use caution to avoid filling your hard drive to capacity.

| **Year**  | **File Size** | **Units** |
|-------|-----------|-------|
| 2001  | 37.07     | GB    |
| 2002  | 37.23     | GB    |
| 2003  | 39.20     | GB    |
| 2004  | 38.95     | GB    |
| 2005  | 30.67     | GB    |
| 2006  | 39.90     | GB    |
| 2007  | 39.98     | GB    |
| 2008  | 42.08     | GB    |
| 2009  | 43.82     | GB    |
| 2010  | 61.65     | GB    |
| 2011  | 68.24     | GB    |
| 2012  | 85.42     | GB    |
| 2013  | 96.55     | GB    |
| 2014  | 111.26    | GB    |
| 2015  | 119.98    | GB    |
| 2016  | 117.42    | GB    |
| 2017  | 124.10    | GB    |
| 2018  | 148.02    | GB    |
| 2019  | 182.28    | GB    |
| 2020  | 192.45    | GB    |
| 2021  | 184.41    | GB    |
| 2022  | 184.1     | GB    |
| 2023  | 80.3+     | GB    |
| **Total** | **2.1**      | **TB**    |

## Notes
* For simplicity, `clean.py` replaces "*" and "Y" substituents from Markush structures with ethyl and O groups. This might not be appropriate for your applications.
<!-- TODO: is the above relevant? -->

## Citation
If you use this code, please cite the following [manuscript](https://arxiv.org/abs/2303.08272):
```
@article{patents-generative2023,
  title={Automated patent extraction powers generative modeling in focused chemical spaces},
  author={Subramanian, Akshay and Greenman, Kevin P. and Gervaix, Alexis and Yang, Tzuhsiung and G{\'{o}}mez-Bombarelli, Rafael},
  journal = {arXiv preprint arXiv:2303.08272},
  doi = {10.48550/ARXIV.2303.08272},
  year={2023}
}
```
