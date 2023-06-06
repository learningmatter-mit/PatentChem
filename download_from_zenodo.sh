#!/bin/bash

# Run this script from the directory you want to become the `data_dir` in future steps
# Download size is ~75GB - make sure you do this on a machine with enough space!

wget https://zenodo.org/record/7992427/files/20230528_chemical_patents_uspto.tar.gz
tar -xzf 20230528_chemical_patents_uspto.tar.gz
rm 20230528_chemical_patents_uspto.tar.gz
cd 20230528_chemical_patents_uspto
for y in {2001..2023}; do(tar -xzf $y.tar.gz $y); done