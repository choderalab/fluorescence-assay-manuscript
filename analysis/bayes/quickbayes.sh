#!/bin/bash

# This script automatically generates pdf file and csv file summaries
# of bayesian predictions of delG from xml data files.
# REQUIREMENTS: assaytools must be installed

# make sure the python path includes local dir
# use env PYTHONPATH="./" before quickmodel

# Run for first spectra assay.

env PYTHONPATH="./" quickmodel --inputs 'inputs_first_spectra' --type 'spectra'

mkdir first_spectra

mv *.pdf first_spectra
mv *.npy first_spectra
mv *.json first_spectra
mv *.pickle first_spectra

# Run for first p38 spectra assay.

env PYTHONPATH="./" quickmodel --inputs 'inputs_p38_spectra1' --type 'spectra'

mkdir first_spectra_p38

mv *.pdf first_spectra_p38
mv *.npy first_spectra_p38
mv *.json first_spectra_p38
mv *.pickle first_spectra_p38

# Run for spectra 2gains assay.

env PYTHONPATH="./" quickmodel --inputs 'inputs_spectra_2gains' --type 'spectra'

mkdir spectra_2gains

mv *.pdf spectra_2gains
mv *.npy spectra_2gains
mv *.json spectra_2gains
mv *.pickle spectra_2gains

# Run for singlet assay with 4 ligands and DMSO backfill

env PYTHONPATH="./" quickmodel --inputs 'inputs_DMSO' 

mkdir DMSO-backfill

mv *.png DMSO-backfill
mv *.npy DMSO-backfill
mv *.json DMSO-backfill
mv *.pickle DMSO-backfill

# Run for singlet assay with 8 different ligands

env PYTHONPATH="./" quickmodel --inputs 'inputs_8ligs'

mkdir 8ligs

mv *.png 8ligs
mv *.npy 8ligs
mv *.json 8ligs
mv *.pickle 8ligs
