#!/bin/bash

# This script automatically generates pdf file and csv file summaries
# of bayesian predictions of delG from xml data files.
# REQUIREMENTS: assaytools must be installed

# Run for first spectra assay.

quickmodel --inputs 'inputs_first_spectra'

mkdir first_spectra

mv *.png first_spectra
mv *.json first_spectra
mv *.pickle first_spectra

# Run for first p38 spectra assay.

quickmodel --inputs 'inputs_p38_spectra1'

mkdir first_spectra_p38

mv *.png first_spectra_p38
mv *.json first_spectra_p38
mv *.pickle first_spectra_p38

# Run for spectra 2gains assay.

quickmodel --inputs 'inputs_spectra_2gains'

mkdir spectra_2gains

mv *.png spectra_2gains
mv *.json spectra_2gains
mv *.pickle spectra_2gains

#Run for p38 spectra analysis for other ligands.

quickmodel --inputs 'inputs_p38_spectra_otherligs'

mkdir first_p38_otherligs

mv *.png first_p38_otherligs
mv *.json first_p38_otherligs
mv *.pickle first_p38_otherligs

# Run for CAII spectra analysis for other ligands.

quickmodel --inputs 'inputs_CAII_spectra_otherligs'

mkdir first_CAII_otherligs

mv *.pdf first_CAII_otherligs
mv *.json first_CAII_otherligs
mv *.pickle first_CAII_otherligs

# Run for singlet assay with 4 ligands and DMSO backfill

quickmodel --inputs 'inputs_DMSO' 

mkdir DMSO-backfill

mv *.png DMSO-backfill
mv *.json DMSO-backfill
mv *.pickle DMSO-backfill

# Run for singlet assay with 8 different ligands

quickmodel --inputs 'inputs_8ligs'

mkdir 8ligs

mv *.png 8ligs
mv *.json 8ligs
mv *.pickle 8ligs

# Run for singlet assay with 8 different ligands and 0.5uM of protein

quickmodel --inputs 'inputs_8ligs_0.5_protein_concentration'

mkdir 8ligs_0.5_protein_concentration

mv *.png 8ligs_0.5_protein_concentration
mv *.json 8ligs_0.5_protein_concentration
mv *.pickle 8ligs_0.5_protein_concentration

# Run for singlet assay with 8 different ligands and 0.25uM of protein

quickmodel --inputs 'inputs_8ligs_0.25_protein_concentration'

mkdir 8ligs_0.25_protein_concentration

mv *.png 8ligs_0.25_protein_concentration
mv *.json 8ligs_0.25_protein_concentration
mv *.pickle 8ligs_0.25_protein_concentration

# Run for singlet assay with 8 different ligands and 0.125uM of protein

quickmodel --inputs 'inputs_8ligs_0.125_protein_concentration'

mkdir 8ligs_0.125_protein_concentration

mv *.png 8ligs_0.125_protein_concentration
mv *.json 8ligs_0.125_protein_concentration
mv *.pickle 8ligs_0.125_protein_concentration
