#!/bin/bash

# This script automatically generates pdf file and csv file summaries
# of bayesian predictions of delG from xml data files.
# REQUIREMENTS: assaytools must be installed

#This is currently incomplete and only analyzes a small selection of the available data

python inputs.py

quickmodel "../../data/singlet/DMSO-backfill"

mv *.png DMSO-backfill
mv *.npy DMSO-backfill
mv *.json DMSO-backfill
mv *.pickle DMSO-backfill
