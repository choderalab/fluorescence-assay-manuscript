#!/bin/bash
#BSUB -J Erl_Perror_0.15_100k
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=8]
#BSUB -W 36:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

echo "LS_SUBCWD" 
echo $LS_SUBCWD
echo "PWD" 
echo $PWD

# make sure job’s submission CWD path can be accessed on the job’s execution host 
#cd $LS_SUBCWD
cd $PWD

# Activate the environment where Assaytools is installed
source activate assaytools

# Export current conda environment requirements
conda list --export > requirements.txt

# Run my program
echo “Started running Assaytools...”
quickmodel --inputs inputs_spectra_Erlotinib --type spectra --nsamples 100000
echo “Done.”
