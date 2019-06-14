# Activate conda environment (Python 3.6)
source activate assaytools_20190220_43a6d53

# Save environment dependencies
conda list --export > requirements.txt

# Plots of 4 ligands with good DeltaG fits
python generate_p38_curve_delG_hist_for_20171119_single_well_exp.py

# Plots of 8 ligands
python generate_p38_curve_delG_hist_for_20171119_single_well_8ligs.py
