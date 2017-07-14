# MIXING BY DIFFUSION EXPERIMENT

This experiment was designed to test if diffusion can help mixing in-384 well plates by letting plate incubate and reading fluorescence signal at time intervals.
 * p38 binding assay with 8 fluorescent ligands (different than default set)
 * After standard protocol, plate was left for diffusion for indicated time and multiple fluorescence reads were performed in time.
 * Performed on June 15-16, 2017
 * 384 well plate, assay volume 50 uL
 * protein: p38 kinase
 * ligands: Bosutinib [AB], Bosutinib Isomer [CD], Erlotinib HCl [EF], Gefitinib [GH], Lapatinib [IJ], Ponatinib [KL], Vandetanib (new stock) [MN], and Vandetanib (old stock) [OP]
 * Without DMSO backfill
 * Time points: 1, 2, 18 and 26 hours after standard plate preparation.

### Quickmodel Analysis
```
quickmodel --inputs 'inputs_p38_singlet_diffusion' --type 'singlet' --nsamples 100000
```

### Manifest
- `20170615_1st_read_1h\` - Contains quickmodel analysis results for spectral data collected after 1 hour waiting time after preparation.
- `20170615_2nd_read_2h\` - Contains quickmodel analysis results for spectral data collected after 2 hour waiting time after preparation.
- `20170616_3rd_read_18h\` - Contains quickmodel analysis results for spectral data collected after 18 hour waiting time after preparation.
- `20170616_4th_read_26h\` - quickmodel analysis results for spectral data collected after 18 hour waiting time after preparation.    
- *inputs_p38_singlet_diffusion.py* - Inputs file for Quickmodel which contains experimental details.
- *run-torque.sh* - Job submission script for computer cluster
- `results-summary\` - Contains grid layout of some figures for easier comparison of fluorescence signal and binding free energy change in time.
