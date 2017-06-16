import cPickle
import numpy as np

import seaborn as sns
cols = sns.color_palette("deep", 3)
sns.set(style='white')
sns.set_context('talk')

from pymbar import timeseries
import matplotlib.pyplot as plt

#from assaytools.assay_simulator import AssaySimulator, predict_assay_error
from assay_simulator import AssaySimulator, predict_assay_error


### IMPORT DATA

data_file = '../../analysis/bayes/first_spectra_p38/p38-Bosutinib-AB_mcmc-2016-09-20 18:56.pickle'
with open(r'%s'%data_file,'rb') as my_file:
    data = cPickle.load(my_file)

DeltaG_samples = data['DeltaG'][0]
(t_equil, g, N_eff) = timeseries.detectEquilibration(DeltaG_samples, fast=True, nskip=1)

# Seeing the mean of the free energy
mean_DeltaG = np.mean(DeltaG_samples[t_equil:])
print('Expected free energy = {0} (thermal units)'.format(mean_DeltaG))

L_total = 10 ** (np.linspace(-10, -5, num=12))
P_total = 1.0e-6* np.ones([12],np.float64)

### SIMULATE ASSAY BASED ON THIS DATA, BUT USING THREE PROTEIN CONCENTRATIONS

nsamples = 1000
protein_concentrations = np.array([0.5E-6, 1E-6, 5E-6])

estimates_set = []
intensities_set = []

for i in range(len(protein_concentrations)):
    sim_assay = AssaySimulator(pymc_data=data, l_total=L_total, p_total=protein_concentrations[i]* np.ones([12],np.float64), inner_filter=False)
    sim_assay.set_mean_parameters(pymc_data=data, t_equil=t_equil)
    estimates, intensities = sim_assay.generate_deltag_estimates(nsamples)
    estimates_set.append(estimates)
    intensities_set.append(intensities)

### PLOT OUR ESTIMATES FROM THESE THREE EXPERIMENTS

fig, ax = plt.subplots(figsize=(6,4))

for i in range(len(protein_concentrations)):
    hist, edges = np.histogram(estimates_set[i], bins=20, normed=True)
    hist = hist/hist.max()
    centers = edges[0:-1] + np.diff(edges) / 2.0
    if i !=2:
        ax.plot(centers, hist, lw=3, ls='--', color=cols[i], label="{} $\mu$M".format(protein_concentrations[i]*1.0E6))
    else:   
        ax.bar(centers, hist, width=0.01, color=cols[i], label="{} $\mu$M".format(protein_concentrations[i]*1.0E6))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
        
#plt.title('Predicted posterior free energy density at \n different protein concentrations', fontsize=14);
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Binding free energy (unitless)', fontsize=16);
plt.ylabel('Density relative to maximum', fontsize=16);
plt.legend(fontsize=15);
plt.tight_layout()

plt.savefig('DeltaG_Estimates_Three_Protein.png', dpi=500)

### PLOT OUR SIMULATED EXPERIMENTS THAT LEAD TO THESE ESTIMATES

fig, ax = plt.subplots(figsize=(6,4))

for k in range(len(intensities_set)):
    for i in range(intensities_set[k].shape[0]):
        if i == 0:
            ax.semilogx(L_total,intensities_set[k][i,:], color=cols[k], lw=3, label="{} $\mu$M".format(protein_concentrations[k]*1.0E6))
        else:
            ax.semilogx(L_total,intensities_set[k][i,:], color=cols[k], alpha=0.01, lw=1)
            
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#plt.title('Simulated fluorescence data for p38-Bosutinib', fontsize=16);
plt.xticks(fontsize=15);
plt.yticks(fontsize=15);
plt.xlabel('Ligand concentration (M)', fontsize=16);
plt.ylabel('Fluorescence', fontsize=16);
plt.legend(fontsize=15);
plt.tight_layout()

plt.savefig('Intensities_Three_Protein.png', dpi=500)



