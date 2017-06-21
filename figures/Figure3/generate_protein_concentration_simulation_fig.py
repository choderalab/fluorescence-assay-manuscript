import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize
import seaborn as sns

sns.set(style='white')
sns.set_context('talk')

# We define a Kd,
Kd = 2e-9 # M

# a protein concentration,
Ptot = 1e-9 * np.ones([12],np.float64) # M

# and a gradient of ligand concentrations for our experiment.
Ltot = 20.0e-6 / np.array([10**(float(i)/2.0) for i in range(12)]) # M

def two_component_binding(Kd, Ptot, Ltot):
    """
    Parameters
    ----------
    Kd : float
        Dissociation constant
    Ptot : float
        Total protein concentration
    Ltot : float
        Total ligand concentration
        
    Returns
    -------
    P : float
        Free protein concentration
    L : float
        Free ligand concentration
    PL : float
        Complex concentration
    """
                                    
    PL = 0.5 * ((Ptot + Ltot + Kd) - np.sqrt((Ptot + Ltot + Kd)**2 - 4*Ptot*Ltot))  # complex concentration (uM)
    P = Ptot - PL; # free protein concentration in sample cell after n injections (uM)                                                                                                                                                                                                                          
    L = Ltot - PL; # free ligand concentration in sample cell after n injections (uM)  
    
    return [P, L, PL]

[L, P, PL] = two_component_binding(Kd, Ptot, Ltot)

# Making max 1400 relative fluorescence units, and scaling all of PL (complex concentration) 
# to that, adding some random noise
npoints = len(Ltot)
sigma = 10.0 # size of noise
F_PL_i = (1400/1e-9)*PL + sigma * np.random.randn(npoints)

#Let's add an F_background just so we don't ever go below zero
F_background = 40
#We also need to model fluorescence for our ligand
F_L_i = F_background + (.4/1e-8)*Ltot + sigma * np.random.randn(npoints)

#Let's also add these to our complex fluorescence readout
F_PL_i = F_background + ((1400/1e-9)*PL + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L + sigma * np.random.randn(npoints))


# We know errors from our pipetting instruments.
P_error = 0.15
L_error = 0.08

assay_volume = 100e-6 # assay volume, L

dPstated = P_error * Ptot
dLstated = L_error * Ltot


from assaytools import pymcmodels
pymc_model = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
    top_complex_fluorescence=F_PL_i,
    top_ligand_fluorescence=F_L_i,
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=assay_volume, DG_prior='uniform')
mcmc = pymcmodels.run_mcmc(pymc_model)

# For this we're just going to define a range of Kd's, with some qualitative variable names.
# Note here that our protein concentration is: Ptot = 1e-9 # M

Kd_high = 0.2e-9 # M
Kd_really_high = 0.02e-9 # M
Kd_low = 20e-9 # M
Kd_really_low = 200e-9 # M


# Will need an F_PL_i for all of these
[L_high, P_high, PL_high] = two_component_binding(Kd_high, Ptot, Ltot)
[L_really_high, P_really_high, PL_really_high] = two_component_binding(Kd_really_high, Ptot, Ltot)
[L_low, P_low, PL_low] = two_component_binding(Kd_low, Ptot, Ltot)
[L_really_low, P_really_low, PL_really_low] = two_component_binding(Kd_really_low, Ptot, Ltot)

F_PL_high_i = F_background + ((1400/1e-9)*PL_high + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L_high + sigma * np.random.randn(npoints))
F_PL_really_high_i = F_background + ((1400/1e-9)*PL_really_high + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L_really_high + sigma * np.random.randn(npoints))
F_PL_low_i = F_background + ((1400/1e-9)*PL_low + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L_low + sigma * np.random.randn(npoints))
F_PL_really_low_i = F_background + ((1400/1e-9)*PL_really_low + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L_really_low + sigma * np.random.randn(npoints))

pymc_model_high = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
    top_complex_fluorescence=F_PL_high_i,
    top_ligand_fluorescence=F_L_i,
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=assay_volume, DG_prior='uniform')
mcmc_high = pymcmodels.run_mcmc(pymc_model_high)

pymc_model_really_high = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
    top_complex_fluorescence=F_PL_really_high_i,
    top_ligand_fluorescence=F_L_i,
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=assay_volume, DG_prior='uniform')
mcmc_really_high = pymcmodels.run_mcmc(pymc_model_really_high)

pymc_model_low = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
    top_complex_fluorescence=F_PL_low_i,
    top_ligand_fluorescence=F_L_i,
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=assay_volume, DG_prior='uniform')
mcmc_low = pymcmodels.run_mcmc(pymc_model_low)

pymc_model_really_low = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
    top_complex_fluorescence=F_PL_really_low_i,
    top_ligand_fluorescence=F_L_i,
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=assay_volume, DG_prior='uniform')
mcmc_really_low = pymcmodels.run_mcmc(pymc_model_really_low)

#Now let's plot our results:

fig, ax = plt.subplots(figsize=(9,4))

property_name = 'top_complex_fluorescence'
complex = getattr(pymc_model, property_name)
complex_high = getattr(pymc_model_high, property_name)
complex_really_high = getattr(pymc_model_really_high, property_name)
complex_low = getattr(pymc_model_low, property_name)
complex_really_low = getattr(pymc_model_really_low, property_name)
property_name = 'top_ligand_fluorescence'
ligand = getattr(pymc_model, property_name)

for top_complex_fluorescence_model in mcmc.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color='lightcoral', alpha=0.2, hold=True)
for top_complex_fluorescence_model in mcmc_high.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color='cornflowerblue', alpha=0.2, hold=True)
for top_complex_fluorescence_model in mcmc_really_high.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color='powderblue', alpha=0.2, hold=True)
for top_complex_fluorescence_model in mcmc_low.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color='darkseagreen', alpha=0.2, hold=True)
for top_complex_fluorescence_model in mcmc_really_low.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color='palegreen', alpha=0.2, hold=True)
for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_ligand_fluorescence_model, marker='.',color='silver', alpha=0.2, hold=True)
    
plt.semilogx(Ltot, complex_really_high.value, color='cyan',marker='o',linestyle='None',label='0.02 nM', hold=True)
plt.semilogx(Ltot, complex_high.value, 'bo',label='0.2 nM', hold=True)
plt.semilogx(Ltot, complex.value, 'ro',label='2 nM', hold=True)
plt.semilogx(Ltot, complex_low.value, color='green',marker='o',linestyle='None',label='20 nM', hold=True)
plt.semilogx(Ltot, complex_really_low.value, color='yellowgreen',marker='o',linestyle='None',label='200 nM', hold=True)
plt.semilogx(Ltot, ligand.value, 'ko',label='ligand', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlabel('$[L]_T$ (M)',fontsize=16);
plt.yticks([])
plt.xticks(fontsize=16)
plt.xlim((1e-11,5e-5))
plt.ylabel('Fluorescence',fontsize=20);
plt.legend(fontsize=16);

plt.tight_layout()

plt.savefig('simulated_fluorescence_P_1nM.png', dpi=500)
plt.savefig('simulated_fluorescence_P_1nM.pdf')

fig, ax = plt.subplots(figsize=(9,4))

binBoundaries = np.linspace(-33,-14,70)

plt.hist(mcmc_really_high.DeltaG.trace(),color='cyan',bins=binBoundaries, hold=True,edgecolor='white')
plt.axvline(x=np.log(Kd_really_high),color='cyan',label='0.02 nM', hold=True)
plt.hist(mcmc_high.DeltaG.trace(),color='blue',bins=binBoundaries, hold=True,edgecolor='white')
plt.axvline(x=np.log(Kd_high),color='blue',label='0.2 nM', hold=True)
plt.hist(mcmc.DeltaG.trace(),color='red',bins=binBoundaries, hold=True,edgecolor='white')
plt.axvline(x=np.log(Kd),color='red',label='2 nM', hold=True)
plt.hist(mcmc_low.DeltaG.trace(),bins=binBoundaries,color='green', hold=True,edgecolor='white')
plt.axvline(x=np.log(Kd_low),color='green',label='20 nM', hold=True)
plt.hist(mcmc_really_low.DeltaG.trace(),bins=binBoundaries,color='yellowgreen', hold=True,edgecolor='white')
plt.axvline(x=np.log(Kd_really_low),color='yellowgreen',label='200 nM', hold=True)
plt.axvline(x=np.log(Ptot[0]),color='0.7',linestyle='--',label='log(protein conc)', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.yticks([])
plt.xticks(fontsize=16)
plt.xlim((-33,-14))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16);
plt.ylabel('$P(\Delta G)$',fontsize=20);
plt.legend(fontsize=16);

plt.tight_layout()

plt.savefig('simulated_fluorescence_delG_result.png', dpi=500)
plt.savefig('simulated_fluorescence_delG_result.pdf')

fig, ax = plt.subplots(figsize=(9,4))

kd_binBoundaries = np.exp(np.arange(-33,-14,0.5))

plt.hist(np.exp(mcmc_really_high.DeltaG.trace()),color='cyan',bins=kd_binBoundaries, edgecolor='white')
plt.axvline(x=Kd_really_high,color='cyan',label='0.02 nM')
plt.hist(np.exp(mcmc_high.DeltaG.trace()),color='blue',bins=kd_binBoundaries, hold=True,edgecolor='white')
plt.axvline(x=Kd_high,color='blue',label='0.2 nM', hold=True)
plt.hist(np.exp(mcmc.DeltaG.trace()),color='red',bins=kd_binBoundaries, hold=True,edgecolor='white')
plt.axvline(x=Kd,color='red',label='2 nM', hold=True)
plt.hist(np.exp(mcmc_low.DeltaG.trace()),bins=kd_binBoundaries,color='green', hold=True,edgecolor='white')
plt.axvline(x=Kd_low,color='green',label='20 nM', hold=True)
plt.hist(np.exp(mcmc_really_low.DeltaG.trace()),bins=kd_binBoundaries,color='yellowgreen', hold=True,edgecolor='white')
plt.axvline(x=Kd_really_low,color='yellowgreen',label='200 nM', hold=True)
plt.axvline(x=Ptot[0],color='0.7',linestyle='--',label='protein conc', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xscale('log')
plt.xlim((1e-13,1e-6))
plt.yticks([])
plt.xticks(fontsize=16)
plt.xlabel('$K_{d}$ ($M$)',fontsize=16);
plt.ylabel('$P(K_{d})$',fontsize=20);
plt.legend(fontsize=16);

plt.tight_layout()

plt.savefig('simulated_fluorescence_Kd_result.png', dpi=500)
plt.savefig('simulated_fluorescence_Kd_result.pdf')







