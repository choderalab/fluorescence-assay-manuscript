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
mcmc = pymcmodels.run_mcmc(pymc_model,db = 'pickle', dbname = '2nM_Kd_mcmc.pickle')

#Let's define a list of Kd's, which will include our original Kd (2e-9 M):
# Note here that our protein concentration is: Ptot = 1e-9 # M

Kd_list = [
           0.02e-9, # M
           0.2e-9, # 
           2e-9, # M
           20e-9, # M
           200e-9, # M
           2000e-9, # M
           20000e-9, # M
          ]

Kd_label = []

for my_Kd in Kd_list:
    if (my_Kd < 1e-6):
        Kd_label.append("%.3g nM" %(my_Kd/1e-9))
    elif (my_Kd < 1e-3):
        Kd_label.append("%.1f uM" %(my_Kd/1e-6))

# Now let's model our results for all of these Kd's

F_PL_i_list = []
F_L_i_list = []
pymc_model_list = []
mcmc_list = []

for i,Kd in enumerate(Kd_list):
    [L, P, PL] = two_component_binding(Kd, Ptot, Ltot)
    F_PL_i = F_background + ((1400/1e-9)*PL + sigma * np.random.randn(npoints)) + ((.4/1e-8)*L + sigma * np.random.randn(npoints))
    F_L_i = F_background + (.4/1e-8)*Ltot + sigma * np.random.randn(npoints)
    
    pymc_model = pymcmodels.make_model(Ptot, dPstated, Ltot, dLstated,
        top_complex_fluorescence=F_PL_i,
        top_ligand_fluorescence=F_L_i,
        use_primary_inner_filter_correction=True,
        use_secondary_inner_filter_correction=True,
        assay_volume=assay_volume, DG_prior='uniform')
    
    mcmc_list.append(pymcmodels.run_mcmc(pymc_model,db = 'pickle', dbname = 'Kdlist_%s_mcmc.pickle'%i))
    
    F_PL_i_list.append(F_PL_i)
    F_L_i_list.append(F_L_i)
    pymc_model_list.append(pymc_model)
    
# Now let's plot our results

# We need to define a few things to use in the plots

property_name = 'top_complex_fluorescence'

complex_list = []

for pymc_model in pymc_model_list:
    complex_list.append(getattr(pymc_model, property_name))

property_name = 'top_ligand_fluorescence'
ligand = getattr(pymc_model_list[2], property_name)

pale_colors = ['powderblue','cornflowerblue','lightcoral','darkseagreen','palegreen','khaki','lightyellow']
dark_colors = ['cyan','blue','red','green','yellowgreen','goldenrod','yellow']

# Now we'll plot our simulated fluorescence

fig, ax = plt.subplots(figsize=(8,4))

for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::10]:
    plt.semilogx(Ltot, top_ligand_fluorescence_model, marker='.',color='silver', alpha=0.2, hold=True)
for i,this_mcmc in enumerate(mcmc_list):
    for top_complex_fluorescence_model in this_mcmc.top_complex_fluorescence_model.trace()[::10]:
        plt.semilogx(Ltot, top_complex_fluorescence_model, marker='.',color=pale_colors[i], alpha=0.2, hold=True)

for i,complex in enumerate(complex_list):
    plt.semilogx(Ltot, complex.value, color=dark_colors[i], marker = 'o', linestyle= 'None', label='%s'%Kd_label[i], hold=True)
plt.semilogx(Ltot, ligand.value, 'ko',label='ligand', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlabel('$[L]_T$ (M)',fontsize=16);
plt.yticks([])
plt.xticks(fontsize=16)
plt.xlim((2e-11,5e-5))
plt.ylabel('Fluorescence',fontsize=20);

legend=plt.legend(fontsize=16,title=r'Kd',loc='center left', bbox_to_anchor=(1, 0.5));
plt.setp(legend.get_title(),fontsize='xx-large')
plt.setp(legend.get_title(),weight='bold')

plt.tight_layout()

plt.savefig('simulated_fluorescence_P_1nM.png', dpi=500, bbox_inches='tight')
plt.savefig('simulated_fluorescence_P_1nM.pdf', bbox_inches='tight')

plt.close()

# Now we'll plot our deltaG histograms

fig, ax = plt.subplots(figsize=(8,4))

binBoundaries = np.linspace(-33,-10,70)

for i,this_mcmc in enumerate(mcmc_list):
    plt.hist(this_mcmc.DeltaG.trace(),color=dark_colors[i],bins=binBoundaries, hold=True,edgecolor='white')
    plt.axvline(x=np.log(Kd_list[i]),color=dark_colors[i],label='%s'%Kd_label[i], hold=True)
plt.axvline(x=np.log(Ptot[0]),color='0.7',linestyle='--',label='[P] = 1 nM', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.yticks([])
plt.xticks(fontsize=16)
plt.xlim((-27,-9))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16);
plt.ylabel('$P(\Delta G)$',fontsize=20);

legend=plt.legend(fontsize=16,title=r'Kd',loc='center left', bbox_to_anchor=(1, 0.5));
plt.setp(legend.get_title(),fontsize='xx-large')
plt.setp(legend.get_title(),weight='bold')

plt.tight_layout()

plt.savefig('simulated_fluorescence_delG_result.png', dpi=500, bbox_inches='tight')
plt.savefig('simulated_fluorescence_delG_result.pdf', bbox_inches='tight')

plt.close()

# Now we'll plot our Kd histograms

fig, ax = plt.subplots(figsize=(8,4))

kd_binBoundaries = np.exp(np.arange(-33,0,0.5))

for i,this_mcmc in enumerate(mcmc_list):
    plt.hist(np.exp(this_mcmc.DeltaG.trace()),color=dark_colors[i],bins=kd_binBoundaries, hold=True,edgecolor='white')
    plt.axvline(x=Kd_list[i],color=dark_colors[i],label='%s'%Kd_label[i], hold=True)
plt.axvline(x=Ptot[0],color='0.7',linestyle='--',label='[P] = 1 nM', hold=True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xscale('log')
plt.xlim((1e-12,1e-4))
plt.yticks([])
plt.xticks(fontsize=16)
plt.xlabel('$K_{d}$ ($M$)',fontsize=16);
plt.ylabel('$P(K_{d})$',fontsize=20);

legend=plt.legend(fontsize=16,title=r'Kd',loc='center left', bbox_to_anchor=(1, 0.5));
plt.setp(legend.get_title(),fontsize='xx-large')
plt.setp(legend.get_title(),weight='bold')

plt.tight_layout()

plt.savefig('simulated_fluorescence_Kd_result.png', dpi=500, bbox_inches='tight')
plt.savefig('simulated_fluorescence_Kd_result.pdf', bbox_inches='tight')

plt.close()

#Now let's calculate the CV and relative bias of these results

CV_list = []
bias_list = []

for i,this_mcmc in enumerate(mcmc_list):
    CV_list.append(np.std(np.exp(this_mcmc.DeltaG.trace()))/np.mean(np.exp(this_mcmc.DeltaG.trace())) * 100)
    bias_list.append((np.mean(np.exp(this_mcmc.DeltaG.trace())) - Kd_list[i]) / Kd_list[i])

f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(8,4.5))

for i,CV in enumerate(CV_list):
    ax1.semilogx(Kd_list[i], CV, dark_colors[i], marker = 'o', linestyle='None')

ax1.set_ylabel('$CV (\%)$',fontsize=20);
ax1.set_ylim((-20,200));

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

for i,bias in enumerate(bias_list):
    ax2.semilogx(Kd_list[i], bias*100, dark_colors[i], marker = 'o', linestyle='None')

ax2.set_ylabel('$RB (\%)$',fontsize=20);
ax2.set_ylim((-80,30));

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

f.subplots_adjust(hspace=0)
plt.xlabel('$K_{d}$ ($M$)',fontsize=20);
plt.xticks(fontsize=16);
plt.yticks(fontsize=16);

plt.tight_layout()

plt.savefig('simulated_fluorescence_CV_bias.png', dpi=500, bbox_inches='tight')
plt.savefig('simulated_fluorescence_CV_bias.pdf', bbox_inches='tight')

plt.close()


