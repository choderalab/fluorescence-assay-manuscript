import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_context("poster")
from matplotlib.colors import LogNorm

#Competitive binding function
def three_component_competitive_binding(Ptot, Ltot, Kd_L, Atot, Kd_A):
    """
    Parameters
    ----------
    Ptot : float
        Total protein concentration
    Ltot : float
        Total tracer(fluorescent) ligand concentration
    Kd_L : float
        Dissociation constant
    Atot : float
        Total competitive ligand concentration
    Kd_A : float
        Dissociation constant
        
    Returns
    -------
    P : float
        Free protein concentration
    L : float
        Free ligand concentration
    A : float
        Free ligand concentration
    PL : float
        Complex concentration
    Kd_L_app : float
        Apparent dissociation constant of L in the presence of A
        
    Usage
    -----
    [P, L, A, PL, Kd_L_app] = three_component_competitive_binding(Ptot, Ltot, Kd_L, Atot, Kd_A)
    """
    Kd_L_app = Kd_L*(1+Atot/Kd_A)                                
    PL = 0.5 * ((Ptot + Ltot + Kd_L_app) - np.sqrt((Ptot + Ltot + Kd_L_app)**2 - 4*Ptot*Ltot))  # complex concentration (uM)
    P = Ptot - PL; # free protein concentration in sample cell after n injections (uM)                                                                                                                                                                                                                          
    L = Ltot - PL; # free tracer ligand concentration in sample cell after n injections (uM)
    A = Atot - PL; # free competitive ligand concentration in sample cell after n injections (uM)
    return [P, L, A, PL, Kd_L_app]

#These are binding affinities to ABL1-nonphosphorylated according to guidetopharmacology.org
Kd_Bos = 0.1e-9 # M
Kd_Gef = 2200e-9 # M
Kd_Ima = 1.1e-9 # M

Ptot = 0.5e-6 # M
Ltot = 20.0e-6 / np.array([10**(float(i)/2.0) for i in range(12)]) # M

concentration_range = [10e-5,10e-6,10e-7,10e-8,10e-9,10e-10,10e-11,10e-12] # M

for i,conc in enumerate(concentration_range):
    [P_gef_ima, L_gef_ima, A_gef_ima, PL_gef_ima, Kd_gef_ima] = three_component_competitive_binding(Ptot, Ltot, Kd_Gef, conc, Kd_Ima)
    if i == 0: 
          competition_grid_gefitinib = PL_gef_ima
    else:
          competition_grid_gefitinib= np.vstack((competition_grid_gefitinib,PL_gef_ima))

#Changing sig figs for Ltot
Ltot_visual = ['%.2g' % a for a in Ltot]
Ltot_visual

plt.figure(figsize=(12, 6))
plt.pcolor(competition_grid_gefitinib, 
           norm=LogNorm(vmin=competition_grid_gefitinib.min(), 
                        vmax=competition_grid_gefitinib.max()), 
           edgecolors='w',
           linewidths=2,
           cmap='PuBu_r')
plt.ticklabel_format(style='plain')
plt.xlabel('fluorescent ligand concentration (M)')
plt.ylabel('non-fluorescent ligand concentration (M)')
plt.xticks(np.arange(0.5, 12.5),Ltot_visual,rotation='vertical');
plt.yticks(np.arange(0.5, 12.5),concentration_range);
plt.ylim((0, len(concentration_range)))
plt.title('Competition Grid for Gefitinib X Imatinib in 0.5 uM Abl')
plt.colorbar(label='complex concentration (M)');
plt.tight_layout()

plt.savefig('competition_grid.png', dpi=500)
plt.savefig('competition_grid.pdf')

plt.clf()

colors = ['cyan','blue','green','yellowgreen','gold','orange','red','maroon']

sns.set(style='white')
sns.set_context('talk')


fig, ax = plt.subplots(figsize=(8,4))
for i in range(len(competition_grid_gefitinib)):
    plt.semilogx(Ltot, competition_grid_gefitinib[i],color=colors[i],label=concentration_range[i]);
handles, labels = ax.get_legend_handles_labels()
plt.legend(reversed(handles), reversed(labels), title='[Ima] (M)');
plt.xlabel('$[L]_T$ (M)',fontsize=16);
plt.yticks([])
plt.xticks(fontsize=16)
plt.xlim((2e-9,5e-5))
plt.ylabel('Fluorescence',fontsize=20);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

plt.savefig('competition_curves.png', dpi=500)
plt.savefig('competition_curves.pdf')
