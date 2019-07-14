# We are interested in making a plot to show that DMSO backfill has little effect on final assay results

# First I will use the inputs.py format to import and parse the xml and overlay results with and without 
#   DMSO backfill for bosutinib isomer

import numpy as np
from glob import glob
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

inputs = {
    'xml_file_path' :  "../../data/single-wavelength/DMSO-backfill/",
    'file_set'      :  {'p38': glob( "../../data/single-wavelength/DMSO-backfill/*.xml")},
    'section'       :  '280_480_TOP_120',
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'Lstated'       :  np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Pstated'       :  0.5e-6 * np.ones([24],np.float64), # protein concentration, M
    'assay_volume'  :  50e-6, # assay volume, L
    'well_area'     :  0.1369, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }
	
from assaytools import parser

[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)  

DMSO_backfill = 'p38-Bosutinib Isomer-CD'
No_DMSO_backfill = 'p38-Bosutinib Isomer-KL'

fig, ax = plt.subplots(figsize=(16,4))

plt.subplot(121)
plt.semilogx(inputs['Lstated'], complex_fluorescence[DMSO_backfill], 'bo',label='with DMSO backfill (complex)')
plt.semilogx(inputs['Lstated'], ligand_fluorescence[DMSO_backfill], marker='o',color='lightblue',linestyle='None',label='with DMSO backfill (ligand)')
plt.semilogx(inputs['Lstated'], complex_fluorescence[No_DMSO_backfill], 'go',linestyle='None',label='without DMSO backfill (complex)')
plt.semilogx(inputs['Lstated'], ligand_fluorescence[No_DMSO_backfill], marker='o',color='palegreen',linestyle='None',label='without DMSO backfill (ligand)')
plt.xlabel('$[L]_T$ (M)');
plt.ylabel('fluorescence units');
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc=0);


# Then I will import the pickle file for the analysis that was already run to plot the histograms

import cPickle

Bsi_backfill_file = '../../analysis/bayes/DMSO-backfill/Bosutinib Isomer-CD_mcmc-2017-04-13 00:31.pickle'
Bsi_no_backfill_file = '../../analysis/bayes/DMSO-backfill/Bosutinib Isomer-KL_mcmc-2017-04-13 01:05.pickle'

with open(r'%s'%Bsi_backfill_file,'rb') as my_file:
    Bsi_backfill_data = cPickle.load(my_file)
with open(r'%s'%Bsi_no_backfill_file,'rb') as my_file:
    Bsi_no_backfill_data = cPickle.load(my_file)
	
binBoundaries = np.linspace(-35,-9,100)

plt.subplot(122)

plt.hist(Bsi_backfill_data['DeltaG'][0],facecolor='blue',bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='with DMSO backfill')
plt.hist(Bsi_no_backfill_data['DeltaG'][0],facecolor='green',bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='without DMSO backfill')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

#plt.title('Src affinities',fontsize=20)
plt.yticks([])
plt.ylim((0,0.9))
plt.ylabel('$P(\Delta G)$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((-25,-8.5))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)

plt.tight_layout()
plt.savefig('DMSO-backfill_combined.png', dpi=300)

interval_backfill = np.percentile(a=Bsi_backfill_data['DeltaG'][0], q=[2.5, 50.0, 97.5])
print('backfill: %s' %interval_backfill)

interval_no_backfill = np.percentile(a=Bsi_no_backfill_data['DeltaG'][0], q=[2.5, 50.0, 97.5])
print('no backfill: %s' %interval_no_backfill)


