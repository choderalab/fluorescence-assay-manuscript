# We are interested in making a plot to show that decreasing protein concentration in this format leads to worse results

# First I will use the inputs.py format to import and parse the xml and overlay results for different protein concentrations

import numpy as np
from glob import glob
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

inputs_p5 = {
    'xml_file_path' :  "../../data/single-wavelength/",
    'file_set'      :  {'p38': glob( "../../data/single-wavelength/p38_0.5*.xml")},
    'section'       :  '280_480_TOP_120',
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Ponatinib','Lapatinib','Saracatinib','Vandetanib'],
    'Lstated'       :  np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Pstated'       :  0.5e-6 * np.ones([24],np.float64), # protein concentration, M
    'assay_volume'  :  50e-6, # assay volume, L
    'well_area'     :  0.1369, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }
    
inputs_p25 = {
    'xml_file_path' :  "../../data/single-wavelength/",
    'file_set'      :  {'p38': glob( "../../data/single-wavelength/p38_0.25*.xml")},
    'section'       :  '280_480_TOP_120',
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Ponatinib','Lapatinib','Saracatinib','Vandetanib'],
    'Lstated'       :  np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Pstated'       :  0.25e-6 * np.ones([24],np.float64), # protein concentration, M
    'assay_volume'  :  50e-6, # assay volume, L
    'well_area'     :  0.1369, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }
    
inputs_p125 = {
    'xml_file_path' :  "../../data/single-wavelength/",
    'file_set'      :  {'p38': glob( "../../data/single-wavelength/p38_0.125*.xml")},
    'section'       :  '280_480_TOP_120',
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Ponatinib','Lapatinib','Saracatinib','Vandetanib'],
    'Lstated'       :  np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Pstated'       :  0.125e-6 * np.ones([24],np.float64), # protein concentration, M
    'assay_volume'  :  50e-6, # assay volume, L
    'well_area'     :  0.1369, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }
	
from assaytools import parser

[complex_fluorescence_p5, ligand_fluorescence_p5] = parser.get_data_using_inputs(inputs_p5)  
[complex_fluorescence_p25, ligand_fluorescence_p25] = parser.get_data_using_inputs(inputs_p25)  
[complex_fluorescence_p125, ligand_fluorescence_p125] = parser.get_data_using_inputs(inputs_p125) 

choice = 'p38-Bosutinib-AB'

fig, ax = plt.subplots(figsize=(16,4))

plt.subplot(121)
plt.semilogx(inputs_p5['Lstated'], complex_fluorescence_p5[choice], 'bo',label='0.5 uM (complex)')
plt.semilogx(inputs_p5['Lstated'], ligand_fluorescence_p5[choice], marker='.',color='cornflowerblue',linestyle='None',label='0.5 uM (ligand)')
plt.semilogx(inputs_p25['Lstated'], complex_fluorescence_p25[choice], 'go',linestyle='None',label='0.25 uM (complex)')
plt.semilogx(inputs_p25['Lstated'], ligand_fluorescence_p25[choice], marker='.',color='mediumseagreen',linestyle='None',label='0.25 uM (ligand)')
plt.semilogx(inputs_p125['Lstated'], complex_fluorescence_p125[choice], marker='o',color='cyan',linestyle='None',label='0.125 uM (complex)')
plt.semilogx(inputs_p125['Lstated'], ligand_fluorescence_p125[choice], marker='.',color='paleturquoise',linestyle='None',label='0.125 uM (ligand)')
plt.xlabel('$[L]_T$ (M)');
plt.ylabel('fluorescence units');
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(loc=0);


# Then I will import the pickle file for the analysis that was already run to plot the histograms

import cPickle

p5_file = '../../analysis/bayes/8ligs_0.5_protein_concentration/Bosutinib-AB_mcmc-2017-03-23 22:14.pickle'
p25_file = '../../analysis/bayes/8ligs_0.25_protein_concentration/Bosutinib-AB_mcmc-2017-03-23 23:42.pickle'
p125_file = '../../analysis/bayes/8ligs_0.125_protein_concentration/Bosutinib-AB_mcmc-2017-03-24 01:05.pickle'


with open(r'%s'%p5_file,'rb') as my_file:
    p5_data = cPickle.load(my_file)
with open(r'%s'%p25_file,'rb') as my_file:
    p25_data = cPickle.load(my_file)
with open(r'%s'%p125_file,'rb') as my_file:
    p125_data = cPickle.load(my_file)
	
binBoundaries = np.linspace(-35,-9,100)

plt.subplot(122)

plt.hist(p5_data['DeltaG'][0],facecolor='blue',bins=binBoundaries,edgecolor='white',normed=1,alpha=0.8,label='0.5 uM')
plt.hist(p25_data['DeltaG'][0],facecolor='green',bins=binBoundaries,edgecolor='white',normed=1,alpha=0.8,label='0.25 uM')
plt.hist(p125_data['DeltaG'][0],facecolor='cyan',bins=binBoundaries,edgecolor='white',normed=1,alpha=0.6,label='0.125 uM')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.yticks([])
plt.ylim((0,0.4))
plt.ylabel('$P(\Delta G)$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((-25,-8.5))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)

plt.tight_layout()
plt.savefig('protein-concentration_combined.png', dpi=300)



