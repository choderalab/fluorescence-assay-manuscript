
import numpy as np
from assaytools import parser
import string
from glob import glob

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

inputs = {
    'xml_file_path' :  "../../data/spectra/",
    'file_set'      :  {'p38': glob("../../data/spectra/p38/2016-01-26/*.xml")},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)

cols = sns.color_palette('PuBu_r', 5)

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib-AB']/complex_fluorescence['p38-Bosutinib-AB'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib-AB']/complex_fluorescence['p38-Bosutinib-AB'].max(),linestyle='--',color=cols[0])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib Isomer-CD']/complex_fluorescence['p38-Bosutinib Isomer-CD'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib Isomer-CD']/complex_fluorescence['p38-Bosutinib Isomer-CD'].max(),linestyle='--',color=cols[1])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Erlotinib-EF']/complex_fluorescence['p38-Erlotinib-EF'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Erlotinib-EF']/complex_fluorescence['p38-Erlotinib-EF'].max(),linestyle='--',color=cols[2])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Gefitinib-GH']/complex_fluorescence['p38-Gefitinib-GH'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Gefitinib-GH']/complex_fluorescence['p38-Gefitinib-GH'].max(),linestyle='--',color=cols[3])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('p38_binding_curve.png', dpi=500)
plt.savefig('p38_binding_curve.pdf')

#Now let's plot p38 deltaG histograms for these. Note that these data are from dP_stated = 0.15 analyses.

import cPickle

p38_Bos_file = 'dPstated/p38-Bosutinib-AB_mcmc-0.pickle'
p38_Bsi_file = 'dPstated/p38-Bosutinib Isomer-CD_mcmc-0.pickle'
p38_Erl_file = 'dPstated/p38-Erlotinib-EF_mcmc-0.pickle'
p38_Gef_file = 'dPstated/p38-Gefitinib-GH_mcmc-0.pickle'

with open(r'%s'%p38_Bos_file,'rb') as my_file:
    p38_Bos_data = cPickle.load(my_file)
with open(r'%s'%p38_Bsi_file,'rb') as my_file:
    p38_Bsi_data = cPickle.load(my_file)
with open(r'%s'%p38_Erl_file,'rb') as my_file:
    p38_Erl_data = cPickle.load(my_file)
with open(r'%s'%p38_Gef_file,'rb') as my_file:
    p38_Gef_data = cPickle.load(my_file)

p38_bosutinib = 3000e-9 # >3000 nM from DiscoverRx screen data 
p38_erlotinib = 3000e-9 # >3000 nM from DiscoverRx screen data 
p38_gefitinib = 3000e-9 # >3000 nM from DiscoverRx screen data 

p38Bos_dG = np.log(p38_bosutinib)
p38Erl_dG = np.log(p38_erlotinib)
p38Gef_dG = np.log(p38_gefitinib)

cols = sns.color_palette('PuBu_r', 5)

binBoundaries = np.linspace(-35,-9,50)

kd_binBoundaries = np.exp(np.arange(-30,-9,0.5))

#Lets make this plot both using or thermal units for delG and molar units for Kd

fig, ax = plt.subplots(figsize=(8,4))

plt.hist(p38_Bos_data['DeltaG'][0],facecolor=cols[0],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Bosutinib')
plt.hist(p38_Bsi_data['DeltaG'][0],facecolor=cols[1],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Bosutinib Isomer')
plt.hist(p38_Erl_data['DeltaG'][0],facecolor=cols[2],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Erlotinib')
plt.hist(p38_Gef_data['DeltaG'][0],facecolor=cols[3],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Gefitinib')

plt.axvline(x=p38Bos_dG,color=cols[0],linestyle='--')
plt.axvline(x=p38Erl_dG,color=cols[2],linestyle='--')
plt.axvline(x=p38Gef_dG,color=cols[3],linestyle='--')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38 affinities',fontsize=20)
plt.yticks([])
plt.ylim((0,0.9))
plt.ylabel('$P(\Delta G)$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((-25,-8.5))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)

plt.tight_layout()

plt.savefig('p38_delG_hist.png', dpi=500)
plt.savefig('p38_delG_hist.pdf')

fig, ax = plt.subplots(figsize=(8,4))

plt.hist(np.exp(p38_Bos_data['DeltaG'][0]),facecolor=cols[0],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Bosutinib')
plt.hist(np.exp(p38_Bsi_data['DeltaG'][0]),facecolor=cols[1],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Bosutinib Isomer')
plt.hist(np.exp(p38_Erl_data['DeltaG'][0]),facecolor=cols[2],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Erlotinib')
plt.hist(np.exp(p38_Gef_data['DeltaG'][0]),facecolor=cols[3],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Gefitinib')

plt.axvline(x=p38_bosutinib,color=cols[0],linestyle='--')
plt.axvline(x=p38_erlotinib,color=cols[2],linestyle='--')
plt.axvline(x=p38_gefitinib,color=cols[3],linestyle='--')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38 affinities',fontsize=20)
plt.yticks([])
#plt.ylim((0,0.9))
plt.ylabel('$P(K_{d})$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((1e-11,2e-4))
plt.xlabel('$K_{d}$ ($M$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)
plt.xscale('log')

plt.tight_layout()

plt.savefig('p38_Kd_hist.png', dpi=500)
plt.savefig('p38_Kd_hist.pdf')










