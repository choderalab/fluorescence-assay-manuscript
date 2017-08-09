
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
    'file_set'      :  {'Src': glob("../../data/spectra/Src/2015-12-15/*.xml"),
                        'SrcGK': glob("../../data/spectra/SrcT338I/*.xml"),
                        'AblGK': glob("../../data/spectra/AblT334I/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml")},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)

cols = sns.color_palette('YlGnBu_r', 5)

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['Src-Bosutinib-AB']/complex_fluorescence['Src-Bosutinib-AB'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['Src-Bosutinib-AB']/complex_fluorescence['Src-Bosutinib-AB'].max(),linestyle='--',color=cols[0])
plt.semilogx(inputs['Lstated'],complex_fluorescence['Src-Bosutinib Isomer-CD']/complex_fluorescence['Src-Bosutinib Isomer-CD'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['Src-Bosutinib Isomer-CD']/complex_fluorescence['Src-Bosutinib Isomer-CD'].max(),linestyle='--',color=cols[1])
plt.semilogx(inputs['Lstated'],complex_fluorescence['Src-Erlotinib-EF']/complex_fluorescence['Src-Erlotinib-EF'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['Src-Erlotinib-EF']/complex_fluorescence['Src-Erlotinib-EF'].max(),linestyle='--',color=cols[2])
plt.semilogx(inputs['Lstated'],complex_fluorescence['Src-Gefitinib-GH']/complex_fluorescence['Src-Gefitinib-GH'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['Src-Gefitinib-GH']/complex_fluorescence['Src-Gefitinib-GH'].max(),linestyle='--',color=cols[3])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('Src',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('Src_binding_curve.png', dpi=500)
plt.savefig('Src_binding_curve.pdf')

#Now let's plot Src deltaG histograms for these. Note that these data are from dP_stated = 0.15 analyses.

import cPickle

Src_Bos_file = 'dPstated/Src-Bosutinib-AB_mcmc-0.pickle'
Src_Bsi_file = 'dPstated/Src-Bosutinib Isomer-CD_mcmc-0.pickle'
Src_Erl_file = 'dPstated/Src-Erlotinib-EF_mcmc-0.pickle'
Src_Gef_file = 'dPstated/Src-Gefitinib-GH_mcmc-0.pickle'

with open(r'%s'%Src_Bos_file,'rb') as my_file:
    Src_Bos_data = cPickle.load(my_file)
with open(r'%s'%Src_Bsi_file,'rb') as my_file:
    Src_Bsi_data = cPickle.load(my_file)
with open(r'%s'%Src_Erl_file,'rb') as my_file:
    Src_Erl_data = cPickle.load(my_file)
with open(r'%s'%Src_Gef_file,'rb') as my_file:
    Src_Gef_data = cPickle.load(my_file)

Src_bosutinib = 1e-9 # 1 nM from DiscoverRx screen data 
Src_erlotinib = 700e-9 # 700 nM from DiscoverRx screen data 
Src_gefitinib = 3800e-9 # 3800 nM from DiscoverRx screen data 

SrcBos_dG = np.log(Src_bosutinib)
SrcErl_dG = np.log(Src_erlotinib)
SrcGef_dG = np.log(Src_gefitinib)

cols = sns.color_palette('YlGnBu_r', 5)

binBoundaries = np.linspace(-35,-9,50)

kd_binBoundaries = np.exp(np.arange(-30,-9,0.5))

#Lets make this plot both using or thermal units for delG and molar units for Kd

fig, ax = plt.subplots(figsize=(8,4))

plt.hist(Src_Bos_data['DeltaG'][0],facecolor=cols[0],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='Src:Bosutinib')
plt.hist(Src_Bsi_data['DeltaG'][0],facecolor=cols[1],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='Src:Bosutinib Isomer')
plt.hist(Src_Erl_data['DeltaG'][0],facecolor=cols[2],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='Src:Erlotinib')
plt.hist(Src_Gef_data['DeltaG'][0],facecolor=cols[3],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='Src:Gefitinib')

plt.axvline(x=SrcBos_dG,color=cols[0],linestyle='--')
plt.axvline(x=SrcErl_dG,color=cols[2],linestyle='--')
plt.axvline(x=SrcGef_dG,color=cols[3],linestyle='--')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('Src affinities',fontsize=20)
plt.yticks([])
plt.ylim((0,0.9))
plt.ylabel('$P(\Delta G)$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((-25,-8.5))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)

plt.tight_layout()

plt.savefig('Src_delG_hist.png', dpi=500)
plt.savefig('Src_delG_hist.pdf')

fig, ax = plt.subplots(figsize=(8,4))

plt.hist(np.exp(Src_Bos_data['DeltaG'][0]),facecolor=cols[0],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='Src:Bosutinib')
plt.hist(np.exp(Src_Bsi_data['DeltaG'][0]),facecolor=cols[1],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='Src:Bosutinib Isomer')
plt.hist(np.exp(Src_Erl_data['DeltaG'][0]),facecolor=cols[2],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='Src:Erlotinib')
plt.hist(np.exp(Src_Gef_data['DeltaG'][0]),facecolor=cols[3],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='Src:Gefitinib')

plt.axvline(x=Src_bosutinib,color=cols[0],linestyle='--')
plt.axvline(x=Src_erlotinib,color=cols[2],linestyle='--')
plt.axvline(x=Src_gefitinib,color=cols[3],linestyle='--')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('Src affinities',fontsize=20)
plt.yticks([])
#plt.ylim((0,0.9))
plt.ylabel('$P(K_{d})$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((1e-11,2e-4))
plt.xlabel('$K_{d}$ ($M$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)
plt.xscale('log')

plt.tight_layout()

plt.savefig('Src_Kd_hist.png', dpi=500)
plt.savefig('Src_Kd_hist.pdf')










