
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
                        'p38': glob("../../data/spectra/p38/2016-01-26/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml")},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       : [ 
        np.array([2.00736e-05,9.179758690368378e-06,4.197950024579236e-06,1.9197437539785283e-06,8.779085171003113e-07,4.014719999999999e-07,1.835951738073675e-07,8.395900049158471e-08,3.839487507957056e-08,1.7558170342006224e-08,8.029439999999995e-09, 0.0], np.float64), 
        np.array([1.9999780000000003e-05,9.146000431435104e-06,4.182512202224779e-06,1.9126839598250786e-06,8.746800375683717e-07,3.9999559999999996e-07,1.8292000862870203e-07,8.365024404449556e-08,3.825367919650156e-08,1.749360075136743e-08,7.999911999999997e-09, 0.0], np.float64),
        np.array([1.838576e-05,8.407900931523358e-06,3.844975572090105e-06,1.7583267536539667e-06,8.040917073849344e-07,3.677151999999999e-07,1.6815801863046714e-07,7.689951144180208e-08,3.516653507307932e-08,1.6081834147698685e-08,7.354303999999996e-09, 0.0], np.float64), 
        np.array([2.0035800000000002e-05,9.162472559405526e-06,4.1900449895616465e-06,1.916128741529322e-06,8.762553536445092e-07,4.0071599999999995e-07,1.8324945118811045e-07,8.38008997912329e-08,3.8322574830586434e-08,1.752510707289018e-08,8.014319999999996e-09, 0.0], np.float64)
         ],# corrected ligand concentrations
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)

cols = sns.color_palette('YlGnBu_r', 5)

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'][0],complex_fluorescence['Src-Bosutinib-AB']/complex_fluorescence['Src-Bosutinib-AB'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'][0],ligand_fluorescence['Src-Bosutinib-AB']/complex_fluorescence['Src-Bosutinib-AB'].max(),linestyle='--',color=cols[0])
plt.semilogx(inputs['Lstated'][1],complex_fluorescence['Src-Bosutinib Isomer-CD']/complex_fluorescence['Src-Bosutinib Isomer-CD'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'][1],ligand_fluorescence['Src-Bosutinib Isomer-CD']/complex_fluorescence['Src-Bosutinib Isomer-CD'].max(),linestyle='--',color=cols[1])
plt.semilogx(inputs['Lstated'][2],complex_fluorescence['Src-Erlotinib-EF']/complex_fluorescence['Src-Erlotinib-EF'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'][2],ligand_fluorescence['Src-Erlotinib-EF']/complex_fluorescence['Src-Erlotinib-EF'].max(),linestyle='--',color=cols[2])
plt.semilogx(inputs['Lstated'][3],complex_fluorescence['Src-Gefitinib-GH']/complex_fluorescence['Src-Gefitinib-GH'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'][3],ligand_fluorescence['Src-Gefitinib-GH']/complex_fluorescence['Src-Gefitinib-GH'].max(),linestyle='--',color=cols[3])

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

import pickle

Src_Bos_file = '../../analysis/bayes/corrected-concentrations/bos/Src-Bosutinib-AB_mcmc-2019-05-17 13:44.pickle'
Src_Bsi_file = '../../analysis/bayes/corrected-concentrations/bos_iso/Src-Bosutinib Isomer-CD_mcmc-2019-05-17 13:45.pickle'
Src_Erl_file = '../../analysis/bayes/corrected-concentrations/erl/Src-Erlotinib-EF_mcmc-2019-05-17 13:48.pickle'
Src_Gef_file = '../../analysis/bayes/corrected-concentrations/gef/Src-Gefitinib-GH_mcmc-2019-05-20 10:30.pickle'

with open(r'%s'%Src_Bos_file,'rb') as my_file:
    Src_Bos_data = pickle.load(my_file)
with open(r'%s'%Src_Bsi_file,'rb') as my_file:
    Src_Bsi_data = pickle.load(my_file)
with open(r'%s'%Src_Erl_file,'rb') as my_file:
    Src_Erl_data = pickle.load(my_file)
with open(r'%s'%Src_Gef_file,'rb') as my_file:
    Src_Gef_data = pickle.load(my_file)

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










