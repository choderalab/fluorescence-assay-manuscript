import os
import numpy as np
from assaytools import parser
import string
from glob import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

xml_files = ['./infinite_results/p38_Abl_WT_GK_Src_WT_GK_conc_0_20190307_111234.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GK_conc_1_20190307_112351.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_2_20190307_113345.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_3_20190307_114336.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_4_20190307_115329.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_5_20190307_120322.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_6_20190307_121315.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_7_20190307_122657.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_8_20190307_123649.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_9_20190307_124642.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_10_20190307_125635.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_11_20190307_130625.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_12_20190307_131618.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_13_20190307_132611.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_14_20190307_133612.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_15_20190307_134607.xml',
             './infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_16_20190307_135607.xml']

ligand_conc = [  0.00000000e+00,   8.00000000e-09,   1.34778097e-08,
         2.27064194e-08,   3.82541000e-08,   6.44476851e-08,
         1.08576705e-07,   1.82922021e-07,   3.08173524e-07,
         5.19188015e-07,   8.74689659e-07,   1.47361260e-06, 2.48263378e-06,
         4.18255821e-06, 7.04646547e-06, 1.118713651e-05, 2.0e-05]

inputs = {
    'single_well'   :  True,
    'xml_files'     :  xml_files,
    'file_set'      :  {'p38': xml_files},
    'protein_wells'  :  {'p38': ['A11', 'C11', 'E11', 'G11']},
    'ligand_order'  :  ['Bosutinib', 'Bosutinib Isomer', 'Erlotinib', 'Gefitinib'],
    'buffer_wells'   :  {'p38': ['A2', 'C2', 'E2', 'G2']},
    'section'       :  'ex280_em480_top_gain100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  0.5e-6 * np.ones([17],np.float64), # protein concentration, M
    'P_error'       :  0.15, # coefficient of protein concentration uncertainity (0.15 for 15% error)
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }
 
[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)

print("complex_fluorescence:\n", complex_fluorescence)
print("ligand_fluorescence:\n", ligand_fluorescence)


# Test matplotlib

#First create some toy data:
#x = np.linspace(0, 2*np.pi, 400)
#y = np.sin(x**2)

#Creates just a figure and only one subplot
#fig, ax = plt.subplots()
#ax.plot(x, y)
#ax.set_title('Simple plot')

#plt.savefig('test.pdf')

# Plot fluorescence binding curve
cols = sns.color_palette('YlGnBu_r', 5)

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib-A11A2']/complex_fluorescence['p38-Bosutinib-A11A2'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib-A11A2']/complex_fluorescence['p38-Bosutinib-A11A2'].max(),linestyle='--',color=cols[0])

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib Isomer-C11C2']/complex_fluorescence['p38-Bosutinib Isomer-C11C2'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib Isomer-C11C2']/complex_fluorescence['p38-Bosutinib Isomer-C11C2'].max(),linestyle='--',color=cols[1])

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Erlotinib-E11E2']/complex_fluorescence['p38-Erlotinib-E11E2'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Erlotinib-E11E2']/complex_fluorescence['p38-Erlotinib-E11E2'].max(),linestyle='--',color=cols[2])

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Gefitinib-G11G2']/complex_fluorescence['p38-Gefitinib-G11G2'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Gefitinib-G11G2']/complex_fluorescence['p38-Gefitinib-G11G2'].max(),linestyle='--',color=cols[3])


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38 in Corning 3651 plate',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('p38_binding_curve_Corning_3651.png', dpi=500)
plt.savefig('p38_binding_curve_Corning_3651.pdf')



# pickle files of quickmodel analysis with P_error 0.35
# /Users/isikm/lab/wetlab_related/EXPERIMENTS/kinase_inhibitors_and_assaytools/20171119_single_well_binding_assay/from_lilac/20190213_dial_quickmodel_nsampl_100k/dial_p38_col_1_2_quickmodel_singlewavelength_nsample_100k

# pickle files of quickmodel analysis with P_error 0.15
# /Users/isikm/lab/wetlab_related/EXPERIMENTS/kinase_inhibitors_and_assaytools/20171119_single_well_binding_assay/from_lilac/20190213_dial_quickmodel_nsampl_100k_P_error_0_15/dial_p38_col_1_2_quickmodel_singlewavelength_nsample_100k_P_error_0_15

# Now let's plot p38 deltaG histograms for these. Note that these data are from dP_stated = 0.15 analyses.

import _pickle as cPickle

p38_Bos_file = './pickle_Perror_0_15/p38-Bosutinib-A11A2_mcmc-2019-03-07 16:05.pickle'
p38_Bsi_file = './pickle_Perror_0_15/p38-Bosutinib Isomer-C11C2_mcmc-2019-03-07 18:22.pickle'
p38_Erl_file = './pickle_Perror_0_15/p38-Erlotinib-E11E2_mcmc-2019-03-07 20:34.pickle'
p38_Gef_file = './pickle_Perror_0_15/p38-Gefitinib-G11G2_mcmc-2019-03-07 22:45.pickle'

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

#cols = sns.color_palette('PuBu_r', 5) # Color of p38 plots in Figure 3
cols = sns.color_palette('YlGnBu_r', 5)


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

plt.title('p38 affinities in Corning 3651 plate',fontsize=20)
plt.yticks([])
plt.ylim((0,0.9))
plt.ylabel('$P(\Delta G)$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((-25,-8.5))
plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)

plt.tight_layout()

plt.savefig('p38_delG_hist_Corning_3651.png', dpi=500)
plt.savefig('p38_delG_hist_Corning_3651.pdf')

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

plt.title('p38 affinities in Corning 3651 plate',fontsize=20)
plt.yticks([])
#plt.ylim((0,0.9))
plt.ylabel('$P(K_{d})$',fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((1e-11,2e-4))
plt.xlabel('$K_{d}$ ($M$)',fontsize=16)
plt.legend(loc=2,fontsize=14,frameon=True,framealpha=0.9)
plt.xscale('log')

plt.tight_layout()

plt.savefig('p38_Kd_hist_Corning_3651.png', dpi=500)
plt.savefig('p38_Kd_hist_Corning_3651.pdf')
