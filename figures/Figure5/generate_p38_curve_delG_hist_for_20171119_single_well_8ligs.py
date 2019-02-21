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

xml_files = ['./infinite_results/2017-11-20 15-48-06_plate_1.xml',
            './infinite_results/2017-11-20 16-17-59_plate_1.xml',
            './infinite_results/2017-11-20 16-41-35_plate_1.xml',
            './infinite_results/2017-11-20 16-59-09_plate_1.xml',
            './infinite_results/2017-11-20 17-19-31_plate_1.xml',
            './infinite_results/2017-11-20 17-37-37_plate_1.xml',
            './infinite_results/2017-11-20 17-58-04_plate_1.xml',
            './infinite_results/2017-11-20 18-16-02_plate_1.xml',
            './infinite_results/2017-11-20 18-35-42_plate_1.xml',
            './infinite_results/2017-11-20 18-54-33_plate_1.xml',
            './infinite_results/2017-11-20 19-12-49_plate_1.xml',
            './infinite_results/2017-11-20 19-31-32_plate_1.xml']


ligand_conc = [ 0.00000000e+00, 8.00000000e-09, 1.74937932e-08, 3.82541000e-08,
                           8.36511642e-08, 1.82922021e-07, 4.00000000e-07, 8.74689659e-07,
                           1.91270500e-06, 4.18255821e-06, 9.14610104e-06, 2.00000000e-05 ]
inputs = {
    'single_well'   :  True,
    'xml_files'     :  xml_files,
    'file_set'      :  {'dialyzed_p38_1': xml_files},
    'protein_wells'  :  {'dialyzed_p38_1': ['A1','B1','C1','D1','E1','F1','G1','H1']},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Gefitinib','Erlotinib','Ponatinib','Lapatinib','Pazopanib','Axitinib'],
    'buffer_wells'   :  {'dialyzed_p38_1': ['A2','B2','C2','D2','E2','F2','G2','H2']},
    'section'       :  '280_480_TOP_100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
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

cols = sns.color_palette('YlGnBu_r', 12)

# Plot lig 1-4

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Bosutinib-A1A2']/complex_fluorescence['dialyzed_p38_1-Bosutinib-A1A2'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Bosutinib-A1A2']/complex_fluorescence['dialyzed_p38_1-Bosutinib-A1A2'].max(),linestyle='--',color=cols[0])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Bosutinib Isomer-B1B2']/complex_fluorescence['dialyzed_p38_1-Bosutinib Isomer-B1B2'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Bosutinib Isomer-B1B2']/complex_fluorescence['dialyzed_p38_1-Bosutinib Isomer-B1B2'].max(),linestyle='--',color=cols[1])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Erlotinib-D1D2']/complex_fluorescence['dialyzed_p38_1-Erlotinib-D1D2'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Erlotinib-D1D2']/complex_fluorescence['dialyzed_p38_1-Erlotinib-D1D2'].max(),linestyle='--',color=cols[2])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Gefitinib-C1C2']/complex_fluorescence['dialyzed_p38_1-Gefitinib-C1C2'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Gefitinib-C1C2']/complex_fluorescence['dialyzed_p38_1-Gefitinib-C1C2'].max(),linestyle='--',color=cols[3])


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('p38_binding_curve_8lig_1st_half.png', dpi=500)
plt.savefig('p38_binding_curve_8lig_1st_half.pdf')

# Plot ligands 5-8

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Ponatinib-E1E2']/complex_fluorescence['dialyzed_p38_1-Ponatinib-E1E2'].max(),marker='o',color=cols[4],label='Ponatinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Ponatinib-E1E2']/complex_fluorescence['dialyzed_p38_1-Ponatinib-E1E2'].max(),linestyle='--',color=cols[4])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Lapatinib-F1F2']/complex_fluorescence['dialyzed_p38_1-Lapatinib-F1F2'].max(),marker='o',color=cols[5],label='Lapatinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Lapatinib-F1F2']/complex_fluorescence['dialyzed_p38_1-Lapatinib-F1F2'].max(),linestyle='--',color=cols[5])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Pazopanib-G1G2']/complex_fluorescence['dialyzed_p38_1-Pazopanib-G1G2'].max(),marker='o',color=cols[6],label='Pazopanib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Pazopanib-G1G2']/complex_fluorescence['dialyzed_p38_1-Pazopanib-G1G2'].max(),linestyle='--',color=cols[6])

plt.semilogx(inputs['Lstated'],complex_fluorescence['dialyzed_p38_1-Axitinib-H1H2']/complex_fluorescence['dialyzed_p38_1-Axitinib-H1H2'].max(),marker='o',color=cols[7],label='Axitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['dialyzed_p38_1-Axitinib-H1H2']/complex_fluorescence['dialyzed_p38_1-Axitinib-H1H2'].max(),linestyle='--',color=cols[7])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('p38_binding_curve_8lig_2nd_half.png', dpi=500)
plt.savefig('p38_binding_curve_8lig_2nd_half.pdf')


# pickle files of quickmodel analysis with P_error 0.35
# /Users/isikm/lab/wetlab_related/EXPERIMENTS/kinase_inhibitors_and_assaytools/20171119_single_well_binding_assay/from_lilac/20190213_dial_quickmodel_nsampl_100k/dial_p38_col_1_2_quickmodel_singlewavelength_nsample_100k

# pickle files of quickmodel analysis with P_error 0.15
# /Users/isikm/lab/wetlab_related/EXPERIMENTS/kinase_inhibitors_and_assaytools/20171119_single_well_binding_assay/from_lilac/20190213_dial_quickmodel_nsampl_100k_P_error_0_15/dial_p38_col_1_2_quickmodel_singlewavelength_nsample_100k_P_error_0_15

# Now let's plot p38 deltaG histograms for these. Note that these data are from dP_stated = 0.15 analyses.

import _pickle as cPickle

p38_Bos_file = './pickle_Perror_0_15/dialyzed_p38_1-Bosutinib-A1A2_mcmc-2019-02-13 13:03.pickle'
p38_Bsi_file = './pickle_Perror_0_15/dialyzed_p38_1-Bosutinib Isomer-B1B2_mcmc-2019-02-13 15:27.pickle'
p38_Erl_file = './pickle_Perror_0_15/dialyzed_p38_1-Erlotinib-D1D2_mcmc-2019-02-13 20:45.pickle'
p38_Gef_file = './pickle_Perror_0_15/dialyzed_p38_1-Gefitinib-C1C2_mcmc-2019-02-13 18:02.pickle'

p38_Pon_file = './pickle_Perror_0_15/dialyzed_p38_1-Ponatinib-E1E2_mcmc-2019-02-13 23:09.pickle'
p38_Lap_file = './pickle_Perror_0_15/dialyzed_p38_1-Lapatinib-F1F2_mcmc-2019-02-14 01:33.pickle'
p38_Paz_file = './pickle_Perror_0_15/dialyzed_p38_1-Pazopanib-G1G2_mcmc-2019-02-14 03:52.pickle'
p38_Axi_file = './pickle_Perror_0_15/dialyzed_p38_1-Axitinib-H1H2_mcmc-2019-02-14 06:11.pickle'

with open(r'%s'%p38_Bos_file,'rb') as my_file:
    p38_Bos_data = cPickle.load(my_file)
with open(r'%s'%p38_Bsi_file,'rb') as my_file:
    p38_Bsi_data = cPickle.load(my_file)
with open(r'%s'%p38_Erl_file,'rb') as my_file:
    p38_Erl_data = cPickle.load(my_file)
with open(r'%s'%p38_Gef_file,'rb') as my_file:
    p38_Gef_data = cPickle.load(my_file)

with open(r'%s'%p38_Pon_file,'rb') as my_file:
    p38_Pon_data = cPickle.load(my_file)
with open(r'%s'%p38_Lap_file,'rb') as my_file:
    p38_Lap_data = cPickle.load(my_file)
with open(r'%s'%p38_Paz_file,'rb') as my_file:
    p38_Paz_data = cPickle.load(my_file)
with open(r'%s'%p38_Axi_file,'rb') as my_file:
    p38_Axi_data = cPickle.load(my_file)


p38_bosutinib = 3000e-9 # >3000 nM from DiscoverRx screen data
p38_erlotinib = 3000e-9 # >3000 nM from DiscoverRx screen data
p38_gefitinib = 3000e-9 # >3000 nM from DiscoverRx screen data

p38Bos_dG = np.log(p38_bosutinib)
p38Erl_dG = np.log(p38_erlotinib)
p38Gef_dG = np.log(p38_gefitinib)

#cols = sns.color_palette('PuBu_r', 5) # Color of p38 plots in Figure 3
cols = sns.color_palette('YlGnBu_r', 12) # Color of Src plots in Figure 3


binBoundaries = np.linspace(-35,-9,50)

kd_binBoundaries = np.exp(np.arange(-30,-9,0.5))

#Lets make this plot both using or thermal units for delG and molar units for Kd

fig, ax = plt.subplots(figsize=(8,4))

plt.hist(p38_Bos_data['DeltaG'][0],facecolor=cols[0],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Bosutinib')
plt.hist(p38_Bsi_data['DeltaG'][0],facecolor=cols[1],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Bosutinib Isomer')
plt.hist(p38_Erl_data['DeltaG'][0],facecolor=cols[2],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Erlotinib')
plt.hist(p38_Gef_data['DeltaG'][0],facecolor=cols[3],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Gefitinib')

plt.hist(p38_Pon_data['DeltaG'][0],facecolor=cols[4],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Ponatinib')
plt.hist(p38_Lap_data['DeltaG'][0],facecolor=cols[5],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Lapatinib')
plt.hist(p38_Paz_data['DeltaG'][0],facecolor=cols[6],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Pazopanib')
plt.hist(p38_Axi_data['DeltaG'][0],facecolor=cols[7],bins=binBoundaries,edgecolor='white',normed=1,alpha=0.9,label='p38:Axitinib')


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

plt.hist(np.exp(p38_Pon_data['DeltaG'][0]),facecolor=cols[4],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Ponabitib')
plt.hist(np.exp(p38_Lap_data['DeltaG'][0]),facecolor=cols[5],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Lapatinib')
plt.hist(np.exp(p38_Paz_data['DeltaG'][0]),facecolor=cols[6],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Pazopanib')
plt.hist(np.exp(p38_Axi_data['DeltaG'][0]),facecolor=cols[7],bins=kd_binBoundaries,edgecolor='white',alpha=0.9,label='p38:Axitinib')


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

plt.savefig('p38_Kd_hist_8lig.png', dpi=500)
plt.savefig('p38_Kd_hist_8lig.pdf')