from glob import glob
from assaytools import parser
import numpy as np

import seaborn as sns
cols = sns.color_palette("deep", 3)
sns.set(style='white')
sns.set_context('talk')

import matplotlib.pyplot as plt

# We are making this figure using Src:Bosutinib Isomer data

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
    
fig, ax = plt.subplots(figsize=(5,3))

plt.semilogx(inputs['Lstated'],complex_fluorescence['Src-Bosutinib Isomer-CD'],'ko',label='protein+ligand')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['Src-Bosutinib Isomer-CD'],'ro',label='buffer+ligand')

plt.xticks(fontsize=15)
plt.yticks([])
plt.xlabel('Ligand Concentration [M]', fontsize=18);
plt.ylabel('Fluorescence', fontsize=18);
plt.legend(fontsize=15);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('example_data.png', dpi=500)
plt.savefig('example_data.pdf')

import cPickle

data_file = './Src-Bosutinib Isomer-CD_mcmc-1.pickle'
with open(r'%s'%data_file,'rb') as my_file:
    data = cPickle.load(my_file)

#First we want to plot our Fluorescence data with the fits to it

fig, ax = plt.subplots(figsize=(5,3))

for top_complex_fluorescence_model in data['top_complex_fluorescence_model'][0][::50]:
    plt.semilogx(inputs['Lstated'], top_complex_fluorescence_model, marker='.',color='silver')
for top_ligand_fluorescence_model in data['top_ligand_fluorescence_model'][0][::50]:
    plt.semilogx(inputs['Lstated'], top_ligand_fluorescence_model, marker='.',color='lightcoral', alpha=0.2)
plt.semilogx(inputs['Lstated'], complex_fluorescence['Src-Bosutinib Isomer-CD'], 'ko',label='complex')
plt.semilogx(inputs['Lstated'], ligand_fluorescence['Src-Bosutinib Isomer-CD'], marker='o',color='firebrick',linestyle='None',label='ligand')

plt.xticks(fontsize=15)
plt.yticks([])
plt.xlabel('Ligand Concentration [M]', fontsize=18);
plt.ylabel('Fluorescence', fontsize=18);
plt.legend(fontsize=15);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('fit_to_data.png', dpi=500)
plt.savefig('fit_to_data.pdf')

#Next we want to plot our DeltaG trace with the histogram and 95% credibility interval

import matplotlib.patches as mpatches
import matplotlib.lines as mlines

interval = np.percentile(a=data['DeltaG'][0], q=[2.5, 50.0, 97.5])
[hist,bin_edges] = np.histogram(data['DeltaG'][0],bins=40,normed=True)
binwidth = np.abs(bin_edges[0]-bin_edges[1])

#set colors for 95% interval
clrs = [(0.7372549019607844, 0.5098039215686274, 0.7411764705882353) for xx in bin_edges]
idxs = bin_edges.argsort()
idxs = idxs[::-1]
gray_before = idxs[bin_edges[idxs] < interval[0]]
gray_after = idxs[bin_edges[idxs] > interval[2]]
for idx in gray_before:
    clrs[idx] = (.5,.5,.5)
for idx in gray_after:
    clrs[idx] = (.5,.5,.5)

hist_legend = mpatches.Patch(color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
    label = '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$'
    %(interval[1],interval[0],interval[2]) )
    
f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))

ax1.plot(data['DeltaG'][0],color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353))
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,99000)

ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)
ax2.barh(bin_edges[:-1],hist,binwidth,color=clrs, edgecolor = "white");
ax2.axhline(y=interval[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.legend(handles=[hist_legend],fontsize=16,loc=0,frameon=True)
ax2.set_xlabel('$P(\Delta G)$',fontsize=16);
plt.xticks([])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('DeltaG_trace_hist.png', dpi=500)
plt.savefig('DeltaG_trace_hist..pdf')
    
#Finally we want to plot our F_PL and F_L traces

fig, ax = plt.subplots(figsize=(9,3))

plt.plot(data['F_PL'][0])

plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('MCMC Sample', fontsize=18);
plt.ylabel('$F\_PL$', fontsize=20);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('F_PL.pdf')
plt.savefig('F_PL.png',dpi=500)

fig, ax = plt.subplots(figsize=(9,3))

plt.plot(data['F_L'][0],color='green')

plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('MCMC Sample', fontsize=18);
plt.ylabel('$F\_L$', fontsize=20);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('F_L.pdf')
plt.savefig('F_L.png',dpi=500)
