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

for top_complex_fluorescence_model in data['top_complex_fluorescence_model'][0][::100]:
    plt.semilogx(inputs['Lstated'], top_complex_fluorescence_model, marker='.',color='silver')
for top_ligand_fluorescence_model in data['top_ligand_fluorescence_model'][0][::100]:
    plt.semilogx(inputs['Lstated'], top_ligand_fluorescence_model, marker='.',color='lightcoral', alpha=0.2)
plt.semilogx(inputs['Lstated'], complex_fluorescence['Src-Bosutinib Isomer-CD'], 'ko',label='protein+ligand')
plt.semilogx(inputs['Lstated'], data['top_complex_fluorescence_model'][0][0], marker='.',color='silver',label='model fit')
plt.semilogx(inputs['Lstated'], ligand_fluorescence['Src-Bosutinib Isomer-CD'], marker='o',color='firebrick',linestyle='None',label='buffer+ligand')
plt.semilogx(inputs['Lstated'], data['top_ligand_fluorescence_model'][0][0], marker='.',color='lightcoral', alpha=0.2, label='model fit')
#These two lines are a repeat of the above to get the data points above the fit points, but still keep the legend in the right order.
plt.semilogx(inputs['Lstated'], complex_fluorescence['Src-Bosutinib Isomer-CD'], 'ko')
plt.semilogx(inputs['Lstated'], ligand_fluorescence['Src-Bosutinib Isomer-CD'], marker='o',color='firebrick',linestyle='None')

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
from matplotlib import gridspec
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

ax1 = plt.subplot(gs[0])
ax1.plot(range(0,len(data['DeltaG'][0]),10),data['DeltaG'][0][::10],color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353))
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,99000)

ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)

ax2 = plt.subplot(gs[1])
ax2.barh(bin_edges[:-1],hist,binwidth,color=clrs, edgecolor = "white");
ax2.axhline(y=interval[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.set_xlabel('$P(\Delta G)$',fontsize=16);
plt.xticks([])
plt.yticks([])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.savefig('DeltaG_trace_hist.png', dpi=500, bbox_inches='tight')
plt.savefig('DeltaG_trace_hist.pdf', bbox_inches='tight')

#Let's make the same plot, but in units of Kd

interval_Kd = np.percentile(a=np.exp(data['DeltaG'][0]), q=[2.5, 50.0, 97.5])
kd_binBoundaries = np.exp(np.arange(-18,-10,0.15))

hist_legend = mpatches.Patch(color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
    label = '$K_{d}$ = %.1f [%.1f,%.1f] uM'
    %(interval_Kd[1]/1e-6,interval_Kd[0]/1e-6,interval_Kd[2]/1e-6) )

f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])

ax1 = plt.subplot(gs[0])
ax1.plot(range(0,len(data['DeltaG'][0]),10),np.exp(data['DeltaG'][0][::10]),color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353))
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$K_{d}$ ($M$)',fontsize=16);
ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,99000)
ax1.set_ylim((1e-8,1e-5))
plt.yscale('log')
ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)

ax2 = plt.subplot(gs[1])
n, bins, patches = ax2.hist(np.exp(data['DeltaG'][0]),color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353),bins=kd_binBoundaries, edgecolor='white',orientation="horizontal")
ax2.axhline(y=interval_Kd[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_Kd[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_Kd[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.set_xlabel('$P(K_{d})$',fontsize=16);
ax2.set_ylim((1e-8,1e-5))
plt.yscale('log')
plt.xticks([])
plt.yticks([])

#set colors for 95% interval
clrs = [(0.7372549019607844, 0.5098039215686274, 0.7411764705882353) for xx in bins]
idxs = bins.argsort()
idxs = idxs[::-1]
gray_before = idxs[bins[idxs] < interval_Kd[0]]
gray_after = idxs[bins[idxs] > interval_Kd[2]]
for idx in gray_before:
    clrs[idx] = (.5,.5,.5)
for idx in gray_after:
    clrs[idx] = (.5,.5,.5)

for i,p in enumerate(patches):
    plt.setp(p, 'facecolor', clrs[i+1])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.savefig('Kd_trace_hist.png', dpi=500, bbox_inches='tight')
plt.savefig('Kd_trace_hist.pdf', bbox_inches='tight')

    
#Finally we want to plot our F_PL and F_L traces

interval_FPL = np.percentile(a=data['F_PL'][0], q=[2.5, 50.0, 97.5])

hist_legend = mpatches.Patch(color='C0',
    label = '$F\_PL$ = %.1e [%.1e,%.1e] '
    %(interval_FPL[1],interval_FPL[0],interval_FPL[2]) )

binBoundaries = np.linspace(1e10,5e11,70)

f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])

ax1 = plt.subplot(gs[0])
ax1.plot(range(0,len(data['F_PL'][0]),10),data['F_PL'][0][::10],color='C0')
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$F\_PL$',fontsize=16);
ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
ax1.set_ylim((1e10,2.5e11))
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,99000)
ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)

ax2 = plt.subplot(gs[1])
n, bins, patches = ax2.hist(data['F_PL'][0],color='C0',bins=binBoundaries,edgecolor='white',orientation="horizontal")
ax2.axhline(y=interval_FPL[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_FPL[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_FPL[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.set_xlabel('$P(F\_PL)$',fontsize=16);
ax2.set_ylim((1e10,2.5e11))
plt.xticks([])
plt.yticks([])

#set colors for 95% interval
clrs = ['C0' for xx in bins]
idxs = bins.argsort()
idxs = idxs[::-1]
gray_before = idxs[bins[idxs] < interval_FPL[0]]
gray_after = idxs[bins[idxs] > interval_FPL[2]]
for idx in gray_before:
    clrs[idx] = (.5,.5,.5)
for idx in gray_after:
    clrs[idx] = (.5,.5,.5)

for i,p in enumerate(patches):
    plt.setp(p, 'facecolor', clrs[i+1])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.savefig('FPL_trace_hist.png', dpi=500, bbox_inches='tight')
plt.savefig('FPL_trace_hist.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(9,3))

plt.plot(range(0,len(data['F_PL'][0]),10),data['F_PL'][0][::10])

plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('MCMC Sample', fontsize=18);
plt.ylabel('$F\_PL$', fontsize=20);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('F_PL.pdf', bbox_inches='tight')
plt.savefig('F_PL.png',dpi=500, bbox_inches='tight')

interval_FL = np.percentile(a=data['F_L'][0], q=[2.5, 50.0, 97.5])

hist_legend = mpatches.Patch(color='C1',
    label = '$F\_L$ = %.1e [%.1e,%.1e] '
    %(interval_FL[1],interval_FL[0],interval_FL[2]) )

binBoundaries = np.linspace(0.5e9,5e9,70)

f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])

ax1 = plt.subplot(gs[0])
ax1.plot(range(0,len(data['F_L'][0]),10),data['F_L'][0][::10],color='C1')
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$F\_L$',fontsize=16);
ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
ax1.set_ylim((0.8e9,3e9))
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,99000)
ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)

ax2 = plt.subplot(gs[1])
n, bins, patches = ax2.hist(data['F_L'][0],bins=binBoundaries,color='C1',edgecolor='white',orientation="horizontal")
ax2.axhline(y=interval_FL[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_FL[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval_FL[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.set_xlabel('$P(F\_L)$',fontsize=16);
ax2.set_ylim((0.8e9,3e9))
plt.xticks([])
plt.yticks([])

#set colors for 95% interval
clrs = ['C1' for xx in bins]
idxs = bins.argsort()
idxs = idxs[::-1]
gray_before = idxs[bins[idxs] < interval_FL[0]]
gray_after = idxs[bins[idxs] > interval_FL[2]]
for idx in gray_before:
    clrs[idx] = (.5,.5,.5)
for idx in gray_after:
    clrs[idx] = (.5,.5,.5)

for i,p in enumerate(patches):
    plt.setp(p, 'facecolor', clrs[i+1])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.savefig('FL_trace_hist.png', dpi=500, bbox_inches='tight')
plt.savefig('FL_trace_hist.pdf', bbox_inches='tight')

fig, ax = plt.subplots(figsize=(9,3))

plt.plot(range(0,len(data['F_L'][0]),10),data['F_L'][0][::10],color='green')

plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('MCMC Sample', fontsize=18);
plt.ylabel('$F\_L$', fontsize=20);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig('F_L.pdf', bbox_inches='tight')
plt.savefig('F_L.png',dpi=500, bbox_inches='tight')
