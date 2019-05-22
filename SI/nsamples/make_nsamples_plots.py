from glob import glob
from assaytools import parser
import numpy as np
import seaborn as sns
cols = sns.color_palette("deep", 3)
sns.set(style='white')
sns.set_context('talk')
import matplotlib.pyplot as plt

import cPickle

import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import gridspec

# Define nthin for this and all subsequent plots
nthin = 10

def plot_me(data,interval,my_color,output_file_single,ylim):
    
    [hist,bin_edges] = np.histogram(data['DeltaG'][0],bins=40,normed=True)
    binwidth = np.abs(bin_edges[0]-bin_edges[1])
    #set colors for 95% interval
    clrs = [my_color for xx in bin_edges]
    idxs = bin_edges.argsort()
    idxs = idxs[::-1]
    gray_before = idxs[bin_edges[idxs] < interval[0]]
    gray_after = idxs[bin_edges[idxs] > interval[2]]
    for idx in gray_before:
        clrs[idx] = (.5,.5,.5)
    for idx in gray_after:
        clrs[idx] = (.5,.5,.5)
    hist_legend = mpatches.Patch(color=my_color,
    label = '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$'
    %(interval[1],interval[0],interval[2]) )

    f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])

    ax1 = plt.subplot(gs[0])
    ax1.plot(range(0,len(data['DeltaG'][0]),nthin),data['DeltaG'][0][::nthin],color=my_color)
    ax1.set_xlabel('MCMC sample',fontsize=16);
    ax1.set_ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
    ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
    ax1.tick_params(labelsize=16)
    #ax1.set_xlim(0,len(data['DeltaG'][0]))
    ax1.set_xlim(0,250000)
    ax1.set_ylim(ylim)
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
    ax2.set_ylim(ylim)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.savefig('%s.png'%output_file_single, dpi=500, bbox_inches='tight')
    plt.savefig('%s.pdf'%output_file_single, bbox_inches='tight')

def plot_all_repeats_two_files(data_files_1,data_files_2, my_colors, output_file_base,ylim):

    #data_files_1 are the long files
    #data_files_2 are the short files
    
    #note this only supports when data_files and my_colors are size 3
    
    with open(r'%s'%data_files_1[0],'rb') as my_file:
        data1 = cPickle.load(my_file)
    with open(r'%s'%data_files_1[1],'rb') as my_file:
        data2 = cPickle.load(my_file)
    with open(r'%s'%data_files_1[2],'rb') as my_file:
        data3 = cPickle.load(my_file)
    
    interval1 = np.percentile(a=data1['DeltaG'][0], q=[2.5, 50.0, 97.5])
    interval2 = np.percentile(a=data2['DeltaG'][0], q=[2.5, 50.0, 97.5])
    interval3 = np.percentile(a=data3['DeltaG'][0], q=[2.5, 50.0, 97.5])
    
    plot_me(data1,interval1,my_colors[0],output_file_base+'_0',ylim)
    plot_me(data2,interval2,my_colors[1],output_file_base+'_1',ylim)
    plot_me(data3,interval3,my_colors[2],output_file_base+'_2',ylim)

    plt.clf()
    plt.figure(figsize=(3,3))
    plt.hist(data1['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[0], normed=True, orientation="horizontal")
    plt.hist(data2['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[1], normed=True, orientation="horizontal")
    plt.hist(data3['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[2], normed=True, orientation="horizontal")
    plt.ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
    plt.xlabel('$P(\Delta G)$',fontsize=16);
    plt.xticks([])
    plt.ylim(ylim)
    plt.savefig('%s_DeltaG_5m.png'%output_file_base, dpi=500, bbox_inches='tight')

    with open(r'%s'%data_files_2[0],'rb') as my_file:
        data1_short = cPickle.load(my_file)
    with open(r'%s'%data_files_2[1],'rb') as my_file:
        data2_short = cPickle.load(my_file)
    with open(r'%s'%data_files_2[2],'rb') as my_file:
        data3_short = cPickle.load(my_file)
    
    interval1_short = np.percentile(a=data1_short['DeltaG'][0], q=[2.5, 50.0, 97.5])
    interval2_short = np.percentile(a=data2_short['DeltaG'][0], q=[2.5, 50.0, 97.5])
    interval3_short = np.percentile(a=data3_short['DeltaG'][0], q=[2.5, 50.0, 97.5])
    
    plot_me(data1_short,interval1_short,my_colors[0],output_file_base+'2m_0',ylim)
    plot_me(data2_short,interval2_short,my_colors[1],output_file_base+'2m_1',ylim)
    plot_me(data3_short,interval3_short,my_colors[2],output_file_base+'2m_2',ylim)
    
    plt.clf()
    plt.figure(figsize=(3,3))
    plt.hist(data1_short['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[0], normed=True, orientation="horizontal")
    plt.hist(data2_short['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[1], normed=True, orientation="horizontal")
    plt.hist(data3_short['DeltaG'][0],40, alpha=0.6, edgecolor='white', color=my_colors[2], normed=True, orientation="horizontal")
    plt.ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
    plt.xlabel('$P(\Delta G)$',fontsize=16);
    plt.xticks([])
    plt.ylim(ylim)
    plt.title('100000')
    plt.savefig('%s_DeltaG_2m.png'%output_file_base, dpi=500, bbox_inches='tight')

    y_G2 = np.array([0,1,2])
    Gef_2m_DelG = np.array([interval1_short[1],interval2_short[1],interval3_short[1]])
    Gef_2m_DelG_error = [np.array([abs(interval1_short[1]-interval1_short[0]),abs(interval2_short[1]-interval2_short[0]),abs(interval3_short[1]-interval3_short[0])]),
                         np.array([abs(interval1_short[1]-interval1_short[2]),abs(interval2_short[1]-interval2_short[2]),abs(interval3_short[1]-interval3_short[2])])]
    
    y_G5 = np.array([3,4,5])
    Gef_5m_DelG = np.array([interval1[1],interval2[1],interval3[1]])
    Gef_5m_DelG_error = [np.array([abs(interval1[1]-interval1[0]),abs(interval2[1]-interval2[0]),abs(interval3[1]-interval3[0])]),
                         np.array([abs(interval1[1]-interval1[2]),abs(interval2[1]-interval2[2]),abs(interval3[1]-interval3[2])])]
    
    plt.clf()
    plt.figure(figsize=(6,3))
    plt.errorbar(Gef_2m_DelG,y_G2,xerr=Gef_2m_DelG_error, fmt='o',color='c',label='Src:Gef w/ iter=2e6')
    plt.errorbar(Gef_5m_DelG,y_G5,xerr=Gef_5m_DelG_error, fmt='o',color='b',label='Src:Gef w/ iter=5e6')
    plt.ylim((-1,8))
    plt.yticks([])
    plt.xlabel('$\Delta G$ ($k_B T$)', fontsize=16)
    plt.xticks(fontsize=12)
    plt.legend(fontsize=12, loc=0)
    plt.savefig('%s_bar_compare.png'%output_file_base, dpi=500, bbox_inches='tight')

SMH_purple = (0.7372549019607844, 0.5098039215686274, 0.7411764705882353)

three_colors = [SMH_purple,'C0','C5']

data_files = ['./5e6/Src-Gefitinib-GH_mcmc-0.pickle','./5e6/Src-Gefitinib-GH_mcmc-1.pickle','./5e6/Src-Gefitinib-GH_mcmc-2.pickle']
data_files_short = ['./2e6/Src-Gefitinib-GH_mcmc-0.pickle','./2e6/Src-Gefitinib-GH_mcmc-1.pickle','./2e6/Src-Gefitinib-GH_mcmc-2.pickle']

Gef_lim = [-13.5,-9.5]

plot_all_repeats_two_files(data_files, data_files_short, three_colors,'Src-Gefitinib',Gef_lim)

Erl_data_files = ['./5e6/Src-Erlotinib-EF_mcmc-0.pickle','./5e6/Src-Erlotinib-EF_mcmc-1.pickle','./5e6/Src-Erlotinib-EF_mcmc-2.pickle']
Erl_data_files_short = ['./2e6/Src-Erlotinib-EF_mcmc-0.pickle','./2e6/Src-Erlotinib-EF_mcmc-1.pickle','./2e6/Src-Erlotinib-EF_mcmc-2.pickle']

Erl_lim = [-35,-10]

plot_all_repeats_two_files(Erl_data_files, Erl_data_files_short, three_colors,'Src-Erlotinib',Erl_lim)

    