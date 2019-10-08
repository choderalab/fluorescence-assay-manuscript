# adapted from Sonya's script for Figure 1 in the paper
# EEG, July 2019

from assaytools import platereader
from glob import glob
from lxml import etree
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

#This function allows us to import a section from an xml formated data file and convert it to a pandas dataframe
def xml2df_section(file,section):

    root = etree.parse(file)

    data = []

    reads = root.xpath("/*/Section[%s]/*/Well"%section)
    Sections = root.xpath("/*/Section")
    #print Sections[(section-1)].attrib['Name']
    section_name = Sections[(section-1)].attrib['Name']

    wellIDs = [read.attrib['Pos'] for read in reads]

    data = [(s.text, float(s.attrib['WL']), r.attrib['Pos'])
        for r in reads
        for s in r]

    dataframe = pd.DataFrame(data, columns=['fluorescence','wavelength (nm)','Well'])

    ### dataframe_rep replaces 'OVER' (when fluorescence signal maxes out) with '3289277', an arbitrarily high number

    #dataframe_rep = dataframe.replace({'OVER':'1.2e5'})

    dataframe_rep = dataframe.replace({'OVER':'NAN'})

    dataframe_rep[['fluorescence']] = dataframe_rep[['fluorescence']].astype('float')

    dataframe_pivot = pd.pivot_table(dataframe_rep, index = 'wavelength (nm)', columns = ['Well'])

    #Rearrange columns so they're in the right order
    cols =  dataframe_pivot['fluorescence'].columns.tolist()
    cols = [cols[0]] + cols[4:12] + cols[1:4] + [cols[12]] + cols[16:24] + cols[13:16]
    dataframe_reindex =  dataframe_pivot.reindex_axis(cols,level='Well',axis=1)

    return [dataframe_reindex,section_name]

#This function allows us to plot spectra choosing ylim and
def plot_spectra_grid_advanced_inset(file,protein,ligands,ligand,section,ylim,lines,Lstated):
    grid = len(protein) + len(ligand)

    # file = file_set[protein]
    file = file

    # make a dataframe
    [df,section_name] = xml2df_section(file,section)

    # pick a title
    title = "%s - %s: %s" %(protein, ligand, section_name)
    print(title)

    # define ylim and lines

    ylim = ylim
    lines = lines

    # plot the spectra
    fig = plt.figure(figsize=(7,6));
    ax = df['fluorescence'].iloc[:,12].plot(ylim=(-10,ylim),legend=False, linewidth=4,color='m',logy=True);
    ax.axvline(x=lines[0],color='black',linestyle='--');
    ax.axvline(x=lines[1],color='blue',linestyle='--');
    for i in range(11):
        df['fluorescence'].iloc[:,i].plot(linewidth=3,c=cm.hsv(i*15), ax = ax,logy=True);
        df['fluorescence'].iloc[:,11+i].plot(legend=False, linewidth=4,c=cm.gray(i*15+50),logy=True,ax = ax, fontsize =20);
    sns.despine()
    plt.yticks([])
    plt.xticks(fontsize=26)
    plt.xlim(310,600)
    plt.xlabel('wavelength (nm)', fontsize=30)
    plt.ylabel('log(Intensity)', fontsize=30)
#    plt.text(550,0.9*ylim,"lines=%s"%lines,color='0.7')
#    plt.title(title)
    plt.tight_layout();

    ## Plot intensity from `lines` wavelengths

    complex_280_340 = df['fluorescence'].loc[lines[0]][:12]
    ligand_280_340 = df['fluorescence'].loc[lines[0]][12:]

    complex_280_480 = df['fluorescence'].loc[lines[1]][:12]
    ligand_280_480 = df['fluorescence'].loc[lines[1]][12:]

    difference_280_480 = complex_280_480.values - ligand_280_480.values
    difference_280_340 = complex_280_340.values - ligand_280_340.values

    difference_280_480_normalized = difference_280_480/np.nanmax(difference_280_480)
    difference_280_340_normalized = difference_280_340/np.nanmax(difference_280_340)

    a = plt.axes([0.7, 0.7, .25, .25])
    plt.semilogx(Lstated,difference_280_480_normalized,color="blue",marker='o',linestyle='None',label='%s nm'%lines[1]);
    plt.semilogx(Lstated,difference_280_340_normalized,color="black",marker='o',linestyle='None',label='%s nm'%lines[0]);
    plt.yticks([])
    x_inset_labels = [0,-8,-7,-6,-5]
    plt.xticks([1e-9,1e-8,1e-7,1e-6,1e-5], x_inset_labels,fontsize=16)
    plt.xlabel('$log_{10}([L])$', fontsize=18);
    plt.ylabel('RFU', fontsize=18);
    plt.legend(loc=0,bbox_to_anchor=(-0.6 , 1.1), fontsize=20);
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)
    a.spines['left'].set_visible(False)
    a.spines['left'].set_smart_bounds(True)

Lstated = [20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0035e-6]

ylim = 500000
lines = [340,480]


file = "infinite_results/Abl_T334I/2019-07-16/Abl_D382N_T334I_Bos_UV_20190716_135721.xml"
ligands = ['Bosutinib']
plot_spectra_grid_advanced_inset(file,'Abl T334I',ligands,'Bosutinib',1,ylim,lines,Lstated)
plt.savefig('Abl T334I-Bosutinib-log-inset.png',dpi=500)
plt.savefig('Abl T334I-Bosutinib-log-inset.pdf')

file = "infinite_results/Abl_T334I/2019-07-16/Abl_D382N_T334I_Bos_Iso_UV_20190716_143313.xml"
ligands = ['Bosutinib Isomer']
plot_spectra_grid_advanced_inset(file,'Abl T334I',ligands,'Bosutinib Isomer',1,ylim,lines,Lstated)
plt.savefig('Abl T334I-Bosutinib Isomer-log-inset.png',dpi=500)
plt.savefig('Abl T334I-Bosutinib Isomer-log-inset.pdf')

file = "infinite_results/Abl_T334I/2019-07-16/Abl_D382N_T334I_Erl_UV_20190716_151055.xml"
ligands = ['Erlotinib']
plot_spectra_grid_advanced_inset(file,'Abl T334I',ligands,'Erlotinib',1,ylim,lines,Lstated)
plt.savefig('Abl T334I-Erlotinib-log-inset.png',dpi=500)
plt.savefig('Abl T334I-Erlotinib-log-inset.pdf')

file = "infinite_results/Abl_T334I/2019-07-16/Abl_D382N_T334I_Gef_UV_20190716_154656.xml"
ligands = ['Gefitinib']
plot_spectra_grid_advanced_inset(file,'Abl T334I',ligands,'Gefitinib',1,ylim,lines,Lstated)
plt.savefig('Abl T334I-Gefitinib-log-inset.png',dpi=500)
plt.savefig('Abl T334I-Gefitinib-log-inset.pdf')


file = "infinite_results/Src_T338I/2019-07-16/Src_T338I_Bos_NB_20190715_125826.xml"
ligands = ['Bosutinib']
plot_spectra_grid_advanced_inset(file,'Src_T338I',ligands,'Bosutinib',1,ylim,lines,Lstated)
plt.savefig('Src_T338I-Bosutinib-log-inset.png',dpi=500)
plt.savefig('Src_T338I-Bosutinib-log-inset.pdf')

file = "infinite_results/Src_T338I/2019-07-16/Src_T338I_Bos_Iso_NB_20190715_131151.xml"
ligands = ['Bosutinib Isomer']
plot_spectra_grid_advanced_inset(file,'Src_T338I',ligands,'Bosutinib Isomer',1,ylim,lines,Lstated)
plt.savefig('Src_T338I-Bosutinib Isomer-log-inset.png',dpi=500)
plt.savefig('Src_T338I-Bosutinib Isomer-log-inset.pdf')

file = "infinite_results/Src_T338I/2019-07-16/Src_T338I_Erl_NB_20190715_132525.xml"
ligands = ['Erlotinib']
plot_spectra_grid_advanced_inset(file,'Src_T338I',ligands,'Erlotinib',1,ylim,lines,Lstated)
plt.savefig('Src_T338I-Erlotinib-log-inset.png',dpi=500)
plt.savefig('Src_T338I-Erlotinib-log-inset.pdf')

file = "infinite_results/Src_T338I/2019-07-16/Src_T338I_Gef_NB_20190715_133852.xml"
ligands = ['Gefitinib']
plot_spectra_grid_advanced_inset(file,'Src_T338I',ligands,'Gefitinib',1,ylim,lines,Lstated)
plt.savefig('Src_T338I-Gefitinib-log-inset.png',dpi=500)
plt.savefig('Src_T338I-Gefitinib-log-inset.pdf')
