import matplotlib
matplotlib.use('Agg')

import numpy as np
from assaytools import parser
import string
from glob import glob
from assaytools import pymcmodels

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import seaborn as sns

import traceback
import cPickle

sns.set(style='white')
sns.set_context('talk')

def mcmc_three_plots(pymc_model,mcmc,Lstated,name,iteration,identifier):
    
    import pymbar
    [t,g,Neff_max] = pymbar.timeseries.detectEquilibration(mcmc.DeltaG.trace())
    
    interval= np.percentile(a=mcmc.DeltaG.trace()[t:], q=[2.5, 50.0, 97.5])
    [hist,bin_edges] = np.histogram(mcmc.DeltaG.trace()[t:],bins=40,normed=True)
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
    
    plt.clf();
    plt.figure(figsize=(12,3));

    plt.subplot(131)
    property_name = 'top_complex_fluorescence'
    complex = getattr(pymc_model, property_name)
    property_name = 'top_ligand_fluorescence'
    ligand = getattr(pymc_model, property_name)
    for top_complex_fluorescence_model in mcmc.top_complex_fluorescence_model.trace()[::10]:
        plt.semilogx(Lstated, top_complex_fluorescence_model, marker='.',color='silver')
    for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::10]:
        plt.semilogx(Lstated, top_ligand_fluorescence_model, marker='.',color='lightcoral', alpha=0.2)
    plt.semilogx(Lstated, complex.value, 'ko',label='complex')
    plt.semilogx(Lstated, ligand.value, marker='o',color='firebrick',linestyle='None',label='ligand')
    plt.xlim(.5e-8,5e-5)
    plt.xlabel('$[L]_T$ (M)');
    plt.yticks([])
    plt.ylabel('fluorescence');
    plt.legend(loc=0);

    plt.subplot(132)
    plt.bar(bin_edges[:-1]+binwidth/2,hist,binwidth,color=clrs, edgecolor = "white");
    sns.kdeplot(mcmc.DeltaG.trace()[t:],bw=.4,color=(0.39215686274509803, 0.7098039215686275, 0.803921568627451),shade=False)
    plt.axvline(x=interval[0],color=(0.5,0.5,0.5),linestyle='--')
    plt.axvline(x=interval[1],color=(0.5,0.5,0.5),linestyle='--')
    plt.axvline(x=interval[2],color=(0.5,0.5,0.5),linestyle='--')
    plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16);
    plt.ylabel('$P(\Delta G)$',fontsize=16);
    plt.xlim(-20,-8)
    hist_legend = mpatches.Patch(color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353), 
        label = '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$' 
        %(interval[1],interval[0],interval[2]) )
    plt.legend(handles=[hist_legend],fontsize=10,loc=0,frameon=True);

    plt.subplot(133)
    plt.plot(range(0,t),mcmc.DeltaG.trace()[:t], 'g.',label='equil. at %s'%t)
    plt.plot(range(t,len(mcmc.DeltaG.trace())),mcmc.DeltaG.trace()[t:], '.')
    plt.xlabel('MCMC sample');
    plt.ylabel('$\Delta G$ ($k_B T$)');
    plt.legend(loc=2);

    plt.tight_layout();
    plt.suptitle('%s binding'%(name),y=1.01);

    plt.savefig('%s_binding_iter%s-%s.png'%(name,iteration,identifier),dpi=1000)
    
    return [t,interval,hist,bin_edges,binwidth,clrs]
    
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

dPstated = 0.35 * inputs['Pstated'] # protein concentration uncertainty
dLstated = 0.08 * inputs['Lstated'] # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)

names = ['Src-Erlotinib-EF']

# Here we are importing the mean value of F_PL from a previous run of assaytools on
# a different Src-Erlotinib experiment from 2016-03-09 using quickmodel, this is
# also analyzed as spectra_2gains in the analysis folder.

data_file = '../../analysis/bayes/spectra_2gains/Src-Erlotinib-EF_mcmc-2017-04-13 14:54.pickle'
# load data from pickle file
with open(r'%s'%data_file,'rb') as my_file:
    data = cPickle.load(my_file)

F_PL_mean = data['F_PL'][0].mean()
F_PL_def = {'Src-Erlotinib-EF': F_PL_mean}
print(F_PL_def)

iterations = 1

#First with no FPL and dPstated of 0.35
t = {}
for name in names:

    for i in range(iterations):

        print('Working on %s iteration %s.' %(name,i))

        pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
           top_complex_fluorescence=complex_fluorescence[name],
           top_ligand_fluorescence=ligand_fluorescence[name],
           use_primary_inner_filter_correction=True,
           use_secondary_inner_filter_correction=True,
           assay_volume=inputs['assay_volume'], DG_prior='uniform')
   
        #mcmc = pymcmodels.run_mcmc(pymc_model, niter=1000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) #just for testing
        mcmc = pymcmodels.run_mcmc(pymc_model, niter=2000000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) 
        
        [t,interval,hist,bin_edges,binwidth,clrs] = mcmc_three_plots(pymc_model,mcmc,inputs['Lstated'],name,i,'normal')

#next with F_PL defined
t = {}
for name in names:

    for i in range(iterations):

        print('Working on %s iteration %s.' %(name,i))

        pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
           top_complex_fluorescence=complex_fluorescence[name],
           top_ligand_fluorescence=ligand_fluorescence[name],
           use_primary_inner_filter_correction=True,
           use_secondary_inner_filter_correction=True,
           assay_volume=inputs['assay_volume'], DG_prior='uniform',
           F_PL=F_PL_def[name],dF_PL=0.50)
   
        #mcmc = pymcmodels.run_mcmc(pymc_model, niter=1000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) #just for testing
        mcmc = pymcmodels.run_mcmc(pymc_model, niter=2000000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) 
        
        [t,interval,hist,bin_edges,binwidth,clrs] = mcmc_three_plots(pymc_model,mcmc,inputs['Lstated'],name,i,'F_PL')

#next with dPstated error of 15%
dPstated = 0.15 * inputs['Pstated'] # protein concentration uncertainty

t = {}
for name in names:

    for i in range(iterations):

        print('Working on %s iteration %s.' %(name,i))

        pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
           top_complex_fluorescence=complex_fluorescence[name],
           top_ligand_fluorescence=ligand_fluorescence[name],
           use_primary_inner_filter_correction=True,
           use_secondary_inner_filter_correction=True,
           assay_volume=inputs['assay_volume'], DG_prior='uniform')
   
        #mcmc = pymcmodels.run_mcmc(pymc_model, niter=1000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) #just for testing
        mcmc = pymcmodels.run_mcmc(pymc_model, niter=2000000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) 
        
        [t,interval,hist,bin_edges,binwidth,clrs] = mcmc_three_plots(pymc_model,mcmc,inputs['Lstated'],name,i,'dP_stated')

#lastly with both modifications
t = {}
for name in names:

    for i in range(iterations):

        print('Working on %s iteration %s.' %(name,i))

        pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
           top_complex_fluorescence=complex_fluorescence[name],
           top_ligand_fluorescence=ligand_fluorescence[name],
           use_primary_inner_filter_correction=True,
           use_secondary_inner_filter_correction=True,
           assay_volume=inputs['assay_volume'], DG_prior='uniform',
           F_PL=F_PL_def[name],dF_PL=0.50)
   
        #mcmc = pymcmodels.run_mcmc(pymc_model, niter=1000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) #just for testing
        mcmc = pymcmodels.run_mcmc(pymc_model, niter=2000000, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,i)) 
        
        [t,interval,hist,bin_edges,binwidth,clrs] = mcmc_three_plots(pymc_model,mcmc,inputs['Lstated'],name,i,'both')