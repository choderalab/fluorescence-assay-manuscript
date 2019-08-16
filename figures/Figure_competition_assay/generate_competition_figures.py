import pymc
from glob import glob
import numpy as np

DG_min = np.log(1e-15) # kT, most favorable (negative) binding free energy possible; 1 fM
DG_max = +0 # kT, least favorable binding free energy possible

#this is just for the lognormal wrapper and inner filter effect and run_mcmc
from assaytools import pymcmodels
from assaytools import bindingmodels

from assaytools import parser
import matplotlib.pyplot as plt

#This is our new competition assay model, which has competition as an option (note just for top fluroescence)
#if the values on the second line are included
def make_comp_model(Pstated, dPstated, Lstated, dLstated,
               Bstated = None, dBstated = None, DG_L_mean = None, DG_L_std = None, # specific for competition assay
               top_complex_fluorescence=None, top_ligand_fluorescence=None,
               DG_prior='uniform',
               concentration_priors='lognormal',
               quantum_yield_priors='lognormal',
               use_primary_inner_filter_correction=True,
               use_secondary_inner_filter_correction=True,
               assay_volume=100e-6, well_area=0.1586,
               epsilon_ex=None, depsilon_ex=None,
               epsilon_em=None, depsilon_em=None,
               ligand_ex_absorbance=None, ligand_em_absorbance=None,
               F_PL=None, dF_PL=None):

    # Compute path length.
    path_length = assay_volume * 1000 / well_area # cm, needed for inner filter effect corrections

    # Compute number of samples.
    N = len(Lstated)
    
    # Check input.
    # TODO: Check fluorescence and absorbance measurements for correct dimensions.
    if (len(Pstated) != N):
        raise Exception('len(Pstated) [%d] must equal len(Lstated) [%d].' % (len(Pstated), len(Lstated)))
    if (len(dPstated) != N):
        raise Exception('len(dPstated) [%d] must equal len(Lstated) [%d].' % (len(dPstated), len(Lstated)))
    if (len(dLstated) != N):
        raise Exception('len(dLstated) [%d] must equal len(Lstated) [%d].' % (len(dLstated), len(Lstated)))

    # Note whether we have top or bottom fluorescence measurements.
    top_fluorescence = (top_complex_fluorescence is not None) or (top_ligand_fluorescence is not None) # True if any top fluorescence measurements provided

    # Create an empty dict to hold the model.
    model = dict()
    
    # Prior on binding free energies.
    if DG_prior == 'uniform':
        DeltaG_B = pymc.Uniform('DeltaG_B', lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
    elif DG_prior == 'chembl':
        DeltaG_B = pymc.Normal('DeltaG_B', mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
    else:
        raise Exception("DG_prior = '%s' unknown. Must be one of 'uniform' or 'chembl'." % DG_prior)
    model['DeltaG_B'] = DeltaG_B
    
    # Prior on known binding free energy.
    if DG_L_mean == None:
        DeltaG_L = pymc.Uniform('DeltaG_L', lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
    else:
        DeltaG_L = pymc.Normal('DeltaG_L', mu=DG_L_mean, tau=DG_L_std) 
    model['DeltaG_L'] = DeltaG_L
        
    if concentration_priors != 'lognormal':
        raise Exception("concentration_priors = '%s' unknown. Must be one of ['lognormal']." % concentration_priors)
    model['log_Ptrue'], model['Ptrue'] = pymcmodels.LogNormalWrapper('Ptrue', mean=Pstated, stddev=dPstated, size=Pstated.shape) # protein concentration (M)
    model['log_Ltrue'], model['Ltrue'] = pymcmodels.LogNormalWrapper('Ltrue', mean=Lstated, stddev=dLstated, size=Lstated.shape) # ligand concentration (M)
    model['log_Btrue'], model['Btrue'] = pymcmodels.LogNormalWrapper('Btrue', mean=Bstated, stddev=dBstated, size=Bstated.shape) # ligand concentration (M)
    model['log_Ltrue_control'], model['Ltrue_control'] = pymcmodels.LogNormalWrapper('Ltrue_control', mean=Lstated, stddev=dLstated, size=Lstated.shape) # ligand concentration in ligand-only wells (M)

    # extinction coefficient
    model['epsilon_ex'] = 0.0
    if use_primary_inner_filter_correction:
        if epsilon_ex:
            model['log_epsilon_ex'], model['epsilon_ex'] =  pymcmodels.LogNormalWrapper('epsilon_ex', mean=epsilon_ex, stddev=depsilon_ex) # prior is centered on measured extinction coefficient
        else:
            model['epsilon_ex'] = pymc.Uniform('epsilon_ex', lower=0.0, upper=1000e3, value=70000.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm

    model['epsilon_em'] = 0.0
    if use_secondary_inner_filter_correction:
        if epsilon_em:
            model['log_epsilon_em'], model['epsilon_em'] =  pymcmodels.LogNormalWrapper('epsilon_em', mean=epsilon_em, stddev=depsilon_em) # prior is centered on measured extinction coefficient
        else:
            model['epsilon_em'] = pymc.Uniform('epsilon_em', lower=0.0, upper=1000e3, value=0.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm

    # Min and max observed fluorescence.
    Fmax = 0.0; Fmin = 1e6;
    if top_complex_fluorescence is not None:
        Fmax = max(Fmax, top_complex_fluorescence.max()); Fmin = min(Fmin, top_complex_fluorescence.min())
    if top_ligand_fluorescence is not None:
        Fmax = max(Fmax, top_ligand_fluorescence.max()); Fmin = min(Fmin, top_ligand_fluorescence.min())

    # Compute initial guesses for fluorescence quantum yield quantities.
    F_plate_guess = Fmin
    F_buffer_guess = Fmin / path_length
    F_L_guess = (Fmax - Fmin) / Lstated.max()
    F_P_guess = Fmin
    F_P_guess = Fmin / Pstated.min()
    F_PL_guess = (Fmax - Fmin) / min(Pstated.max(), Lstated.max())

    # Priors on fluorescence intensities of complexes (later divided by a factor of Pstated for scale).

    if quantum_yield_priors == 'lognormal':
        stddev = 1.0 # relative factor for stddev guess
        model['log_F_plate'], model['F_plate'] =  pymcmodels.LogNormalWrapper('F_plate', mean=F_plate_guess, stddev=stddev*F_plate_guess) # plate fluorescence
        model['log_F_buffer'], model['F_buffer'] =  pymcmodels.LogNormalWrapper('F_buffer', mean=F_buffer_guess, stddev=stddev*F_buffer_guess) # buffer fluorescence
        model['log_F_buffer_control'], model['F_buffer_control'] =  pymcmodels.LogNormalWrapper('F_buffer_control', mean=F_buffer_guess, stddev=stddev*F_buffer_guess) # buffer fluorescence
        if (F_PL is not None) and (dF_PL is not None):
            model['log_F_PL'], model['F_PL'] =  pymcmodels.LogNormalWrapper('F_PL', mean=F_PL, stddev=dF_PL)
        else:
            model['log_F_PL'], model['F_PL'] =  pymcmodels.LogNormalWrapper('F_PL', mean=F_PL_guess, stddev=stddev*F_PL_guess) # complex fluorescence
        model['log_F_P'], model['F_P'] =  pymcmodels.LogNormalWrapper('F_P', mean=F_P_guess, stddev=stddev*F_P_guess) # protein fluorescence
        model['log_F_L'], model['F_L'] =  pymcmodels.LogNormalWrapper('F_L', mean=F_L_guess, stddev=stddev*F_L_guess) # ligand fluorescence
    else:
        raise Exception("quantum_yield_priors = '%s' unknown. Must be one of ['lognormal', 'uniform']." % quantum_yield_priors)

    # Unknown experimental measurement error.
    if top_fluorescence:
        model['log_sigma_top'] = pymc.Uniform('log_sigma_top', lower=-10, upper=np.log(Fmax), value=np.log(5))
        model['sigma_top'] = pymc.Lambda('sigma_top', lambda log_sigma=model['log_sigma_top'] : np.exp(log_sigma) )
        model['precision_top'] = pymc.Lambda('precision_top', lambda log_sigma=model['log_sigma_top'] : np.exp(-2*log_sigma) )

    if top_fluorescence:
        model['log_sigma_abs'] = pymc.Uniform('log_sigma_abs', lower=-10, upper=0, value=np.log(0.01))
        model['sigma_abs'] = pymc.Lambda('sigma_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(log_sigma) )
        model['precision_abs'] = pymc.Lambda('precision_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(-2*log_sigma) )

    if top_complex_fluorescence is not None:
        @pymc.deterministic
        def top_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                           F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                           Ptrue=model['Ptrue'], Ltrue=model['Ltrue'], Btrue=model['Btrue'], DeltaG_L=model['DeltaG_L'], DeltaG_B=model['DeltaG_B'],
                                           epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            [P_i, L_i, PL_i, B_i, PB_i] = bindingmodels.CompetitionBindingModel.equilibrium_concentrations(Ptrue[:], Ltrue[:], DeltaG_L, Btrue[:], DeltaG_B)
            IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_complex_fluorescence_model'] = top_complex_fluorescence_model
        model['log_top_complex_fluorescence'], model['top_complex_fluorescence'] = pymcmodels.LogNormalWrapper('top_complex_fluorescence',
            mean=model['top_complex_fluorescence_model'], stddev=model['sigma_top'],
            size=[N], observed=True, value=top_complex_fluorescence) # observed data
        
    if top_ligand_fluorescence is not None:
        @pymc.deterministic
        def top_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer_control=model['F_buffer_control'],
                                          F_L=model['F_L'],
                                          Ltrue_control=model['Ltrue_control'],
                                          epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue_control, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*Ltrue_control) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue_control + F_buffer_control*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_ligand_fluorescence_model'] = top_ligand_fluorescence_model
        model['log_top_ligand_fluorescence'], model['top_ligand_fluorescence'] = pymcmodels.LogNormalWrapper('top_ligand_fluorescence',
                                                       mean=model['top_ligand_fluorescence_model'], stddev=model['sigma_top'],
                                                       size=[N], observed=True, value=top_ligand_fluorescence) # observed data
        
    #Below this is competition model stuff!    
    
    # Compute number of samples.
    if Bstated is not None:
        N_B = len(Bstated)

        # Check input.
        # TODO: Check fluorescence and absorbance measurements for correct dimensions.
        if (len(Lstated) != N_B):
            raise Exception('len(Lstated) [%d] must equal len(Bstated) [%d].' % (len(Lstated), len(Bstated)))
        if (len(Pstated) != N_B):
            raise Exception('len(Pstated) [%d] must equal len(Bstated) [%d].' % (len(Pstated), len(Bstated)))
        if (len(dPstated) != N_B):
            raise Exception('len(dPstated) [%d] must equal len(Bstated) [%d].' % (len(dPstated), len(Bstated)))
        if (len(dBstated) != N_B):
            raise Exception('len(dBstated) [%d] must equal len(Bstated) [%d].' % (len(dBstated), len(Bstated)))

            
    if Bstated is not None:
        @pymc.deterministic
        def top_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                           F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                           Ptrue=model['Ptrue'], Ltrue=model['Ltrue'], Btrue=model['Btrue'],
                                           DeltaG_L=model['DeltaG_L'], DeltaG_B=model['DeltaG_B'],
                                           epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            [P_i, L_i, PL_i, B_i, PB_i] = bindingmodels.CompetitionBindingModel.equilibrium_concentrations(Ptrue[:], Ltrue[:], DeltaG_L, Btrue[:], DeltaG_B)
            IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_complex_fluorescence_model'] = top_complex_fluorescence_model
        model['log_top_complex_fluorescence'], model['top_complex_fluorescence'] = pymcmodels.LogNormalWrapper('top_complex_fluorescence',
            mean=model['top_complex_fluorescence_model'], stddev=model['sigma_top'],
            size=[N], observed=True, value=top_complex_fluorescence) # observed data

        @pymc.deterministic
        def top_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer_control=model['F_buffer_control'],
                                          F_L=model['F_L'],
                                          Ltrue_control=model['Ltrue_control'],
                                          epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue_control, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*Ltrue_control) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue_control + F_buffer_control*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_ligand_fluorescence_model'] = top_ligand_fluorescence_model
        model['log_top_ligand_fluorescence'], model['top_ligand_fluorescence'] = pymcmodels.LogNormalWrapper('top_ligand_fluorescence',
                                                       mean=model['top_ligand_fluorescence_model'], stddev=model['sigma_top'],
                                                       size=[N], observed=True, value=top_ligand_fluorescence) # observed data
    # Promote this to a full-fledged PyMC model.
    pymc_model = pymc.Model(model)

    # Return the pymc model
    return pymc_model
	
	
inputs = {
    'xml_file_path' :  "../../data/",
    'file_set'      :  {'Abl-IMA': glob("../../data/single_wavelength/competition-assay/Abl*16-22-45_plate*.xml")},
    'ligand_order'  :  ['Gefitinib','Gefitinib','Gefitinib','Gefitinib'],
    'competitive_ligand'  :  'Imatinib',
    'section'       :  '280_BottomRead',
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Bstated'       :  10.0e-6 * np.ones([12],np.float64), # competitive ligand concentration, M
    'Pstated'       :  0.5e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf],
    'DG_L_mean'     :  -12.5, #DeltaG of Gefitinib mean estimate (kT)
    'DG_L_std'      :  1      #DeltaG of Gefitinib standard deviation estimate (kT)
    }

dPstated = 0.35 * inputs['Pstated'] # protein concentration uncertainty
dLstated = 0.08 * inputs['Lstated'] # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)
dBstated = 0.08 * inputs['Bstated'] # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)

[complex_fluorescence_comp, ligand_fluorescence_comp] = parser.get_data_using_inputs(inputs)  

name = 'Abl-IMA-Gefitinib-AB'

pymc_comp_model = make_comp_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
    inputs['Bstated'], dBstated, inputs['DG_L_mean'], inputs['DG_L_std'],# specific for competition assay
    top_complex_fluorescence=complex_fluorescence_comp[name],
    top_ligand_fluorescence=ligand_fluorescence_comp[name],
    use_primary_inner_filter_correction=True,
    use_secondary_inner_filter_correction=True,
    assay_volume=inputs['assay_volume'], DG_prior='uniform')

comp_mcmc = pymc.MCMC(pymc_comp_model)

nburn = 5000
niter = 200000
nthin = 20

comp_mcmc.sample(iter=(nburn+niter), burn=nburn, thin=nthin, progress_bar=False, tune_throughout=True)

# I want to be able to plot the plain Abl-Gef data on this.

inputs = {
    'xml_file_path' :  "../../data/",
    'file_set'      :  {'Abl': glob("../../data/single_wavelength/competition-assay/Abl*15-59-53_plate*.xml")},
    'ligand_order'  :  ['Gefitinib','Gefitinib','Gefitinib','Gefitinib'],
    'competitive_ligand'  :  'NONE',
    'section'       :  '280_BottomRead',
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'Bstated'       :  0.0 * np.ones([12],np.float64), # competitive ligand concentration, M
    'Pstated'       :  0.5e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf],
    'DG_L_mean'     :  -12.5, #DeltaG of Gefitinib mean estimate (kT)
    'DG_L_std'      :  1      #DeltaG of Gefitinib standard deviation estimate (kT)
    }
	
[complex_fluorescence_comp_noIMA, ligand_fluorescence_comp_noIMA] = parser.get_data_using_inputs(inputs)  
	
name_noIMA = 'Abl-Gefitinib-AB'
name_noIMA2 = 'Abl-Gefitinib-CD'
name_noIMA3 = 'Abl-Gefitinib-EF'
name_noIMA4 = 'Abl-Gefitinib-GH'

name2 = 'Abl-IMA-Gefitinib-CD'
name3 = 'Abl-IMA-Gefitinib-EF'
name4 = 'Abl-IMA-Gefitinib-GH'

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'], complex_fluorescence_comp_noIMA[name_noIMA],'ko',label='Abl:Gefitinib')

plt.semilogx(inputs['Lstated'], complex_fluorescence_comp[name],'ro',label='Abl:Gefitinib:Imatinib')

plt.semilogx(inputs['Lstated'],ligand_fluorescence_comp_noIMA[name_noIMA],'o',color='dimgrey',label='Gefitinib')

plt.semilogx(inputs['Lstated'],ligand_fluorescence_comp[name],'o',color='lightsalmon',label='Gefitinib:Imatinib')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlabel('$[L]_T$ (M)',fontsize=16);
plt.yticks([])
plt.xticks(fontsize=16)
plt.ylabel('Fluorescence',fontsize=20);
plt.legend()

plt.tight_layout()

plt.savefig('Abl_Gef_Ima_fig.png')

fig, ax = plt.subplots(figsize=(8,4))

for top_complex_fluorescence_model in comp_mcmc.top_complex_fluorescence_model.trace()[::10]:
    plt.semilogx(inputs['Lstated'], top_complex_fluorescence_model, color='tomato',marker='.', alpha=0.2, hold=True)

for top_ligand_fluorescence_model in comp_mcmc.top_ligand_fluorescence_model.trace()[::10]:
    plt.semilogx(inputs['Lstated'], top_ligand_fluorescence_model, color='peachpuff',marker='.', alpha=0.2, hold=True)
    
plt.semilogx(inputs['Lstated'], complex_fluorescence_comp[name],'o',color='darkred',label='complex')

plt.semilogx(inputs['Lstated'],ligand_fluorescence_comp[name],'o',color='lightsalmon',label='ligand')
    
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlabel('$[L]_T$ (M)',fontsize=16);
plt.yticks([])
plt.xticks(fontsize=16)
plt.ylabel('Fluorescence',fontsize=20);

plt.tight_layout()

plt.savefig('Ima_trace.png', dpi=500, bbox_inches='tight')

import matplotlib.patches as mpatches
import matplotlib.lines as mlines

interval = np.percentile(a=comp_mcmc.DeltaG_B.trace(), q=[2.5, 50.0, 97.5])
[hist,bin_edges] = np.histogram(comp_mcmc.DeltaG_B.trace(),bins=40,normed=True)
binwidth = np.abs(bin_edges[0]-bin_edges[1])

#set colors for 95% interval
clrs = [('tomato') for xx in bin_edges]
idxs = bin_edges.argsort()
idxs = idxs[::-1]
gray_before = idxs[bin_edges[idxs] < interval[0]]
gray_after = idxs[bin_edges[idxs] > interval[2]]
for idx in gray_before:
    clrs[idx] = (.5,.5,.5)
for idx in gray_after:
    clrs[idx] = (.5,.5,.5)

hist_legend = mpatches.Patch(color=('tomato'),
    label = '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$'
    %(interval[1],interval[0],interval[2]) )
	
f, (ax1, ax2) = plt.subplots(1,2, sharey=True,figsize=(10,3))
from matplotlib import gridspec
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 

ax1 = plt.subplot(gs[0])
ax1.plot(range(0,len(comp_mcmc.DeltaG_B.trace())),comp_mcmc.DeltaG_B.trace(),color=('tomato'))
ax1.set_xlabel('MCMC sample',fontsize=16);
ax1.set_ylabel('$\Delta G$ ($k_B T$)',fontsize=16);
ax1.legend(handles=[hist_legend],fontsize=14,loc=4,frameon=True)
ax1.tick_params(labelsize=16)
ax1.set_xlim(0,11000)

ax1.spines['top'].set_visible(False)

f.subplots_adjust(wspace=0)

ax2 = plt.subplot(gs[1])
ax2.barh(bin_edges[:-1],hist,binwidth,color=clrs, edgecolor = "white");
ax2.axhline(y=AblIma_dG,color='magenta')
ax2.axhline(y=interval[0],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[1],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.axhline(y=interval[2],color=(0.5,0.5,0.5),linestyle='--')
ax2.set_xlabel('$P(\Delta G)$',fontsize=16);
plt.xticks([])
plt.yticks([])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.savefig('Ima_DelG_trace.png', dpi=500, bbox_inches='tight')





