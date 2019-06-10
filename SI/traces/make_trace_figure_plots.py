import matplotlib.pyplot as plt
import pickle
import matplotlib
import numpy as np

Src_Erl_file = '../../analysis/bayes/corrected-concentrations/erl/Src-Erlotinib-EF_mcmc-2019-05-17 13:48.pickle'

with open(r'%s'%Src_Erl_file,'rb') as my_file:
    Src_Erl_data = pickle.load(my_file)
	
priority_keys = ['Ptrue','Ltrue','F_buffer','F_plate','F_P','F_L','F_PL','DeltaG','epsilon_ex','epsilon_em','sigma_top','sigma_abs']

colors = matplotlib.cm.viridis(np.linspace(0, 0.85, len(priority_keys)))

fig, axs = plt.subplots(int(len(priority_keys)/2),2, figsize=(13, 20), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, key in enumerate(priority_keys):
    try:
        axs[i].plot(Src_Erl_data[key][0],color=colors[i],marker='.',linestyle = 'None');
        axs[i].set_title('%s'%key,fontsize=14);
        axs[i].set_xlabel('MCMC sample')
        axs[i].set_ylabel('%s'%key)
        plt.tight_layout()
    except Exception as e:
        print('%s has no [0]'%key)
plt.savefig('priority_keys_traces.png',dpi=300)

plt.clf()

fig, axs = plt.subplots(int(len(priority_keys)/2),2, figsize=(13, 20), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, key in enumerate(priority_keys):
    try:
        axs[i].plot(Src_Erl_data['DeltaG'][0],Src_Erl_data[key][0],color=colors[i],marker='.',linestyle = 'None');
        axs[i].set_title('DeltaG v. %s'%key,fontsize=14);
        axs[i].set_xlabel('DeltaG')
        axs[i].set_ylabel('%s'%key)
        plt.tight_layout()
    except Exception as e:
        print('%s has no [0]'%key)
plt.savefig('priority_keys_traces_vs_DeltaG.png',dpi=300)

