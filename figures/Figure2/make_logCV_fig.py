import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize
import seaborn as sns

sns.set(style='white')
sns.set_context('talk')

import pickle

Kd_list = [
           0.02e-9, # M
           0.2e-9, # 
           2e-9, # M
           20e-9, # M
           200e-9, # M
           2000e-9, # M
           20000e-9, # M
          ]
dark_colors = ['cyan','blue','red','green','yellowgreen','goldenrod','yellow']
          
trace_list = []
for i in range(len(Kd_list)):
    data_file = 'Kdlist_%s_mcmc.pickle'%i
    with open(r'%s'%data_file,'rb') as my_file:
        data = pickle.load(my_file)
    trace_list.append(data['DeltaG'][0])

#Now let's calculate the CV and relative bias of the delG instead results

CV_list = []
bias_list = []

for i,trace in enumerate(trace_list):
    CV_list.append(np.std(trace)/np.mean(trace) * 100)
    bias_list.append((np.mean(trace) - np.log(Kd_list[i])) / np.log(Kd_list[i]))

f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(8,4.5))

ax1.axhline(y=0,color='0.8',linestyle='--')
for i,CV in enumerate(CV_list):
    ax1.semilogx(Kd_list[i], CV, dark_colors[i], marker = 'o', linestyle='None')
    
ax1.set_ylabel('%CV ($log_{10}(K_{d})$)',fontsize=14);
#ax1.set_ylim((-15,0));

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax2.axhline(y=0,color='0.8',linestyle='--')
for i,bias in enumerate(bias_list):
    ax2.semilogx(Kd_list[i], bias*100, dark_colors[i], marker = 'o', linestyle='None')

ax2.set_ylabel('%RB ($log_{10}(K_{d})$)',fontsize=14);
#ax2.set_ylim((-5,30));

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

f.subplots_adjust(hspace=0)
plt.xlabel('$K_{d}$ ($M$)',fontsize=20);
plt.xticks(fontsize=16);
plt.yticks(fontsize=16);

plt.tight_layout()

plt.savefig('simulated_fluorescence_log_CV_bias.png', dpi=500, bbox_inches='tight')
plt.savefig('simulated_fluorescence_log_CV_bias.pdf', bbox_inches='tight')

plt.close()