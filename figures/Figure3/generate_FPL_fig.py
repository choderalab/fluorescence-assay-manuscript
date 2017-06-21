import cPickle
import numpy as np
import matplotlib.pyplot as plt

#Note that these data are from dP_stated = 0.15 analyses.

Src_Bos_file = '../Figure4/dPstated/Src-Bosutinib-AB_mcmc-0.pickle'
Src_Bsi_file = '../Figure4/dPstated/Src-Bosutinib Isomer-CD_mcmc-0.pickle'
Src_Erl_file = '../Figure4/dPstated/Src-Erlotinib-EF_mcmc-0.pickle'
Src_Gef_file = '../Figure4/dPstated/Src-Gefitinib-GH_mcmc-0.pickle'

with open(r'%s'%Src_Bos_file,'rb') as my_file:
    Src_Bos_data = cPickle.load(my_file)
with open(r'%s'%Src_Bsi_file,'rb') as my_file:
    Src_Bsi_data = cPickle.load(my_file)
with open(r'%s'%Src_Erl_file,'rb') as my_file:
    Src_Erl_data = cPickle.load(my_file)
with open(r'%s'%Src_Gef_file,'rb') as my_file:
    Src_Gef_data = cPickle.load(my_file)

#Next we need to define 95% credibility intervals
FPL_interval_Src_Bos = np.percentile(a=Src_Bos_data['F_PL'][0], q=[2.5, 50.0, 97.5])
FPL_interval_Src_Bsi = np.percentile(a=Src_Bsi_data['F_PL'][0], q=[2.5, 50.0, 97.5])
FPL_interval_Src_Erl = np.percentile(a=Src_Erl_data['F_PL'][0], q=[2.5, 50.0, 97.5])
FPL_interval_Src_Gef = np.percentile(a=Src_Gef_data['F_PL'][0], q=[2.5, 50.0, 97.5])

FL_interval_Src_Bos = np.percentile(a=Src_Bos_data['F_L'][0], q=[2.5, 50.0, 97.5])
FL_interval_Src_Bsi = np.percentile(a=Src_Bsi_data['F_L'][0], q=[2.5, 50.0, 97.5])
FL_interval_Src_Erl = np.percentile(a=Src_Erl_data['F_L'][0], q=[2.5, 50.0, 97.5])
FL_interval_Src_Gef = np.percentile(a=Src_Gef_data['F_L'][0], q=[2.5, 50.0, 97.5])

#Next we need to write this in a form so we can take it into our bar plot fun
F_PL_Src = [FPL_interval_Src_Bos[1],FPL_interval_Src_Bsi[1],FPL_interval_Src_Erl[1],FPL_interval_Src_Gef[1]]
F_PL_Src_err = [[FPL_interval_Src_Bos[1]-FPL_interval_Src_Bos[0],
                 FPL_interval_Src_Bsi[1]-FPL_interval_Src_Bsi[0],
                 FPL_interval_Src_Erl[1]-FPL_interval_Src_Erl[0],
                 FPL_interval_Src_Gef[1]-FPL_interval_Src_Gef[0]],
                [FPL_interval_Src_Bos[2]-FPL_interval_Src_Bos[1],
                 FPL_interval_Src_Bsi[2]-FPL_interval_Src_Bsi[1],
                 FPL_interval_Src_Erl[2]-FPL_interval_Src_Erl[1],
                 FPL_interval_Src_Gef[2]-FPL_interval_Src_Gef[1]]]

F_L_Src = [FL_interval_Src_Bos[1],FL_interval_Src_Bsi[1],FL_interval_Src_Erl[1],FL_interval_Src_Gef[1]]
F_L_Src_err = [[FL_interval_Src_Bos[1]-FL_interval_Src_Bos[0],
                FL_interval_Src_Bsi[1]-FL_interval_Src_Bsi[0],
                FL_interval_Src_Erl[1]-FL_interval_Src_Erl[0],
                FL_interval_Src_Gef[1]-FL_interval_Src_Gef[0]],
               [FL_interval_Src_Bos[2]-FL_interval_Src_Bos[1],
                FL_interval_Src_Bsi[2]-FL_interval_Src_Bsi[1],
                FL_interval_Src_Erl[2]-FL_interval_Src_Erl[1],
                FL_interval_Src_Gef[2]-FL_interval_Src_Gef[1]]]


# define coords for plotting
x_FPL = range(len(F_PL_Src))
x_FL = [x + 0.2 for x in x_FPL]

labels = 'Src:Bos','Src:Bsi','Src:Erl','Src:Gef'

fig, ax = plt.subplots(figsize=(6,4))

plt.bar(x_FPL,F_PL_Src,width=0.2,yerr=F_PL_Src_err,align='center',color='C0',label='F_PL',log=True)
plt.bar(x_FL,F_L_Src,width=0.2,yerr=F_L_Src_err,align='center',color='C1',label='F_L',log=True)
plt.yticks(fontsize=16)
plt.ylabel('$F$', fontsize=20)
plt.xticks(range(len(F_PL_Src)),labels, rotation=45, fontsize=16);
plt.legend(fontsize=14);

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

plt.savefig('Src_FPL_FL_barplot.png', dpi=500)
plt.savefig('Src_FPL_FL_barplot.pdf')


