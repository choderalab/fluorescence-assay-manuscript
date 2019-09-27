#Src WT and GK
#Analysis conducted on lilac, here results are taken from output json files
#quickmodel --inputs inputs_multiple_well_single_wv_Src_4ligs --type singlet --nsamples 100000
# Just to note: 'P_error'       :  0.35,

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

import matplotlib.lines as mlines

# Src "DeltaG_cred_int": "$\\Delta G$ =  -26.6 [-34.2,-19.2] $k_B T$"
Bos = [-34.2,-26.6,-19.2]

# SrcGK "DeltaG_cred_int": "$\\Delta G$ =  -26.6 [-34.1,-19.3] $k_B T$",
GK_Bos = [-34.1,-26.6,-19.3]

# "DeltaG_cred_int": "$\\Delta G$ =  -27 [-34.2,-19.7] $k_B T$"
Bsi = [-34.2,-27,-19.7]

# "DeltaG_cred_int": "$\\Delta G$ =  -15.1 [-20.7,-13.9] $k_B T$",
GK_Bsi = [-20.7,-15.1,-13.9]

# "DeltaG_cred_int": "$\\Delta G$ =  -14 [-15.3,-12.8] $k_B T$",
Erl = [-15.3,-14,-12.8]

# "DeltaG_cred_int": "$\\Delta G$ =  -14 [-15.9,-12.6] $k_B T$",
GK_Erl = [-15.9,-14,-12.6]

# "DeltaG_cred_int": "$\\Delta G$ =  -13.5 [-14.6,-12.6] $k_B T$",
Gef = [-14.6,-13.5,-12.6]

# "DeltaG_cred_int": "$\\Delta G$ =  -12.2 [-12.7,-11.8] $k_B T$"
GK_Gef = [-12.7,-12.2,-11.8]

labels = ['Bos','Bsi','Erl','Gef']

fig, ax = plt.subplots(figsize=(8,4))

plt.errorbar(0.75,Bos[1],yerr=[[abs(Bos[1]-Bos[2])],[abs(Bos[1]-Bos[0])]], fmt='o',color='c',markersize=12)
plt.errorbar(1.25,GK_Bos[1],yerr=[[abs(GK_Bos[1]-GK_Bos[2])],[abs(GK_Bos[1]-GK_Bos[0])]], fmt='X',color='c',markersize=12)
plt.errorbar(1.75,Bsi[1],yerr=[[abs(Bsi[1]-Bsi[2])],[abs(Bsi[1]-Bsi[0])]], fmt='o',color='b',markersize=12)
plt.errorbar(2.25,GK_Bsi[1],yerr=[[abs(GK_Bsi[1]-GK_Bsi[2])],[abs(GK_Bsi[1]-GK_Bsi[0])]], fmt='X',color='b',markersize=12)
plt.errorbar(2.75,Erl[1],yerr=[[abs(Erl[1]-Erl[2])],[abs(Erl[1]-Erl[0])]], fmt='o',color='g',markersize=12)
plt.errorbar(3.25,GK_Erl[1],yerr=[[abs(GK_Erl[1]-GK_Erl[2])],[abs(GK_Erl[1]-GK_Erl[0])]], fmt='X',color='g',markersize=12)
plt.errorbar(3.75,Gef[1],yerr=[[abs(Gef[1]-Gef[2])],[abs(Gef[1]-Gef[0])]], fmt='o',color='purple',markersize=12)
plt.errorbar(4.25,GK_Gef[1],yerr=[[abs(GK_Gef[1]-GK_Gef[2])],[abs(GK_Gef[1]-GK_Gef[0])]], fmt='X',color='purple',markersize=12)

plt.ylim([-8,-37])
plt.xticks([1,2,3,4],labels)
plt.ylabel('$\Delta G$ ($k_B T$)', fontsize=16)
plt.yticks(fontsize=12)
plt.xticks(fontsize=16)

plt.title('Src',fontsize=20)

WT = mlines.Line2D([0], [0], marker='o', color='k', label='WT',
                          markerfacecolor='k', markersize=12)
GK = mlines.Line2D([0], [0], marker='X', color='k', label='GK',
                          markerfacecolor='k', markersize=12)
plt.legend(handles=[WT,GK],fontsize=12,loc=0,frameon=True)

plt.tight_layout()

plt.savefig('Src_delG_compare.png', dpi=500)
plt.savefig('Src_delG_compare.pdf')

