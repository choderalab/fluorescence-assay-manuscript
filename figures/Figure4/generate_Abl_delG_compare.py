#Abl WT and GK
#Analysis conducted on lilac, here results are taken from output json files
#quickmodel --inputs inputs_multiple_well_single_wv_Abl_4ligs --type singlet --nsamples 100000
# Just to note: 'P_error'       :  0.35,

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

import matplotlib.lines as mlines

# Abl "DeltaG_cred_int": "$\\Delta G$ =  -27.4 [-34.2,-20.4] $k_B T$",
Bos = [-34.2,-27.4,-20.4] 

# AblGK "DeltaG_cred_int": "$\\Delta G$ =  -26.9 [-34.2,-19.3] $k_B T$",
GK_Bos = [-34.2,-26.9,-19.3]

# "DeltaG_cred_int": "$\\Delta G$ =  -26.6 [-34.1,-18.9] $k_B T$",
Bsi = [-34.1,-26.6,-18.9] 

# "DeltaG_cred_int": "$\\Delta G$ =  -26.6 [-34.1,-18.6] $k_B T$",
GK_Bsi =[-34.1,-26.6 ,-18.6]

# "DeltaG_cred_int": "$\\Delta G$ =  -22 [-33.9,-15.4] $k_B T$",
Erl = [-33.9,-22 ,-15.4] 

# "DeltaG_cred_int": "$\\Delta G$ =  -24.7 [-34,-17.3] $k_B T$",
GK_Erl = [-34,-24.7 ,-17.3]

# "DeltaG_cred_int": "$\\Delta G$ =  -14.1 [-15.4,-12.9] $k_B T$"
Gef = [-15.4,-14.1 ,-12.9]

# "DeltaG_cred_int": "$\\Delta G$ =  -13.5 [-13.5,-13.5] $k_B T$",
GK_Gef = [-13.5,-13.5,-13.5] 

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

plt.title('Abl',fontsize=20)

WT = mlines.Line2D([0], [0], marker='o', color='k', label='WT',
                          markerfacecolor='k', markersize=12)
GK = mlines.Line2D([0], [0], marker='X', color='k', label='GK',
                          markerfacecolor='k', markersize=12)
plt.legend(handles=[WT,GK],fontsize=12,loc=0,frameon=True)

plt.tight_layout()

plt.savefig('Abl_delG_compare.png', dpi=500)
plt.savefig('Abl_delG_compare.pdf')

