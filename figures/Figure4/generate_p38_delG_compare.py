#p38
#Analysis conducted on lilac, here results are taken from output json files
#quickmodel --inputs inputs_multiple_well_single_wv_p38_4ligs --type singlet --nsamples 100000
# Just to note: 'P_error'       :  0.35,

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

import matplotlib.lines as mlines

# p38  "DeltaG_cred_int": "$\\Delta G$ =  -14 [-15.5,-12.6] $k_B T$",
Bos = [-15.5,-14 ,-12.6] 

# "DeltaG_cred_int": "$\\Delta G$ =  -11.9 [-13,-10.4] $k_B T$",
Bsi = [-13,-11.9 ,-10.4] 

# "DeltaG_cred_int": "$\\Delta G$ =  -10.8 [-12,-9.48] $k_B T$",
Erl = [-12,-10.8 ,-9.48] 

#  "DeltaG_cred_int": "$\\Delta G$ =  -12.9 [-14.2,-11.5] $k_B T$",
Gef = [-14.2,-12.9 ,-11.5]


labels = ['Bos','Bsi','Erl','Gef']

fig, ax = plt.subplots(figsize=(8,4))

plt.errorbar(1,Bos[1],yerr=[[abs(Bos[1]-Bos[2])],[abs(Bos[1]-Bos[0])]], fmt='o',color='c',markersize=12)
plt.errorbar(2,Bsi[1],yerr=[[abs(Bsi[1]-Bsi[2])],[abs(Bsi[1]-Bsi[0])]], fmt='o',color='b',markersize=12)
plt.errorbar(3,Erl[1],yerr=[[abs(Erl[1]-Erl[2])],[abs(Erl[1]-Erl[0])]], fmt='o',color='g',markersize=12)
plt.errorbar(4,Gef[1],yerr=[[abs(Gef[1]-Gef[2])],[abs(Gef[1]-Gef[0])]], fmt='o',color='purple',markersize=12)
plt.xticks([1,2,3,4],labels)
plt.ylim([-8,-37])
plt.ylabel('$\Delta G$ ($k_B T$)', fontsize=16)
plt.yticks(fontsize=12)
plt.xticks(fontsize=16)

plt.title('p38',fontsize=20)

plt.tight_layout()

plt.savefig('p38_delG_compare.png', dpi=500)
plt.savefig('p38_delG_compare.pdf')

