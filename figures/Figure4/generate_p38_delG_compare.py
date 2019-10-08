#p38
#Analysis conducted on lilac, here results are taken from output json files
#quickmodel --inputs inputs_multiple_well_single_wv_p38_4ligs --type singlet --nsamples 100000
# Just to note: 'P_error'       :  0.15,

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

import matplotlib.lines as mlines

# Import 95% credibility intervals from json file quickmodel results
import json
Bos_file = 'mwsw/p38-Bosutinib-AB-2019-10-04 09:06.json'
Bsi_file = 'mwsw/p38-Bosutinib Isomer-CD-2019-10-04 12:13.json'
Erl_file = 'mwsw/p38-Erlotinib-EF-2019-10-04 15:18.json'
Gef_file = 'mwsw/p38-Gefitinib-GH-2019-10-04 18:51.json'

def json_parse(file):
    with open(file) as json_file:
        data = json.load(json_file)
    three = data['DeltaG_cred_int'].split("=")[1].split(" $")[0].strip().split("[")
    first = float(three[0].strip())
    inner = three[1].split(",")
    second = float(inner[0])
    third = float(inner[1][:-1])
    lig = [third,
           first,
           second]

    print(lig)
    return(lig)

Bos = json_parse(Bos_file)
Bsi = json_parse(Bsi_file)
Erl = json_parse(Erl_file)
Gef = json_parse(Gef_file)

labels = ['bosutinib','bosutinib isomer','erlotinib','gefitinib']


fig, ax = plt.subplots(figsize=(6,4))

plt.errorbar(Bos[1],1,xerr=[[abs(Bos[1]-Bos[2])],[abs(Bos[1]-Bos[0])]], fmt='o',color='c',markersize=12)
plt.errorbar(Bsi[1],2,xerr=[[abs(Bsi[1]-Bsi[2])],[abs(Bsi[1]-Bsi[0])]], fmt='o',color='b',markersize=12)
plt.errorbar(Erl[1],3,xerr=[[abs(Erl[1]-Erl[2])],[abs(Erl[1]-Erl[0])]], fmt='o',color='g',markersize=12)
plt.errorbar(Gef[1],4,xerr=[[abs(Gef[1]-Gef[2])],[abs(Gef[1]-Gef[0])]], fmt='o',color='purple',markersize=12)

plt.xlim([-8,-37])
plt.yticks([1,2,3,4],labels)
plt.xlabel('$\Delta G$ ($k_B T$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=20)

plt.title('p38',fontsize=20)

plt.tight_layout()

plt.savefig('p38_delG_compare.png', dpi=500)
plt.savefig('p38_delG_compare.pdf')

