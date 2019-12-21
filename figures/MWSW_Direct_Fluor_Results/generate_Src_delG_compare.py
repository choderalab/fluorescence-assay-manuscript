#Src WT and GK
#Analysis conducted on lilac, here results are taken from output json files
#quickmodel --inputs inputs_multiple_well_single_wv_Src_4ligs --type singlet --nsamples 100000
# Just to note: 'P_error'       :  0.15,

import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

import matplotlib.lines as mlines

# Import 95% credibility intervals from json file quickmodel results
import json
Bos_file = 'mwsw/Src-Bosutinib-AB-2019-10-04 09:04.json'
GK_Bos_file = 'mwsw/Src_T338I-Bosutinib-AB-2019-10-04 09:05.json'
Bsi_file = 'mwsw/Src-Bosutinib Isomer-CD-2019-10-04 14:42.json'
GK_Bsi_file = 'mwsw/Src_T338I-Bosutinib Isomer-CD-2019-10-04 15:09.json'
Erl_file = 'mwsw/Src-Erlotinib-EF-2019-10-04 20:22.json'
GK_Erl_file = 'mwsw/Src_T338I-Erlotinib-EF-2019-10-04 21:00.json'
Gef_file = 'mwsw/Src-Gefitinib-GH-2019-10-05 01:56.json'
GK_Gef_file = 'mwsw/Src_T338I-Gefitinib-GH-2019-10-05 02:22.json'

Bos_spectra_file = '../../analysis/bayes/corrected-concentrations-Perror15/Bos/Src-Bosutinib-AB-2019-10-07 11:37.json'
Bsi_spectra_file = '../../analysis/bayes/corrected-concentrations-Perror15/Bsi/Src-Bosutinib Isomer-CD-2019-10-07 11:38.json'
Erl_spectra_file = '../../analysis/bayes/corrected-concentrations-Perror15/Erl/Src-Erlotinib-EF-2019-10-07 11:40.json'
Gef_spectra_file = '../../analysis/bayes/corrected-concentrations-Perror15/Gef/Src-Gefitinib-GH-2019-10-07 11:40.json'


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
GK_Bos = json_parse(GK_Bos_file)
Bsi = json_parse(Bsi_file)
GK_Bsi = json_parse(GK_Bsi_file)
Erl = json_parse(Erl_file)
GK_Erl = json_parse(GK_Erl_file)
Gef = json_parse(Gef_file)
GK_Gef = json_parse(GK_Gef_file)

Bos_spectra = json_parse(Bos_spectra_file)
Bsi_spectra = json_parse(Bsi_spectra_file)
Erl_spectra = json_parse(Erl_spectra_file)
Gef_spectra = json_parse(Gef_spectra_file)

labels = ['bosutinib','bosutinib isomer','erlotinib','gefitinib']

fig, ax = plt.subplots(figsize=(6,8))

plt.errorbar(Bos_spectra[1],0.75,xerr=[[abs(Bos_spectra[1]-Bos_spectra[2])],[abs(Bos_spectra[1]-Bos_spectra[0])]], fmt='.',color='c',markersize=12)
plt.errorbar(Bos[1],1,xerr=[[abs(Bos[1]-Bos[2])],[abs(Bos[1]-Bos[0])]], fmt='o',color='c',markersize=12)
plt.errorbar(GK_Bos[1],1.25,xerr=[[abs(GK_Bos[1]-GK_Bos[2])],[abs(GK_Bos[1]-GK_Bos[0])]], fmt='X',color='c',markersize=12)

plt.errorbar(Bsi_spectra[1],1.75,xerr=[[abs(Bsi_spectra[1]-Bsi_spectra[2])],[abs(Bsi_spectra[1]-Bsi_spectra[0])]], fmt='.',color='b',markersize=12)
plt.errorbar(Bsi[1],2,xerr=[[abs(Bsi[1]-Bsi[2])],[abs(Bsi[1]-Bsi[0])]], fmt='o',color='b',markersize=12)
plt.errorbar(GK_Bsi[1],2.25,xerr=[[abs(GK_Bsi[1]-GK_Bsi[2])],[abs(GK_Bsi[1]-GK_Bsi[0])]], fmt='X',color='b',markersize=12)

plt.errorbar(Erl_spectra[1],2.75,xerr=[[abs(Erl_spectra[1]-Erl_spectra[2])],[abs(Erl_spectra[1]-Erl_spectra[0])]], fmt='.',color='g',markersize=12)
plt.errorbar(Erl[1],3,xerr=[[abs(Erl[1]-Erl[2])],[abs(Erl[1]-Erl[0])]], fmt='o',color='g',markersize=12)
plt.errorbar(GK_Erl[1],3.25,xerr=[[abs(GK_Erl[1]-GK_Erl[2])],[abs(GK_Erl[1]-GK_Erl[0])]], fmt='X',color='g',markersize=12)

plt.errorbar(Gef_spectra[1],3.75,xerr=[[abs(Gef_spectra[1]-Gef_spectra[2])],[abs(Gef_spectra[1]-Gef_spectra[0])]], fmt='.',color='purple',markersize=12)
plt.errorbar(Gef[1],4,xerr=[[abs(Gef[1]-Gef[2])],[abs(Gef[1]-Gef[0])]], fmt='o',color='purple',markersize=12)
plt.errorbar(GK_Gef[1],4.25,xerr=[[abs(GK_Gef[1]-GK_Gef[2])],[abs(GK_Gef[1]-GK_Gef[0])]], fmt='X',color='purple',markersize=12)

plt.xlim([-8,-37])
plt.yticks([1,2,3,4],labels)
plt.xlabel('$\Delta G$ ($k_B T$)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=20)

plt.title('Src',fontsize=20)

Spectra = mlines.Line2D([0], [0], marker='.', color='k', label='spectra',
                          markerfacecolor='k', markersize=12)
WT = mlines.Line2D([0], [0], marker='o', color='k', label='WT',
                          markerfacecolor='k', markersize=12)
GK = mlines.Line2D([0], [0], marker='X', color='k', label='GK',
                          markerfacecolor='k', markersize=12)
plt.legend(handles=[Spectra,WT,GK],fontsize=12,loc=0,frameon=True)

plt.tight_layout()

plt.savefig('Src_delG_compare.png', dpi=500)
plt.savefig('Src_delG_compare.pdf')
