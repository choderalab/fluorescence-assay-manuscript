
import numpy as np
from assaytools import parser
import string
from glob import glob

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches
import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

ligand_conc = [ 2.00000000e-05, 9.14610104e-06, 4.18255821e-06, 1.91270500e-06, 8.74689659e-07,
                4.00000000e-07, 1.82922021e-07, 8.36511642e-08, 3.82541000e-08, 1.74937932e-08, 8.00000000e-09, 0.00000000e+00]

inputs = {
    'xml_file_path' :  "../../data/single-wavelength/",
    'file_set'      :  {'Src': glob("../../data/single-wavelength/Src/2019-07-10/*.xml"),                        
                        'p38': glob("../../data/single-wavelength/p38/2019-07-09/*.xml"),
                        'Abl': glob("../../data/single-wavelength/Abl/2019-07-10/*.xml")},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  '280_480_TOP_100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  0.5e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

[complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)

cols = sns.color_palette('PuBu_r', 5)

fig, ax = plt.subplots(figsize=(8,4))

plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib-AB']/complex_fluorescence['p38-Bosutinib-AB'].max(),marker='o',color=cols[0],label='Bosutinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib-AB']/complex_fluorescence['p38-Bosutinib-AB'].max(),linestyle='--',color=cols[0])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Bosutinib Isomer-CD']/complex_fluorescence['p38-Bosutinib Isomer-CD'].max(),marker='o',color=cols[1],label='Bosutinib Isomer')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Bosutinib Isomer-CD']/complex_fluorescence['p38-Bosutinib Isomer-CD'].max(),linestyle='--',color=cols[1])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Erlotinib-EF']/complex_fluorescence['p38-Erlotinib-EF'].max(),marker='o',color=cols[2],label='Erlotinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Erlotinib-EF']/complex_fluorescence['p38-Erlotinib-EF'].max(),linestyle='--',color=cols[2])
plt.semilogx(inputs['Lstated'],complex_fluorescence['p38-Gefitinib-GH']/complex_fluorescence['p38-Gefitinib-GH'].max(),marker='o',color=cols[3],label='Gefitinib')
plt.semilogx(inputs['Lstated'],ligand_fluorescence['p38-Gefitinib-GH']/complex_fluorescence['p38-Gefitinib-GH'].max(),linestyle='--',color=cols[3])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title('p38',fontsize=20)
plt.yticks([])
plt.ylabel('Normalized Fluorescence',fontsize=16)
plt.xticks(fontsize=16)
plt.xlabel('$[L] (M)$',fontsize=16)
plt.legend(loc=2)

plt.tight_layout()

plt.savefig('p38_binding_curve.png', dpi=500)
plt.savefig('p38_binding_curve.pdf')

