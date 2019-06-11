import json
import numpy as np
from glob import glob

xml_files = ['../infinite_results/p38_Abl_WT_GK_Src_WT_GK_conc_0_20190307_111234.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GK_conc_1_20190307_112351.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_2_20190307_113345.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_3_20190307_114336.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_4_20190307_115329.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_5_20190307_120322.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_6_20190307_121315.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_7_20190307_122657.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_8_20190307_123649.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_9_20190307_124642.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_10_20190307_125635.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_11_20190307_130625.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_12_20190307_131618.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_13_20190307_132611.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_14_20190307_133612.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_15_20190307_134607.xml',
             '../infinite_results/p38_Abl_WT_GK_Src_WT_GKconc_16_20190307_135607.xml']

ligand_conc = [  0.00000000e+00,   8.00000000e-09,   1.34778097e-08,
         2.27064194e-08,   3.82541000e-08,   6.44476851e-08,
         1.08576705e-07,   1.82922021e-07,   3.08173524e-07,
         5.19188015e-07,   8.74689659e-07,   1.47361260e-06, 2.48263378e-06,
         4.18255821e-06, 7.04646547e-06, 1.118713651e-05, 2.0e-05]

inputs = {
    'single_well'   :  True,
    'xml_files'     :  xml_files,
    'file_set'      :  {'Src_WT': xml_files},
    'protein_wells'  :  {'Src_WT': ['B7', 'D7', 'F7', 'H7']},
    'ligand_order'  :  ['Bosutinib', 'Bosutinib Isomer', 'Erlotinib', 'Gefitinib'],
    'buffer_wells'   :  {'Src_WT': ['A2', 'C2', 'E2', 'G2']},
    'section'       :  'ex280_em480_top_gain100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  0.5e-6 * np.ones([17],np.float64), # protein concentration, M
    'P_error'       :  0.15, # coefficient of protein concentration uncertainity (0.15 for 15% error)
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
