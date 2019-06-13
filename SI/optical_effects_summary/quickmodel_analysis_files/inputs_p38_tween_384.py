import json
import numpy as np
from glob import glob


ligand_conc = [2e-05, 1.40145183674e-05, 9.82033625345e-06, 6.88136413989e-06, 4.82195020655e-06, 3.37886548681e-06, 2.36765862128e-06, 1.65907976178e-06, 1.16256018972e-06, 8.146360566e-07, 5.70836598897e-07, 4e-07, 2.80290367347e-07, 1.96406725069e-07, 1.37627282798e-07, 9.64390041309e-08, 6.75773097361e-08, 4.73531724257e-08, 3.31815952356e-08, 2.32512037944e-08, 1.6292721132e-08, 1.14167319779e-08, 8e-09, 0.00000000e+00,]

xml_files = ["infinite_results/p38_4lig_tween_20170627_150509.xml"]

inputs = {
    'xml_file_path'     :  "infinite_results/",
    'file_set'      :  {'p38' : xml_files},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  '280_480_TOP_100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  0.5e-6 * np.ones([24],np.float64), # protein concentration, M
    'P_error'       :  0.15, # coefficient of protein concentration uncertainity (0.15 for 15% error)
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
