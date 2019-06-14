import json
import numpy as np
from glob import glob


ligand_conc = [ 2.00000000e-05, 9.14610104e-06, 4.18255821e-06, 1.91270500e-06, 8.74689659e-07,
                4.00000000e-07, 1.82922021e-07, 8.36511642e-08, 3.82541000e-08, 1.74937932e-08, 8.00000000e-09, 0.00000000e+00]

inputs = {
    'xml_file_path'     :  "infinite_results",
    'file_set'      :  {'p38': glob("infinite_results/*.xml")},
    'ligand_order'  :  [None,'Bosutinib','Bosutinib Isomer','Erlotinib'],
    'section'       :  'ex280_scan_top_gain100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  0.5e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
