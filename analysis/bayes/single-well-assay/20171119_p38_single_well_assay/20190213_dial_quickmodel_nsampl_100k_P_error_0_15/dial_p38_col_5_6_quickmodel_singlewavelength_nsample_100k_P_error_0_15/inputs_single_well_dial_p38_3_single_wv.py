import json
import numpy as np
from glob import glob

xml_files = ['../infinite_results/2017-11-20 15-48-06_plate_1.xml',
            '../infinite_results/2017-11-20 16-17-59_plate_1.xml',
            '../infinite_results/2017-11-20 16-41-35_plate_1.xml',
            '../infinite_results/2017-11-20 16-59-09_plate_1.xml',
            '../infinite_results/2017-11-20 17-19-31_plate_1.xml',
            '../infinite_results/2017-11-20 17-37-37_plate_1.xml',
            '../infinite_results/2017-11-20 17-58-04_plate_1.xml',
            '../infinite_results/2017-11-20 18-16-02_plate_1.xml',
            '../infinite_results/2017-11-20 18-35-42_plate_1.xml',
            '../infinite_results/2017-11-20 18-54-33_plate_1.xml',
            '../infinite_results/2017-11-20 19-12-49_plate_1.xml',
            '../infinite_results/2017-11-20 19-31-32_plate_1.xml']

ligand_conc = [ 0.00000000e+00, 8.00000000e-09, 1.74937932e-08, 3.82541000e-08,
                           8.36511642e-08, 1.82922021e-07, 4.00000000e-07, 8.74689659e-07,
                           1.91270500e-06, 4.18255821e-06, 9.14610104e-06, 2.00000000e-05 ]

inputs = {
    'single_well'   :  True,
    'xml_files'     :  xml_files,
    'file_set'      :  {'dialyzed_p38_3': xml_files},
    'protein_wells'  :  {'dialyzed_p38_3': ['A5','B5','C5','D5','E5','F5','G5','H5']},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Gefitinib','Erlotinib','Ponatinib','Lapatinib','Pazopanib','Axitinib'],
    'buffer_wells'   :  {'dialyzed_p38_3': ['A6','B6','C6','D6','E6','F6','G6','H6']},
    'section'       :  '280_480_TOP_100',
    'wavelength'    :  '480',
    'Lstated'       :  np.array(ligand_conc, np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'P_error'       :  0.15, # coefficient of protein concentration uncertainty (0.15 for 15% error)
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
