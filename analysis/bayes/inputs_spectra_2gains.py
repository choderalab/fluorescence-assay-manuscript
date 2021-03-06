import json
import numpy as np
from glob import glob

inputs = {
    'xml_file_path' :  "../../data/spectra/",
    'file_set'      :  {'Src': glob("../../data/spectra/Src/2016-03-09/*.xml"),
                        'p38': glob("../../data/spectra/p38/2016-03-07/*.xml"),
                        'CAII': glob("../../data/spectra/CAII/2016-03-15/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2016-03-11/*.xml")},
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib'],
    'section'       :  'em280_Copy2',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
