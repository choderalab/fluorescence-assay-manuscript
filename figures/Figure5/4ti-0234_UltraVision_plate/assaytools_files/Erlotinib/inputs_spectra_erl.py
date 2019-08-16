import json
import numpy as np
from glob import glob

inputs = {
    'xml_file_path' :  "../../data/spectra/",
    'file_set'      :  {'Src': glob("../../data/spectra/Src/2015-12-15/*.xml"),
                        'p38': glob("../../data/spectra/p38/2016-01-26/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml")},
    'ligand_order'  :  [None, None, 'Erlotinib', None],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([1.838576e-05,8.407900931523358e-06,3.844975572090105e-06,1.7583267536539667e-06,8.040917073849344e-07,3.677151999999999e-07,1.6815801863046714e-07,7.689951144180208e-08,3.516653507307932e-08,1.6081834147698685e-08,7.354303999999996e-09, 0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
