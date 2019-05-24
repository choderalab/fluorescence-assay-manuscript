import json
import numpy as np
from glob import glob

inputs = {
    'xml_file_path' :  "../../data/spectra/",
    'file_set'      :  {'Src': glob("../../data/spectra/Src/2015-12-15/*.xml"),
                        'p38': glob("../../data/spectra/p38/2016-01-26/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml")},
    'ligand_order'  :  [None,'Bosutinib Isomer', None, None],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([1.9999780000000003e-05,9.146000431435104e-06,4.182512202224779e-06,1.9126839598250786e-06,8.746800375683717e-07,3.9999559999999996e-07,1.8292000862870203e-07,8.365024404449556e-08,3.825367919650156e-08,1.749360075136743e-08,7.999911999999997e-09, 0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
