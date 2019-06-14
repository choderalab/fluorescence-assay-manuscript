import json
import numpy as np
from glob import glob

inputs = {
    'xml_file_path' :  "../../data/spectra/",
    'file_set'      :  {'Src': glob("../../data/spectra/Src/2015-12-15/*.xml"),
                        'p38': glob("../../data/spectra/p38/2016-01-26/*.xml"),
                        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml")},
    'ligand_order'  :  ['Bosutinib', None, None, None],
    'section'       :  'em280',
    'wavelength'    :  '480',
    'Lstated'       :  np.array([2.00736e-05,9.179758690368378e-06,4.197950024579236e-06,1.9197437539785283e-06,8.779085171003113e-07,4.014719999999999e-07,1.835951738073675e-07,8.395900049158471e-08,3.839487507957056e-08,1.7558170342006224e-08,8.029439999999995e-09, 0.0], np.float64), # ligand concentration
    'Pstated'       :  1.0e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
