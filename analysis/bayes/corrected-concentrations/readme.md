### Ligand concentration arrays corrected to actual dispensed vial concentrations.

#### Table of ligand information. 
Vials were made using a Quantos QB5 automated liquid and powder dispensing system.

| Ligand | Supplier | CAT No. | Lot No. | Conc. (mg/g) | Conc. (mg/mL) | Conc. (M) |
| ------ | -------- | ------- | ------- | ------------ | ------------- | --------- |
| Bosutinib | LC Labs | B-1788 | BSB-102 | 4.84 | 5.324 | 0.0100368 | 
| Bosutinib Isomer | LC Labs | B-1722 | BWC-102 | 4.82  | 5.304 | 0.00999989 |
| Erlotinib | LC Labs | E-4007 | BBE-107 | 3.59 | 3.953 | 0.00919288 |
| Gefitinib | LC Labs | G-4408 | BGF-106 | 4.07 | 4.477 | 0.0100179 |

#### Inputs into _calculate_Lstated_array.py_ from Assaytools.

| Ligand | Script Input | Script Output |
| ------ | ------------ | ------------- |
| Bosutinib | python calculate_Lstated_array.py --n_wells 11 --h_conc 20e-6 --l_conc 8e-9 --target_stock_conc 0.010 --true_stock_conc 0.0100368 --dilution logarithmic | np.array([2.00736e-05,9.179758690368378e-06,4.197950024579236e-06,1.9197437539785283e-06,8.779085171003113e-07,4.014719999999999e-07,1.835951738073675e-07,8.395900049158471e-08,3.839487507957056e-08,1.7558170342006224e-08,8.029439999999995e-09, 0.0], np.float64) |
| Bosutinib Isomer | python calculate_Lstated_array.py --n_wells 11 --h_conc 20e-6 --l_conc 8e-9 --target_stock_conc 0.010 --true_stock_conc 0.00999989 --dilution logarithmic | np.array([1.9999780000000003e-05,9.146000431435104e-06,4.182512202224779e-06,1.9126839598250786e-06,8.746800375683717e-07,3.9999559999999996e-07,1.8292000862870203e-07,8.365024404449556e-08,3.825367919650156e-08,1.749360075136743e-08,7.999911999999997e-09, 0.0], np.float64)
| Erlotinib | python calculate_Lstated_array.py --n_wells 11 --h_conc 20e-6 --l_conc 8e-9 --target_stock_conc 0.010 --true_stock_conc 0.00919288 --dilution logarithmic | np.array([1.838576e-05,8.407900931523358e-06,3.844975572090105e-06,1.7583267536539667e-06,8.040917073849344e-07,3.677151999999999e-07,1.6815801863046714e-07,7.689951144180208e-08,3.516653507307932e-08,1.6081834147698685e-08,7.354303999999996e-09, 0.0], np.float64)
| Gefitinib | python calculate_Lstated_array.py --n_wells 11 --h_conc 20e-6 --l_conc 8e-9 --target_stock_conc 0.010 --true_stock_conc 0.0100179 --dilution logarithmic | np.array([2.0035800000000002e-05,9.162472559405526e-06,4.1900449895616465e-06,1.916128741529322e-06,8.762553536445092e-07,4.0071599999999995e-07,1.8324945118811045e-07,8.38008997912329e-08,3.8322574830586434e-08,1.752510707289018e-08,8.014319999999996e-09, 0.0], np.float64)

Note: A 0.0 concentration was manually appended to the end of the concentration arrays to mirror the actual experiment in which no ligand was dispensed into one of the wells. 
