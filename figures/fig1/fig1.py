from assaytools import grant
from glob import glob
import matplotlib.pyplot as plt

file_set = {'Src': glob("../../data/spectra/Src/2015-12-15/*.xml"),
        'Abl': glob("../../data/spectra/Abl/2015-12-18/*.xml"),
        'p38': glob("../../data/spectra/p38/2016-01-26/*.xml")}
ligands = ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib']

grant.plot_spectra_grid(file_set,ligands,output='Figure1')