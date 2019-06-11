# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 19:51:46 2019

@author: Vall
"""

import os
import iv_save_module as ivs

rods_filename = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\Análisis\Rods_LIGO1.txt'
folder_base = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\Mediciones\{}\Ajustes'

filenames = []
rods = []
with open(rods_filename, 'r') as file:
    for line in file:
        if line[0]!='#':
            filenames.append(line.split('\t')[0])
            rods.append(line.split('\t')[1].split('\n')[0])

folders = []
for file in filenames:
    day = file.split('_')[1]
    day = day[:4] + '-' + day[4:6] + '-' + day[6:]
    folders.append(folder_base.format(day))
del day

fits_data = []
fits_footer = []
for fil, fol in zip(filenames, folders):
    data, fits_header, footer = ivs.loadTxt(os.path.join(fol, fil+'.txt'))
    fits_data.append(data)
    fits_footer.append(footer)
del fil, fol