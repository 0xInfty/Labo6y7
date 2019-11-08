# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import numpy as np
import iv_save_module as ivs
import iv_analysis_module as iva
import matplotlib.pyplot as plt

#%%

# Parameters
density = 19300 # kg/m3
file = r'C:\Users\Valeria\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1_PostUSA.txt'

# Load data
data, header, footer = ivs.loadTxt(file)
length = data[:, 2] * 1e-9 # from nm to m
damping_time = data[:, 7] * 1e-12 # from ps to s
this_data = np.array([length, damping_time]).T

# Order data
index = np.argsort(length)
length = length[index]
damping_time = damping_time[index]
this_data = this_data[index, :]

#%%

# Define function to use while fitting
def tau_function(length, viscosity):
    tau = 2 * (length * density / (np.pi * viscosity))**2
    return tau

# Fit
fit_viscosity = iva.nonLinearFit(length, damping_time, tau_function, 
                                 showplot=False)[1][0]

#%%

x1=data[:, 2] * 1e-9
y1=data[:,6]

# Plot
plt.figure()
ax = plt.subplot()
plt.plot(x1, y1, '.')
#plt.plot(length, tau_function(length, fit_viscosity[0]), '-r')
plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'Frecuencia $GHz$ (s)')
plt.title(r'f vs. L')
plt.legend(['Datos'])
plt.grid(axis='x')