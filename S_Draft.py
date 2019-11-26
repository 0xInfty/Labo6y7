# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import numpy as np
import iv_save_module as ivs
import iv_analysis_module as iva
import matplotlib.pyplot as plt



# Parameters
density = 19300 # kg/m3
file = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1.txt'
file2 = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1_PostUSA.txt'
# Load data
data, header, footer = ivs.loadTxt(file)
data2, header2, footer2 = ivs.loadTxt(file2)

length = data[:, 2] * 1e-9 # from nm to m
length2 = data2[:, 2] * 1e-9 # from nm to m

#damping_time = data[:, 7] * 1e-12 # from ps to s
#this_data = np.array([length, damping_time]).T

# Order data
index = np.argsort(length)
index2 = np.argsort(length2)

length = length[index]
length2 = length2[index2]

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

x1 = length
y1 = data[:, 6]
x2 = length2
y2 = data2[:, 6]

#x1 = length
#y1 = damping_time

#x2 = x1
#y2 = tau_function(length, fit_viscosity[0])



# Plot
plt.figure()
ax = plt.subplot()
plt.plot(x1 , y1 , 'o''r')
plt.plot(x2, y2 , 'o''b')
plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'frecuencia $GHz$ (s)')
plt.title(r'Frecuencia vs Longitud')
plt.legend(['En aire', 'En Ta2O5'])
plt.grid(axis='x', which = 'both')

ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')

