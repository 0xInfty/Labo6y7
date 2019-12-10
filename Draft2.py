 # -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import numpy as np
import iv_save_module as ivs
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#%% DATA

file = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO5bis.txt'
file2 = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1.txt'

# Load data
data, header, footer = ivs.loadTxt(file)
data2, header, footer = ivs.loadTxt(file2)

# Parameters

rhoAu = 19.3e3 # kg/m3
rhoTa = 8.18e3 # kg/m3
Viscosity = 2e-3 # Pa/s for gold
young = 63.942079 * 1e9 #Pa/s

r = data[:,0] * 1e-9 / 2
A  = np.pi*(r**2)
L  = data[:, 2] * 1e-9 # from nm to m
L2 = data2[:, 2] * 1e-9 # from nm to m

w0 = data[:, 6] * 1e9 # from ps to s
w  = data2[:,6] * 1e9 # from ps to s

index = np.argsort(w)
L2 = L2[index[4:]]
w = w[index[4:]]

newL = np.array(list(L)+list(L2))
L = newL

neww = np.array(list(w0)+list(w))
w0 = neww

# Order data
index = np.argsort(L)
#index2 = np.argsort(L2)

L = L[index]
#L2 = L2[index2]

w0 = w0[index]
#w  = w[index2]

#r = r[index]

#%% Expresions

G = (w**2 - w0**2) / ( 2.75/(rhoAu*A) - (np.pi*r/(rhoAu*A))**2 * rhoTa )        #surrounded rod for gamma = 0
#w0 = (1/(2*L))*((young/rhoAu)**(1/2))                                           #free rod


#%% FIT

# Define function to use while fitting
def freerod(L,young):
    return (1/(2*L))*((young/rhoAu)**(1/2))

#def surroundedrod(w0,K1):
#    return 
    
# Fit
popt, pcov = curve_fit(freerod,L,w0,p0=78)
print (popt *1e-9)

#%% PLOT

x=L
y=w0

# Plot
plt.figure()
plt.plot(x * 1e9 , y * 1e-9 , 'o''r')
plt.plot(L2 * 1e9 , w * 1e-9 , 'o''b')
plt.plot(L * 1e9, freerod(L,popt) * 1e-9, 'k')
#plt.plot(x , freerod(x,popt))
plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'Frecuencia (GHz)')
plt.title(r'Frecuencia vs Longitud')
plt.legend(['Ta2O5', 'SiO2','ajuste'])
plt.grid()

