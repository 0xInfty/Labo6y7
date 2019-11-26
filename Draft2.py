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

file = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1.txt'
file2 = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1_PostUSA.txt'

# Load data
data, header, footer = ivs.loadTxt(file)
data2, header, footer2 = ivs.loadTxt(file2)

#this_data = np.array([x, y]).T

# Parameters

rhoAu = 19300       # kg/m3
rhoTa = 8180        # kg/m3
gammaAu = 2e-3      # Pa/s

d  = data[:, 0]
r  = d/2   
A  = np.pi*(r**2)
L  = data[:, 2] * 1e-9 # from nm to m
L2 = data2[:, 2] * 1e-9 # from nm to m

w0 = data[:, 6] * 1e9 # from ps to s
w  = data2[:,6] * 1e9 # from ps to s

A = np.mean(A)
d = np.mean(A)
# Order data
index = np.argsort(L)
index2 = np.argsort(L2)

L = L[index]
L2 = L2[index2]

w0 = w0[index]
w  = w[index2]

young=78.43854513*1e9
sigmayoung=3.18683486e+09

G=4.32155921e+28
YG=43215.5921
sigmaG=1.90824023e+27
#%% FIT

# Define function to use while fitting
def freerod(L, young):
    return ((1/(2*L)))*(young/rhoAu)**(1/2)

def surroundedrod(L2,G):
    return (1/(2*np.pi))*np.sqrt((((1/(2*L2)))**2)*(young/rhoAu)+((G*2.75)/(rhoAu*A))-(d*np.pi*np.sqrt(rhoTa*G)+(np.pi**2*gammaAu/(2*L**2*rhoAu)))**2)
    
# Fit
popt, pcov = curve_fit(surroundedrod,L2,w,p0=3.66931003e+28)
print (popt)

sigma=np.sqrt(np.diag(pcov))
print (sigma)
#%% PLOT
# Plot

plt.figure()
ax = plt.subplot()

plt.plot(L2 , w , 'o''r')
plt.plot(L2 , surroundedrod(L2,G),'r')
plt.plot(L, w0 , 'o''b')
plt.plot(L , freerod(L,young),'b')

plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'frecuencia $GHz$ (s)')
plt.title(r'Frecuencia vs Longitud')
plt.legend(['En Ta2O5', 'ajuste','En aire', 'ajuste'])
plt.grid(axis='x', which = 'both')
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')