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

file = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\ComparedAnalysis_FusedSilica\Resultados_Comparados_LIGO1.txt'
file2 = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7\Análisis\ComparedAnalysis_FusedSilica\Resultados_Comparados_LIGO1_PostUSA.txt'

# Load data
data, header, footer = ivs.loadTxt(file) # In air
data2, header, footer2 = ivs.loadTxt(file2) # In Ta2O5

filter_outliers = False
filter_notcommon = False

# Parameters
rhoAu = 19300       # kg/m3
rhoTa = 8180        # kg/m3
gammaAu = 2e-3      # Pa/s

r = data[:,0] * 1e-9 / 2
A  = np.pi*(r**2)
L  = data[:, 2] * 1e-9 # from nm to m
L2 = data2[:, 2] * 1e-9 # from nm to m

w0 = data[:, 6] * 1e9 # from ps to s
w  = data2[:,6] * 1e9 # from ps to s


if filter_notcommon:
    
    notcommon_index = footer['rods'].index('9,10')
    
    index = list(range(len(L)))
    index.remove(notcommon_index)
    
    w = w[index]
    L = L[index]
    
    r = r[index]
    A = A[index]
    

if filter_outliers:
    
    # Order data
    
    index = np.argsort(L)#[2:]         elimina los 2 primeros números
    index2 = np.argsort(L2)
    
    L = L[index]
    w0 = w0[index]
    
    w  = w[index2]
    L2 = L2[index2]
    
    r = r[index]
    A = A[index]

#RESULTS
youngAu = 82.20e+9     #Pa/s      (Popt)
sigmayoungAu = 1.2e+09 #Young error [Pa/s]

youngTa = 63.942e+9     #Pa/s      (Popt)
sigmayoungTa = 0.94e9   #Young error [Pa/s]


#%% Expresions

G = (w**2 - w0**2) / ( 2.75/(rhoAu*A) - (np.pi*r/(rhoAu*A))**2 * rhoTa )        #surrounded rod for gamma = 0
w0 = (1/(2*L))*((youngAu/rhoAu)**(1/2))                                           #free rod
w = (1/(2*np.pi))*np.sqrt((((1/(2*L2)))**2)*(young/rhoAu)+((G*2.75)/(rhoAu*A))-(2*r*np.pi*np.sqrt(rhoTa*G)/(2*rhoAu*A)+(np.pi**2*gammaAu/(2*L2**2*rhoAu)))**2)
                                                                                #surrounded rod


#%% FIT

def freerod(L, young):
    return ((1/(2*L)))*(young/rhoAu)**(1/2)

def surroundedrod(L,G):
    return (1/(2*np.pi))*np.sqrt((((1/(2*L)))**2)*(youngAu/rhoAu)+((G*2.75)/(rhoAu*A))-((np.pi**2*gammaAu/(2*L**2*rhoAu)))**2)
    
# Fit
popt, pcov = curve_fit(freerod,L,w0)
sigma = np.sqrt(np.diag(pcov))
print (popt *1e-9,sigma *1e-9)

#%% PLOT

x = np.linspace(L[0],L2[-1],1000)
x2 = np.linspace(L2[0],L2[-1],1000)

# Plot
plt.figure()

plt.plot(L * 1e9 , w0 * 1e-9 , 'x''r')
plt.plot(L2 * 1e9 , w *  1e-9 , 'x''b')
plt.plot(x * 1e9 , freerod(x,youngAu) * 1e-9, 'r')
plt.plot(x * 1e9 , freerod(x,youngTa) * 1e-9, 'b')
 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia (GHz)')
plt.title(r'Frecuencia vs Longitud')
plt.legend(['En SiO2', 'ajuste','En Ta2O5', 'ajuste'])

ax = plt.subplot()
plt.xticks()
plt.yticks()
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')