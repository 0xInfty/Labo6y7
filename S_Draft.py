 # -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import numpy as np
import iv_save_module as ivs
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

<<<<<<< HEAD


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
=======
#%% DATA

new_file = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Análisis\ComparedAnalysis_FusedSilica/Resultados_Comparados.txt'
#file = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1_bis.txt'
#file2 = r'F:\Pump-Probe\Iván y Valeria\OneDrive\Labo 6 y 7\Análisis\Resultados_Totales_LIGO1_PostUSA.txt'

# Load data
new_data, header, footer = ivs.loadTxt(new_file)
#data, header, footer = ivs.loadTxt(file)
#data2, header, footer = ivs.loadTxt(file2)

#this_data = np.array([x, y]).T

# Parameters

rhoAu = 19.3e3 # kg/m3
rhoTa = 8.18e3 # kg/m3
Viscosity = 2e-3 # Pa/s for gold

r = new_data[:,0] * 1e-9 / 2
#r  = data[:, 0] * 1e-9 / 2    #radius
A  = np.pi*(r**2)
L = new_data[:,2]
#L  = data[:, 2] * 1e-9 # from nm to m
#L2 = data2[:, 2] * 1e-9 # from nm to m

#w0 = data[:, 6] * 2 * np.pi * 1e9 # from ps to s
#w  = data2[:,6] * 2 * np.pi * 1e9 # from ps to s
w0 = new_data[:, 6] * 2 * np.pi * 1e9 # from ps to s
w  = new_data[:,7] * 2 * np.pi * 1e9 # from ps to s

# Order data
index = np.argsort(L)
#index2 = np.argsort(L2)

L = L[index]
#L2 = L2[index2]

w0 = w0[index]
w = w[index]
#w  = w[index2]

r = r[index]

#%% Expresions

G = (w**2 - w0**2) / ( 2.75/(rhoAu*A) - (np.pi*r/(rhoAu*A))**2 * rhoTa )

>>>>>>> f6d51826c79ecf4d1460fbecbe1da8ba138be943

#%% FIT

# Define function to use while fitting
def freerod(L, young):
    return (1/(2*L))*((young/rhoAu)**(1/2))

def surroundedrod(w0,K1):
    return 
    
# Fit
popt, pcov = curve_fit(surroundedrod,L,w,p0=10)
print (popt *1e-9)

#%% PLOT

<<<<<<< HEAD
x1 = length
y1 = data[:, 6]
x2 = length2
y2 = data2[:, 6]

#x1 = length
#y1 = damping_time

#x2 = x1
#y2 = tau_function(length, fit_viscosity[0])


=======
x=
y=
>>>>>>> f6d51826c79ecf4d1460fbecbe1da8ba138be943

# Plot
plt.figure()
ax = plt.subplot()
<<<<<<< HEAD
plt.plot(x1 , y1 , 'o''r')
plt.plot(x2, y2 , 'o''b')
=======
plt.plot(x , y , 'o''r')
plt.plot(x , freerod(x,popt))
>>>>>>> f6d51826c79ecf4d1460fbecbe1da8ba138be943
plt.xlabel('Longitud $L$ (m)')
plt.ylabel(r'frecuencia $GHz$ (s)')
plt.title(r'Frecuencia vs Longitud')
plt.legend(['En aire', 'En Ta2O5'])
plt.grid(axis='x', which = 'both')
<<<<<<< HEAD

=======
>>>>>>> f6d51826c79ecf4d1460fbecbe1da8ba138be943
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
<<<<<<< HEAD

=======
>>>>>>> f6d51826c79ecf4d1460fbecbe1da8ba138be943
