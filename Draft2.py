 # -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import os
import numpy as np
import iv_save_module as ivs
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
#import scipy.stats as st
from scipy.optimize import curve_fit

#%% DATA

home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'

figs_folder = 'Informe L7\Figuras'
data_folder = 'Informe L7\Datos Iván'

file = os.path.join(home, data_folder, 
                    'Resultados_Comparados_LIGO1 sin outl.txt')
file2 = os.path.join(home, data_folder,
                     'Resultados_Comparados_LIGO1_PostUSA sin outl.txt')

# Load data
data, header, footer = ivs.loadTxt(file) # In air
data2, header, footer2 = ivs.loadTxt(file2) # In Ta2O5

# Parameters
rhoAu = 19300       # kg/m3
rhoTa = 8180        # kg/m3
gammaAu = 2e-3      # Pa/s
cLTa = 4920         # m/s

r = data[:,0] * 1e-9 / 2
A  = np.pi*(r**2)
L  = data[:, 2] * 1e-9 # from nm to m
L2 = data2[:, 2] * 1e-9 # from nm to m

f0 = data[:, 6] * 1e9 # from ps to s
f  = data2[:,6] * 1e9 # from ps to s

#RESULTS
youngAu = 82.20e+9     #Pa/s      (Popt)
stdyoungAu = 1.2e+09 #Young error [Pa/s]

youngTa = 63.942e+9     #Pa/s      (Popt)
stdyoungTa = 0.94e9   #Young error [Pa/s]

G = 33.82e9
stdG= 16.33e9
# Order data
    
index = np.argsort(L)#[2:]         elimina los 2 primeros números
index2 = np.argsort(L2)

L = L[index]
f0 = f0[index]

f  = f[index2]
L2 = L2[index2]

r = r[index]
A = A[index]

#%% Expresions

G = ((f*2*np.pi)**2 - (f0*2*np.pi)**2) / ( 2.75/(rhoAu*A) - (((np.pi*r)/(rhoAu*A))**2)*rhoTa )        #surrounded rod for gamma = 0
Gmeaned = ((np.mean(f)*2*np.pi)**2 - (np.mean(f0)*2*np.pi)**2) / ( 2.75/(rhoAu*np.mean(A)) - (((np.pi*np.mean(r))/(rhoAu*np.mean(A)))**2)*rhoTa )        #surrounded rod for gamma = 0

f0 = (1/(2*L))*((youngAu/rhoAu)**(1/2))                                         #free rod
f = (1/(2*np.pi))*np.sqrt((((1/(2*L2)))**2)*(youngAu/rhoAu)+((G*2.75)/(rhoAu*A))-(2*r*np.pi*np.sqrt(rhoTa*G)/(2*rhoAu*A)+(np.pi**2*gammaAu/(2*L2**2*rhoAu)))**2)
                                                                                #surrounded rod
print(np.mean(G)*1e-9)
print(np.std(G)*1e-9)
print(Gmeaned*1e-9)

print('Módulo de corte G: ' + 
      ivu.errorValueLatex(np.mean(G), np.std(G), units='Pa', symbol='±'))


#%%

alpha = rhoTa * (cLTa**2) / np.mean(G)
Dalpha = rhoTa * (cLTa**2) * np.std(G) / (np.mean(G)**2)
youngTa_c = np.mean(G) * ( 3 - ( 1 / (alpha-1) ) )
DyoungTa_c = np.std(G) * ( 3 - ( 1 / (alpha-1) ) )

#print(youngTa_c*1e-9)

print('Young usando c: ' + 
      ivu.errorValueLatex(youngTa_c, DyoungTa_c, units='Pa', symbol='±'))

#%% FIT

def freerod(L, young):
    return ((1/(2*L)))*(young/rhoAu)**(1/2)

def surroundedrod(L,G):
    return (1/(2*np.pi))*np.sqrt((((1/(2*L)))**2)*(youngAu/rhoAu)+((G*2.75)/(rhoAu*A))-((np.pi**2*gammaAu/(2*L**2*rhoAu)))**2)
    
# Fit
popt, pcov = curve_fit(freerod,L,f0)
sigma = np.sqrt(np.diag(pcov))
print (popt *1e-9,sigma *1e-9)

#%% PLOT

x = np.linspace(L[0],L[-1],1000)
x2 = np.linspace(L2[0],L2[-1],1000)

# Plot
plt.figure()
ax = plt.subplot()

plt.plot(L * 1e9 , f0 * 1e-9 , 'x''r', 
         label='Fused Silica + Aire')
plt.plot(L2 * 1e9 , f *  1e-9 , 'x''b', 
         label=r'Fused Silica + Ta$_2$O$_5$')
plt.plot(x * 1e9 , freerod(x,youngAu) * 1e-9, 'r', 
         label=('Ajuste Aire $E_{ef}$=' +
                ivu.errorValueLatex(youngAu, stdyoungAu, units='Pa')))
#plt.plot(x * 1e9 , freerod(x,youngTa) * 1e-9, 'b',
#         label=('Ajuste Ta$_2$O$_5$ $E_{ef}$=' +
#                ivu.errorValueLatex(youngTa, stdyoungTa, units='Pa')))
 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia $F$ (GHz)')
plt.title(r'Frecuencia vs Longitud')
plt.legend()

ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')

#%%

ivs.saveFig

#%% HISTOGRAM
plt.figure()
ax = plt.subplot()
n, bins, patches = plt.hist(G*1e-9, bins = 5, density=True,alpha=0.8, facecolor='blue', rwidth=0.95)
del patches

# Add curve over it
#x = np.linspace(np.min(bins), np.max(bins), 50)
#plt.plot(x,st.norm.pdf(x,np.mean(G)*1e-9,np.std(G)*1e-9),'k')
plt.vlines(x=32.6341121726834,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],linestyles='dashed')
plt.vlines(x=np.mean(G)*1e-9,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],linestyles='solid')

# Format plot
plt.xlabel("Módulo de corte G (GPa)")
plt.ylabel(r"Densidad de probabilidad $\int f(F) dF = 1$")
plt.legend(['Gaussiana asociada', 'G promediado','Valor medio'])