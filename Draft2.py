 # -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 14:34:41 2019

@author: Lec
"""

import os
import numpy as np
import iv_save_module as ivs
import iv_analysis_module as iva
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
#import scipy.stats as st
#from scipy.optimize import curve_fit

#%% DATA

home = r'C:\Users\Usuario\OneDrive\Labo 6 y 7\OneDrive\Labo 6 y 7'

figs_folder = 'Informe L7\Figuras\Figuras análisis\Modelos L1 (G, E, etc)'
data_folder = 'Informe L7\Datos Iván'

file = os.path.join(home, data_folder, 
                    'Resultados_Comparados_LIGO1 sin outl.txt')
file2 = os.path.join(home, data_folder,
                     'Resultados_Comparados_LIGO1_PostUSA sin outl.txt')
file3 = os.path.join(home, data_folder,
                     'Resultados_Comparados_LIGO5bis.txt')

# Load data
data, header, footer = ivs.loadTxt(file) # Fused Silica + Air
data2, header, footer2 = ivs.loadTxt(file2) # Fused Silica + Ta2O5
data3, header, footer3 = ivs.loadTxt(file3) # Ta2O5 + Air

# Parameters
rhoAu = 19300       # kg/m3
rhoTa = 8180        # kg/m3
gammaAu = 2e-3      # Pa/s
cLTa = 4920         # m/s

f0 = data[:, 6] * 1e9 # from GHz to Hz
d = data[:,0] * 1e-9
L  = data[:, 2] * 1e-9 # from nm to m

f  = data2[:,6] * 1e9
d2 = data2[:,0] * 1e-9
L2 = data2[:, 2] * 1e-9

fL5 = data3[:,6] * 1e9
dL5 = data3[:,0] * 1e-9
LL5 = data3[:,2] * 1e-9

#Results
youngAu = 82.20e+9     #Pa/s      (Popt)
stdyoungAu = 1.2e+09 #Young error [Pa/s]

youngTa = 63.942e+9     #Pa/s      (Popt)
stdyoungTa = 0.94e9  #Young error [Pa/s]

meanG = 33.82e9         #Pa/s
stdG = 16.33e9       #G error [Pa/s]

# Order data  
index = np.argsort(L)#[2:]         elimina los 2 primeros números
index2 = np.argsort(L2)
index3 = np.argsort(LL5)

L = L[index]
d = d[index]
f0 = f0[index]

f  = f[index2]
L2 = L2[index2]

fL5 = fL5[index3]
LL5 = LL5[index3]

#%% Expresions

def Gsurroundedrod(f, f0, d, rhos):
    aux = rhoAu * d**2 * np.pi**2
    aux2 = ( 2.75 / np.pi ) - ( rhos / rhoAu )
    G = ( aux / aux2 ) * ( f**2 - f0**2 )
    return G

#G = ((f*2*np.pi)**2 - (f0*2*np.pi)**2) / ( 2.75/(rhoAu*A) - (((np.pi*r)/(rhoAu*A))**2)*rhoTa )        #surrounded rod for gamma = 0
#Gmeaned = ((np.mean(f)*2*np.pi)**2 - (np.mean(f0)*2*np.pi)**2) / ( 2.75/(rhoAu*np.mean(A)) - (((np.pi*np.mean(r))/(rhoAu*np.mean(A)))**2)*rhoTa )        #surrounded rod for gamma = 0

G = Gsurroundedrod(f, f0, d, rhoTa)
Gmeaned = Gsurroundedrod(np.mean(f), np.mean(f0), np.mean(d), rhoTa)

#print(np.mean(G)*1e-9)
#print(np.std(G)*1e-9)
print('Módulo de corte G: ' + 
      ivu.errorValueLatex(np.mean(G), np.std(G), units='Pa', symbol='±'))
print('Módulo de corte usando valores medios <G>: {:.2f} GPa'.format(Gmeaned*1e-9))

# Try out a fit
def Fsurroundedrod_fit(L, G):
    aux = youngAuFS / (4 * rhoAu * L**2)
    aux2 = G / ( rhoAu * np.pi**2 * np.mean(d)**2 )
    aux3 = ( 2.75 / np.pi ) -  ( rhoTa / rhoAu )
    return np.sqrt( aux + aux2 * aux3 )

Gfit, stdGfit = iva.nonLinearFit(
        L2, f, Fsurroundedrod_fit, initial_guess=(np.mean(G)))[1][0]
print('Módulo de corte G: ' + 
      ivu.errorValueLatex(Gfit, stdGfit, units='Pa', symbol='±'))

#%% CALCULUS OF E USING G AND CL

alpha = rhoTa * (cLTa**2) / np.mean(G)
Dalpha = rhoTa * (cLTa**2) * np.std(G) / (np.mean(G)**2)
youngTa_c = np.mean(G) * ( 3 - ( 1 / (alpha-1) ) )
DyoungTa_c = np.std(G) * ( 3 - ( 1 / (alpha-1) ) )

#print(youngTa_c*1e-9)

print('Young Ta2O5 usando c: ' + 
      ivu.errorValueLatex(youngTa_c, DyoungTa_c, units='Pa', symbol='±'))

#%% FIT

def Ffreerod(L, young):
    return ( 1 / (2*L) ) * np.sqrt( young / rhoAu )

def Fsurroundedrod_young(L, d, young, G, rhos):
    aux = young / (4 * rhoAu * L**2)
    aux2 = G / ( rhoAu * np.pi**2 * d**2 )
    aux3 = ( 2.75 / np.pi ) -  ( rhos / rhoAu )
    return np.sqrt( aux + aux2 * aux3 )

def Fsurroundedrod_fair(f0, d, G, rhos):
    aux2 = G / ( rhoAu * np.pi**2 * d**2 )
    aux3 = ( 2.75 / np.pi ) -  ( rhos / rhoAu )
    return np.sqrt( f0**2 + aux2 * aux3 )

# (1/(2*np.pi))*np.sqrt((((1/(2*L2)))**2)*(youngAu/rhoAu)+ ((G*2.75)/(rhoAu*A)) - (2*r*np.pi*np.sqrt(rhoTa*G)/(2*rhoAu*A) + (np.pi**2*gammaAu/(2*L2**2*rhoAu)))**2)

# Fit
youngAuFS, stdyoungAuFS = iva.nonLinearFit(L, f0, Ffreerod, showplot=False)[1][0]
print('Young Au usando Fused-Silica: ' + 
      ivu.errorValueLatex(youngAuFS, stdyoungAuFS, units='Pa', symbol='±'))

youngAuTa2O5, stdyoungAuTa2O5 = iva.nonLinearFit(
        LL5, fL5, Ffreerod, showplot=False)[1][0]
print('Young Au usando Ta2O5: ' + 
      ivu.errorValueLatex(youngAuTa2O5, stdyoungAuTa2O5, units='Pa', symbol='±'))

#%% PLOT

x = np.linspace(L[0],L[-1],1000)
xd= np.linspace(L[-1],LL5[-1],1000)

x2 = np.linspace(LL5[0],LL5[-1],1000)
x2d= np.linspace(L[0],LL5[0],1000)

# Plot
plt.figure()
ax = plt.subplot()

data1, = plt.plot(L * 1e9, f0 * 1e-9 , 'x''r', markersize=7,
                  label='Fused Silica + Aire')
line1, = plt.plot(x * 1e9, Ffreerod(x, youngAu) * 1e-9, 'r', linewidth=1.5,
                  label=('Ajuste Aire $E_{ef}$ = ' + ivu.errorValueLatex(youngAuFS, stdyoungAuFS, units='Pa')))
dashline1, = plt.plot(xd * 1e9, Ffreerod(xd, youngAu) * 1e-9, '--''r', linewidth=1.5)
                      
#line2, = plt.plot(x * 1e9, Fsurroundedrod_young(x, np.mean(d), youngAuFS, 
#                                                np.mean(G), rhoTa) * 1e-9, 
#                 'b--',
#                 label=('Modelo Ta$_2$O$_5$ $G$ = ' +
#                        ivu.errorValueLatex(np.mean(G), np.std(G), units='Pa')))

data3, = plt.plot(LL5 * 1e9, fL5 *  1e-9 , 'x''b', markersize=7, label=r'Ta$_2$O$_5$ + Aire')
line3, = plt.plot(x2 * 1e9, Ffreerod(x2, youngTa) * 1e-9, 'b', linewidth=1.5,
                 label=('Ajuste Ta$_2$O$_5$ $E_{ef}$ = ' + ivu.errorValueLatex(youngTa, stdyoungTa, units='Pa')))
dashline3, = plt.plot(x2d * 1e9, Ffreerod(x2d, youngTa) * 1e-9, '--''b', linewidth=1.5)

#ax.fill_between(x * 1e9, 
#                Fsurroundedrod_young(x, np.mean(d), youngAuFS, 
#                                     np.mean(G) - np.std(G), rhoTa) * 1e-9,
#                Fsurroundedrod_young(x, np.mean(d), youngAuFS, 
#                                     np.mean(G) + np.std(G), rhoTa) * 1e-9,
#                color='b',
#                alpha=0.1)
 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia $F$ (GHz)')
plt.title(r'Frecuencia vs Longitud')
leg = plt.legend(handles=[line1, line3], loc='lower left')
ax2 = plt.gca().add_artist(leg)
plt.legend(handles=[data1, data3], loc='upper right')

ax = plt.subplot()
plt.xticks()
plt.yticks()
ax.minorticks_on()
#ax.tick_params(axis='y', which='minor', left=False)
#ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')

ivs.saveFig(os.path.join(home, figs_folder, 'F L1 pre y post USA.png'), 
            overwrite=True)

#%%



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