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
import scipy.stats as st
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as mcolors
#from scipy.optimize import curve_fit

#%% DATA

home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'

figs_folder = 'Informe L7\Figuras\Figuras análisis\Modelos (G, E, etc)'
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

# Results
#youngAu = 82.20e+9     #Pa/s      (Popt)
#stdyoungAu = 1.2e+09 #Young error [Pa/s]
#
#youngTa = 63.942e+9     #Pa/s      (Popt)
#stdyoungTa = 0.94e9   #Young error [Pa/s]
#
#meanG = 33.82e9
#stdG = 16.33e9

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

#%% EFFECTIVE YOUNG FITS

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
chiyoungAuFS = sum( (Ffreerod(L, youngAuFS) - f0)**2 ) / len(f0)

print('Young Au efectivo usando Fused-Silica: ' + 
      ivu.errorValueLatex(youngAuFS, stdyoungAuFS, #symbol='±',
                          units='Pa'))
print('Chi cuadrado: {:.2e}'.format(chiyoungAuFS))

youngAuTa2O5, stdyoungAuTa2O5 = iva.nonLinearFit(
        LL5, fL5, Ffreerod, showplot=False)[1][0]
chiyoungAuTa2O5 = sum( (Ffreerod(LL5, youngAuTa2O5) - fL5)**2 ) / len(fL5)

print('Young Au efectivo usando Ta2O5: ' + 
      ivu.errorValueLatex(youngAuTa2O5, stdyoungAuTa2O5, #symbol='±',
                          units='Pa'))
print('Chi cuadrado: {:.2e}'.format(chiyoungAuTa2O5))

#%% SHEAR MODULUS

def Gsurroundedrod_f0(f, d, f0, rhos):
    aux = rhoAu * d**2 * np.pi**2
    aux2 = ( 2.75 / np.pi ) - ( rhos / rhoAu )
    G = ( aux / aux2 ) * ( f**2 - f0**2 )
    return G

def Gsurroundedrod_young(f, d, L, youngAuFS, rhos):
    aux = rhoAu * d**2 * np.pi**2
    aux2 = ( 2.75 / np.pi ) - ( rhos / rhoAu )
    G = ( aux / aux2 ) * ( f**2 - youngAuFS / (4 * rhoAu * L**2) )
    return G

#G = ((f*2*np.pi)**2 - (f0*2*np.pi)**2) / ( 2.75/(rhoAu*A) - (((np.pi*r)/(rhoAu*A))**2)*rhoTa )        #surrounded rod for gamma = 0
#Gmeaned = ((np.mean(f)*2*np.pi)**2 - (np.mean(f0)*2*np.pi)**2) / ( 2.75/(rhoAu*np.mean(A)) - (((np.pi*np.mean(r))/(rhoAu*np.mean(A)))**2)*rhoTa )        #surrounded rod for gamma = 0

G = Gsurroundedrod_f0(f, d, f0, rhoTa)
G_young = Gsurroundedrod_young(f, d, L, youngAuFS, rhoTa)

Gmeaned = Gsurroundedrod_f0(np.mean(f), np.mean(d), np.mean(f0), rhoTa)
Gmeaned_young = Gsurroundedrod_young(np.mean(f), np.mean(d), np.mean(L), 
                                     youngAuFS, rhoTa)

#print(np.mean(G)*1e-9)
#print(np.std(G)*1e-9)
print('Módulo de corte G: ' + 
      ivu.errorValueLatex(np.mean(G), np.std(G), #symbol='±',
                          units='Pa', error_digits=1,))
print('Módulo de corte usando valores medios <G>: {:.2f} GPa'.format(
        Gmeaned*1e-9))
print('Módulo de corte G usando Eef: ' + 
      ivu.errorValueLatex(np.mean(G_young), np.std(G_young), #symbol='±',
                          units='Pa', error_digits=1,))
print('Módulo de corte usando Eef y valores medios <G>: {:.2f} GPa'.format(
        Gmeaned_young*1e-9))

# Try out a fit
def Fsurroundedrod_fit(L, G):
    aux = youngAuFS / (4 * rhoAu * L**2)
    aux2 = G / ( rhoAu * np.pi**2 * np.mean(d)**2 )
    aux3 = ( 2.75 / np.pi ) -  ( rhoTa / rhoAu )
    return np.sqrt( aux + aux2 * aux3 )

Gfit, stdGfit = iva.nonLinearFit(L2, f, Fsurroundedrod_fit, 
                                 initial_guess=(np.mean(G)),
                                 showplot=False)[1][0]
chiGfit = sum( (Fsurroundedrod_fit(L2, Gfit) - f)**2 ) / len(f)

print('Módulo de corte G ajustado: ' + 
      ivu.errorValueLatex(Gfit, stdGfit, #symbol='±',
                          units='Pa'))
print('Chi cuadrado: {:.2e}'.format(chiGfit))

# Let's study if the approximation was valid
rod_term = np.pi**2 * gammaAu * d / (4 * L**2)
medium_term = np.sqrt(rhoTa * Gfit)
print("External medium term is at least {:.2f} times bigger than \
inner medium term".format(np.min(medium_term)/np.max(rod_term)))
print("Using mean values, external medium term is {:.2f} times bigger than \
inner medium term".format(np.min(medium_term/rod_term)))

#%% YOUNG E USING G AND CL

alpha = rhoTa * (cLTa**2) / Gfit
youngTa_c = Gfit * ( 3 - ( 1 / (alpha-1) ) )

aux = youngTa_c - ( Gfit / ( ( alpha - 1 )**2 ) )
DyoungTa_c = aux * stdGfit

#aux = 3 - 3*Gfit - ( (rhoTa * cLTa**2) + Gfit) / ( (rhoTa * cLTa**2) - Gfit)
#DyoungTa_c = np.abs(np.mean(aux)) * stdGfit

#print(youngTa_c*1e-9)
print('Young Ta2O5 usando c: ' + 
      ivu.errorValueLatex(youngTa_c, DyoungTa_c,# symbol='±',
                          units='Pa'))

#%%

alpha_2 = rhoTa * (cLTa**2) / G
youngTa_c_2 = G * ( 3 - ( 1 / (alpha_2-1) ) )

#print(youngTa_c*1e-9)
print('Young Ta2O5 usando c: ' + 
      ivu.errorValueLatex(youngTa_c, np.std(youngTa_c_2),# symbol='±',
                          units='Pa'))

#%% G HISTOGRAM

nbins = (max(max(G), max(G_young)) - min(min(G), min(G_young))) / 2
nbins = nbins * np.cbrt(len(G) + len(G_young)) / np.mean([st.iqr(G), st.iqr(G_young)])
nbins = int(round(nbins))

if False: print('https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule')

fig, [axh, axbp] = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace':0})

# Histogram
n, bins, patches = axh.hist([G*1e-9, G_young*1e-9], nbins, density=100, 
                             alpha=.4, rwidth=.8, color=['b', 'r'])
del patches
axh.set_frame_on(False)
axh.yaxis.set_visible(False)

# Boxplot
bplot = axbp.boxplot(
    [G*1e-9, G_young*1e-9], 
    showmeans=True, meanline=True, 
    meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
    medianprops={'color':'k', 'linewidth':2},
    flierprops={'markersize':7},
    patch_artist=True,
    widths=.6,
    labels=[r'', r''],
    vert=False)
for p, c in zip(bplot['boxes'], ['blue','red']):
    p.set_facecolor(c) # paint white boxes
    p.set_alpha(.4)
del p, bplot

# Format
base_height = .2
base_width = .6
label_right_space = .1
label_left_space = .08
alpha = .25

# Axes size
box = axbp.get_position()
box.y1 = box.y0 + base_height
box.y0 = box.y0 + alpha
box.y1 = box.y1 + alpha
axbp.set_position(box)

box = axh.get_position()
box.y1 = box.y0 + .25
axh.set_position(box)

# Scale format
ylimsh = axh.get_ylim()
axh.autoscale(False)
ylimsbp = axbp.get_ylim()
for t in axbp.yaxis.get_major_ticks():
    t.set_visible(False)
box = axh.get_position()
box.x0 = box.x0 - .05 * box.width
axh.set_position(box)
box = axbp.get_position()
box.x0 = box.x0 - .05 * box.width
axbp.set_position(box)

# Add curves over it
axh.vlines(Gfit*1e-9, *ylimsh, linestyles='-.')
#axbp.vlines(np.mean(G)*1e-9, *ylimsbp, linestyles='solid',
#           label=r'$\overline{G}$')
#axbp.vlines(Gmeaned*1e-9, *ylimsbp, linestyles='dashed', 
#           label=r'$<G>$')

# Format plot
plt.xlabel("Módulo de corte G (GPa)")
axh.legend(['Ajuste', r'$G \rightarrow f_i$', r'$G \rightarrow E_{ef}$'],
           loc=(.75,.65))
#plt.ylabel(r"Densidad de probabilidad $\int f(F) dF = 1$")
#plt.legend()

# Grid's format
axbp.xaxis.set_minor_locator(AutoMinorLocator())
axbp.grid(which='major', axis='x')
axbp.grid(which='minor', axis='x', linestyle=':')
axbp.grid(which='major', axis='y')
axbp.yaxis.tick_right()
axbp.yaxis.set_label_position('right')

ivs.saveFig(os.path.join(home, figs_folder, 'HistogramasG.png'), 
            overwrite=True)

#%% PLOT L1 (FUSED SILICA): AIR VS TA2O5

x = np.linspace(L[0], L[-1], 50)

# Plot
plt.figure()
ax = plt.subplot()

data_plot = []
data_plot.append( plt.plot(L * 1e9, f0 * 1e-9 , '+''r', markersize=8,
                           label='FS + Aire')[0] )
data_plot.append( plt.plot(L2 * 1e9, f * 1e-9 , '.''b', markersize=7,
                           label='FS + Ta$_2$O$_5$')[0] )

lines = []
lines.append( plt.plot(x * 1e9, Ffreerod(x, youngAuFS) * 1e-9, 'r', 
                       label=('Ajuste Aire $E_{ef}$ = ' +
                              ivu.errorValueLatex(youngAuFS, 
                                                  stdyoungAuFS, 
                                                  units='Pa')))[0] )
#lines.append( plt.plot(x * 1e9, 
#                       Fsurroundedrod_young(x, np.mean(d), 
#                                            youngAuFS, 
#                                            np.mean(G), rhoTa) * 1e-9, 
#                        'b--',
#                        label=('Modelo Ta$_2$O$_5$ $G$ = ' +
#                               ivu.errorValueLatex(np.mean(G), np.std(G), 
#                                                   units='Pa')))[0] )
lines.append( plt.plot(x * 1e9, Fsurroundedrod_fit(x, Gfit) * 1e-9, 
                       'b',
                       label=('Ajuste Ta$_2$O$_5$ $G$ = ' +
                              ivu.errorValueLatex(Gfit, stdGfit, 
                                                  units='Pa')))[0] )
#ax.fill_between(x * 1e9, 
#                Fsurroundedrod_young(x, np.mean(d), youngAuFS, 
#                                     np.mean(G) - np.std(G), rhoTa) * 1e-9,
#                Fsurroundedrod_young(x, np.mean(d), youngAuFS, 
#                                     np.mean(G) + np.std(G), rhoTa) * 1e-9,
#                color='b',
#                alpha=0.1)

# Format 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia $F$ (GHz)')
leg = plt.legend(handles=lines, loc='lower left')
ax2 = plt.gca().add_artist(leg)
plt.legend(handles=data_plot, loc='upper right')
ylims = ax.get_ylim()
#new_ylims = (ylims[0] - .1 * (ylims[1] - ylims[0]), ylims[1])
#ax.set_ylim(new_ylims)

# Grid format
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')

ivs.saveFig(os.path.join(home, figs_folder, 'FvsLL1preypostUSA.png'), # _Calculus
            overwrite=True)

#%% PLOT L1 y L5 (YOUNG EFECTIVO): AIR VS TA2O5

x = np.linspace(L[0], L[-1], 50)
x2 = np.linspace(LL5[0], LL5[-1], 50)
xf = np.linspace(min(L[0],LL5[0]), max(L[-1],LL5[-1]), 50)

# Plot
plt.figure()
ax = plt.subplot()

data_plot = []
data_plot.append( plt.plot(L * 1e9, f0 * 1e-9 , '+''r', markersize=8,
                           label='FS + Aire')[0] )
data_plot.append( plt.plot(LL5 * 1e9, fL5 * 1e-9 , 'x', markersize=7,
                           color=mcolors.CSS4_COLORS['forestgreen'],
                           label='Ta$_2$O$_5$ + Aire')[0] )

lines = []
lines.append( plt.plot(x * 1e9, Ffreerod(x, youngAuFS) * 1e-9, 'r', 
                       label=(r'Ajuste Aire $E_{ef}$ = ' +
                              ivu.errorValueLatex(youngAuFS, 
                                                  stdyoungAuFS, 
                                                  units='Pa')))[0] )
plt.plot( xf * 1e9, Ffreerod(xf, youngAuFS) * 1e-9, 'r--' )
lines.append( plt.plot(x2 * 1e9, Ffreerod(x2, youngAuTa2O5) * 1e-9, 
                       color=mcolors.CSS4_COLORS['forestgreen'],
                       label=(r'Ajuste Ta$_2$O$_5$ $E_{ef}$ = ' +
                              ivu.errorValueLatex(youngAuTa2O5, 
                                                  stdyoungAuTa2O5, 
                                                  units='Pa')))[0] )
plt.plot( xf * 1e9, Ffreerod(xf, youngAuTa2O5) * 1e-9, '--',
          color=mcolors.CSS4_COLORS['forestgreen'])
#ax.fill_between(x * 1e9, 
#                Ffreerod(x, youngAuFS - stdyoungAuFS) * 1e-9,
#                Ffreerod(x, youngAuFS + stdyoungAuFS) * 1e-9,
#                color='r',
#                alpha=0.1)
#ax.fill_between(x * 1e9, 
#                Ffreerod(x, youngAuTa2O5 - stdyoungAuTa2O5) * 1e-9,
#                Ffreerod(x, youngAuTa2O5 + stdyoungAuTa2O5) * 1e-9,
#                color='b',
#                alpha=0.1)

# Format 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia $F$ (GHz)')
leg = plt.legend(handles=lines, loc='lower left')
ax2 = plt.gca().add_artist(leg)
plt.legend(handles=data_plot, loc='upper right')

# Grid format
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')

ivs.saveFig(os.path.join(home, figs_folder, 'FvsLL1yL5aire.png'), 
            overwrite=True)

#%% PLOT ALL

x = np.linspace(L[0], L[-1], 50)
x2 = np.linspace(LL5[0], LL5[-1], 50)
xf = np.linspace(min(L[0],LL5[0]), max(L[-1],LL5[-1]), 50)

# Plot
plt.figure()
ax = plt.subplot()

data_plot = []
data_plot.append( plt.plot(L * 1e9, f0 * 1e-9 , 'x''r', markersize=7,
                           label='FS + Aire')[0] )
data_plot.append( plt.plot(L2 * 1e9, f * 1e-9 , 'x''b', markersize=7,
                           label='FS + Ta$_2$O$_5$')[0] )
data_plot.append( plt.plot(LL5 * 1e9, fL5 * 1e-9 , 'x', markersize=7,
                           color=mcolors.CSS4_COLORS['forestgreen'],
                           label='Ta$_2$O$_5$ + Aire')[0] )

lines = []
lines.append( plt.plot(x * 1e9, Ffreerod(x, youngAuFS) * 1e-9, 'r', 
                       label=(r'Ajuste Aire $E_{ef}$ = ' +
                              ivu.errorValueLatex(youngAuFS, 
                                                  stdyoungAuFS, 
                                                  units='Pa')))[0] )
plt.plot( xf * 1e9, Ffreerod(xf, youngAuFS) * 1e-9, 'r--' )
lines.append( plt.plot(x * 1e9, Fsurroundedrod_fit(x, Gfit) * 1e-9, 
                       'b',
                       label=('Ajuste Ta$_2$O$_5$ $G$ = ' +
                              ivu.errorValueLatex(Gfit, stdGfit, 
                                                  units='Pa')))[0] )
lines.append( plt.plot(x2 * 1e9, Ffreerod(x2, youngAuTa2O5) * 1e-9, 
                       color=mcolors.CSS4_COLORS['forestgreen'],
                       label=(r'Ajuste Ta$_2$O$_5$ $E_{ef}$ = ' +
                              ivu.errorValueLatex(youngAuTa2O5, 
                                                  stdyoungAuTa2O5, 
                                                  units='Pa')))[0] )
plt.plot( xf * 1e9, Ffreerod(xf, youngAuTa2O5) * 1e-9, '--',
          color=mcolors.CSS4_COLORS['forestgreen'])
#ax.fill_between(x * 1e9, 
#                Ffreerod(x, youngAuFS - stdyoungAuFS) * 1e-9,
#                Ffreerod(x, youngAuFS + stdyoungAuFS) * 1e-9,
#                color='r',
#                alpha=0.1)
#ax.fill_between(x * 1e9, 
#                Ffreerod(x, youngAuTa2O5 - stdyoungAuTa2O5) * 1e-9,
#                Ffreerod(x, youngAuTa2O5 + stdyoungAuTa2O5) * 1e-9,
#                color='b',
#                alpha=0.1)

# Format 
plt.xlabel('Longitud $L$ (nm)')
plt.ylabel(r'Frecuencia $F$ (GHz)')
leg = plt.legend(handles=lines, loc='upper right')
ax2 = plt.gca().add_artist(leg)
plt.legend(handles=data_plot, loc='lower left')

# Grid format
ax.minorticks_on()
ax.tick_params(axis='y', which='minor', left=False)
ax.tick_params(length=5)
ax.grid(axis='x', which='both')
plt.grid(axis='y', which = 'both')

ivs.saveFig(os.path.join(home, figs_folder, 'FvsL.png'), 
            overwrite=True)

