# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: LEC
"""

import iv_save_module as ivs
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import os
from tkinter import Tk, messagebox

def roundMatlab(x, Matlab_round_needed=True):
    """Returns round value in Matlab 2014's style.
    
    In Pyhon 3.7.3...
    >> round(80.5) = 80
    >> round(81.5) = 82
    
    But in Matlab 2014...
    >> round(80.5) = 81
    >> round(81.5) = 82
    
    Parameters
    ----------
    x : float
        Number to be rounded.
    Matlab_round_needed=True : bool
        Whether your Python version needs this function to round like Matlab 
        or not. Python 3.7.3 needs it.
    
    Returns
    -------
    y : int
        Rounded number.
    """
    
    xround = int(x)
    even = xround/2 == int(xround/2) # True if multiple of 2
    
    if even:
        y = round(x) + 1
    else:
        y = round(x)
    
    if Matlab_round_needed:
        return y
    else:
        return round(x)

#%%

# Parameters
dt = 2 # in ps
filename = os.path.join(os.getcwd(), 'Datos.txt')
autosave = True
Matlab_round_needed = True # Python 3.7.3 needs it

#%%

# Get data
t, x = np.loadtxt(filename)
T = max(t) - min(t)

## Kind of answer we want...
#x = c1.*exp(-b1.*t).*cos(w1*t+p1) + c2.*exp(-b2.*t).*cos(w2.*t+p2) + 
#    + c3.*exp(-b3.*t.^2) + noise + coherent_artifact;

## To do that we could define...
#cn = 8e-15
#noise = []
#for i in range(50):
#    noise.append(np.cos(2*pi*np.random.rand() * (t+dt) / (2*dt)))
#noise = np.array(noise)
#noise = cn * sum(noise)
#coherent_artifact = 1e-15 * np.exp(-t**2 / 2 / 0.07**2)
# Proporcional to cross-correlation

#%% ---------------------------------------------------------------------------
# FREQUENCIES AND DAMPING FACTORS
# -----------------------------------------------------------------------------

# Create data matrix
N = len(x)
M = round(0.75 * N)
X = np.array([x[i+j+1] for j in range(N-M) for i in range(M)]).reshape((N-M,M))

# Diagonalize square data matrix
[eigenvalues, eigenvectors] = np.linalg.eig( np.matmul(X, X.T) )
ordered_index = eigenvalues.argsort() # From smallest to largest value 
eigenvalues = eigenvalues[ordered_index] # Eigenvalues
eigenvectors = eigenvectors[:, ordered_index] # Eigenvectors on columns
#eigenvectors = np.array([l/l[0] for l in eigenvectors.T]).T # Normalize
rank = np.linalg.matrix_rank(np.diag(eigenvalues)) # Size measure

"""
La segunda columna de eV la tengo con signo invertido respecto a Matlab.
==> Obtengo resultados distintos para U.
"""

# Choose number of significant values
Nsignificant = 4
#plt.figure()
#plt.semilogy(eigenvalues, linestyle='none', marker='o', 
#             fillstyle='none', markersize=10)
#plt.xlabel("Número")
#plt.ylabel("Autovalores")
#Nsignificant = int(input('¿Número de valores singulares?\t'))

# Crop data according to it
F = np.zeros((N-M, N-M))
F[-Nsignificant:,-Nsignificant:] = np.diag(1/np.sqrt(eigenvalues[-Nsignificant:]))
auxiliar = np.matmul(eigenvectors, F)
U = np.matmul(X.T, auxiliar) # Xmatrix.T * eigenvectors * F

# Define polinomyal
auxiliar = np.matmul(eigenvectors.T, x[:N-M])
auxiliar = np.matmul(F.T, auxiliar)
A = np.matmul(U, auxiliar) # U * F.T * eigenvectors.T * xvector 
# |--> Least-Squares?
coeficients = np.array([1, *list(-A)])

# Solve and find its roots
roots = np.roots(coeficients)
ordered_index = abs(roots).argsort() 
roots = roots[ordered_index][::-1] # From largest to smallest absolute value

# Calculate damping constants 'b' and frequencies 'omega'
damping_constants = (np.log(abs(roots)) / dt)[:rank] # Crop them accordingly
angular_frequencies = (np.angle(roots) / dt)[:rank]

#%%

# Sort them
ordered_index = angular_frequencies.argsort() # From smallest to largest freq
angular_frequencies = angular_frequencies[ordered_index]
damping_constants = damping_constants[ordered_index]

# Crop them according to number of real roots and rank of diagonalized matrix
Nzeros = len(angular_frequencies) - np.count_nonzero(angular_frequencies)
angular_frequencies = abs(angular_frequencies)[:int(roundMatlab(
        (rank-Nzeros)/2+Nzeros,
        Matlab_round_needed))]
damping_constants = damping_constants[:int(roundMatlab(
        (rank-Nzeros)/2+Nzeros,
        Matlab_round_needed))]

# Then crop them according to the number of positive or zero damping constants
Npositives = len(damping_constants[damping_constants>=0])
ordered_index = damping_constants.argsort()[::-1] # From largest to smallest
damping_constants = damping_constants[ordered_index][:Npositives]
angular_frequencies = angular_frequencies[ordered_index][:Npositives]

# Now I have the smallest frequencies and largest damping constants
# Then I calculate the characteristic time tau and the quality factor Q
quality_factors = angular_frequencies / (2*damping_constants)
characteristic_times = 1 / damping_constants # in ps

#%% ---------------------------------------------------------------------------
# AMPLITUDES AND PHASES
# -----------------------------------------------------------------------------

# Create modelled data matrix
Nfit_terms = len(angular_frequencies)
t2 = np.arange(0, N*dt, dt) # Time starting on zero
X2 = np.zeros((N, 2*Nfit_terms))
for i, b, omega in zip(range(Nfit_terms), 
                       damping_constants, 
                       angular_frequencies):
    X2[:, 2*i] = np.exp(-b*t2) * np.cos(omega*t2)
    X2[:, 2*i+1] = -np.exp(-b*t2) * np.sin(omega*t2)

# Diagonalize square Hermitian modelled data matrix
[eigenvalues2, eigenvectors2] = np.linalg.eigh( np.matmul(X2, X2.T) )
ordered_index = eigenvalues2.argsort() # From smallest to largest absolute
eigenvalues2 = eigenvalues2[ordered_index] # Eigenvalues
eigenvectors2 = eigenvectors2[:, ordered_index] # Eigenvectors on columns
#eigenvectors2 = np.array([l/l[0] for l in eigenvectors2.T]).T # Normalize
rank2 = np.linalg.matrix_rank(np.diag(eigenvalues2)) # Size measure

# Choose number of significant values
Nsignificant2 = np.linalg.matrix_rank( np.matmul(X2, X2.T) )

# Crop data according to it
F2 = np.zeros((N, N))
F2[-Nsignificant2:,-Nsignificant2:] = np.diag(
        1/np.sqrt(eigenvalues2[-Nsignificant2:]))
auxiliar = np.matmul(eigenvectors2, F2)
U2 = np.matmul(X2.T, auxiliar) # Xmatrix.T * eigenvectors * F

# Get defining vector
auxiliar = np.matmul(eigenvectors2.T, x)
auxiliar = np.matmul(F2.T, auxiliar)
A2 = np.matmul(U2, auxiliar) # U * F.T * eigenvectors.T * xvector 
# |--> Least-Squares?

# Calculate phases 'phi' and amplitudes 'C'
amplitudes = []
phases = []
for i in range(Nfit_terms):
    
    if A2[2*i]==0 and A2[2*i+1]==0:
        amplitudes.append( 0 )
        phases.append( 0 )
    elif A2[2*i]==0:
        amplitudes.append( abs(A2[2*i+1]) )
        phases.append( np.sign(A2[2*i+1]) * pi/2 )
    elif A2[2*i+1]==0:
        amplitudes.append( abs(A2[2*i]) )
        phases.append( (1-np.sign(A2[2*i])) * pi/2 )
    else:
        amplitudes.append( np.sqrt(A2[2*i+1]**2 + A2[2*i]**2) )
        phases.append( np.arctan2(A2[2*i+1], A2[2*i]) )
amplitudes = np.array(amplitudes)
phases = np.array(phases)
pi_phases = phases / pi # in radians written as multiples of pi

#%% ---------------------------------------------------------------------------
# FIT, PLOTS AND STATISTICS
# -----------------------------------------------------------------------------

fit_terms = np.array([a * np.exp(-b*(t-t[0])) * np.cos(omega*(t-t[0]) + phi)
                     for a, b, omega, phi in zip(amplitudes,
                                                 damping_constants,
                                                 angular_frequencies,
                                                 phases)]).T
fit = sum(fit_terms.T)
chi_squared = sum( (fit - x)**2 ) / N # Best if absolute is smaller

## Statistics of the residue
#residue = x - fit
#residue_transform = abs(np.fft.rfft(residue))
#residue_frequencies = 1000 * np.fft.rfftfreq(N, d=dt) # in GHz
#plt.plot(residue_frequencies, residue_transform)

# Raman-like Spectrum parameters
frequencies_damping = 1000 * damping_constants / (2*pi) # in GHz
frequencies = 1000 * angular_frequencies / (2*pi) # in GHz
max_frequency = max(frequencies)
raman_frequencies = np.arange(0, 1.5*max_frequency, max_frequency/1000)

# Raman-like Spectrum per se
raman_spectrum_terms = np.zeros((len(raman_frequencies), Nfit_terms))
for i in range(Nfit_terms):
   if angular_frequencies[i]==0:
      raman_spectrum_terms[:,i] = 0
   else:
      raman_spectrum_terms[:,i] = amplitudes[i] * np.imag( 
         frequencies[i] / 
         (frequencies[i]**2 - raman_frequencies**2 - 
         2j * raman_frequencies * frequencies_damping[i]))
raman_spectrum = np.sum(raman_spectrum_terms, axis=1)

# What I would like this function to return
results = np.array([frequencies,
                    characteristic_times,
                    quality_factors,
                    amplitudes,
                    pi_phases]).T

others = dict(
        fit = np.array([t, x, fit, *list(fit_terms.T)]).T,
        raman = np.array([raman_frequencies,
                          raman_spectrum,
                          *list(raman_spectrum_terms.T)]).T,
        chi_squared = chi_squared
        )

#%%

# First I deglose data
fit = others['fit']
raman = others['raman']
Nfit_terms = fit.shape[1] - 3

# In order to save, if needed, I will need...
name = os.path.split(os.path.splitext(filename)[0])[-1]
path = os.path.split(filename)[0]

# Then, to plot, I first start a figure
plt.figure()
grid = plt.GridSpec(3, 5, hspace=0.1)

# In the upper subplot, I put the Raman-like spectrum
ax_spectrum = plt.subplot(grid[0,:4])
plt.plot(raman[:,0], raman[:,1], linewidth=2)
lspectrum_terms = plt.plot(raman[:,0], raman[:,2:], 
                           linewidth=2)
for l in lspectrum_terms: l.set_visible(False)
plt.xlabel("Frecuencia (GHz)")
plt.ylabel("Amplitud (u.a.)")
ax_spectrum.xaxis.tick_top()
ax_spectrum.xaxis.set_label_position('top')

# In the lower subplot, I put the data and fit
ax_data = plt.subplot(grid[1:,:])
ldata, = plt.plot(fit[:,0], fit[:,1], 'k', linewidth=0.5)
ax_data.autoscale(False)
lfit, = plt.plot(fit[:,0], fit[:,2], linewidth=2)
lfit_terms = plt.plot(fit[:,0], fit[:,3:], linewidth=2)
for l in lfit_terms: l.set_visible(False)
plt.xlabel("Tiempo (ps)")
plt.ylabel(r"Voltaje ($\mu$V)")

# Because it's pretty, I make an interactive legend
ax_legend = plt.axes([0.75, 0.642, 0.155, 0.24])
check_legend = wid.CheckButtons(ax_legend, ('Data', 
                                 'Ajuste', 
                                 *['Término {:.0f}'.format(i+1) 
                                 for i in range(Nfit_terms)]), 
                    (True, True, *[False for i in range(Nfit_terms)]))
check_legend.labels[1].set_color(lfit.get_color())
for leg, l in zip(check_legend.labels[2:], lfit_terms):
    leg.set_color(l.get_color())

# For that, I'll need a callback function  
def check_legend_callback(label):
    if label == 'Data':
        ldata.set_visible(not ldata.get_visible())
    elif label == 'Ajuste':
        lfit.set_visible(not lfit.get_visible())
    else:
        for i in range(Nfit_terms):
            if label == 'Término {:.0f}'.format(i+1):
                lfit_terms[i].set_visible(not lfit_terms[i].get_visible())
                lspectrum_terms[i].set_visible(
                        not lspectrum_terms[i].get_visible())
    plt.draw()
check_legend.on_clicked(check_legend_callback)

# Since I can, I would also like an interactive 'Save' button
ax_save = plt.axes([0.8, 0.01, 0.1, 0.04])
check_save = wid.Button(ax_save, 'Guardar')

# For that, I'll need another callback function
def check_save_callback(event):
    Tk().withdraw()
#        tk.newfilename = askopenfilename()
    newpath = os.path.join(path, 'Figuras')
    if not os.path.isdir(newpath):
        os.makedirs(newpath)
    newfilename = ivs.freeFile(os.path.join(newpath, name+'_fit.png'),
                               newformat='{}__{}')
    ax_save.set_visible(False)
    plt.savefig(newfilename, bbox_inches='tight')
    ax_save.set_visible(True)
    messagebox.showinfo('¡Listo!',
                        'Imagen guardada como {}.png'.format(
                os.path.split(os.path.splitext(newfilename)[0])[-1]))
check_save.on_clicked(check_save_callback)

# Once I have all that, I'll show the plot
plt.show()

# Like it is shown for the first time, autosave if configured that way
if autosave:
    newpath = os.path.join(path, 'Figuras')
    if not os.path.isdir(newpath):
        os.makedirs(newpath)
    ax_save.set_visible(False)
    plt.savefig(os.path.join(newpath, name+'_fit.png'), bbox_inches='tight')
    ax_save.set_visible(True)