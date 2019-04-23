# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: LEC
"""

from numpy import pi
import numpy as np
#import matplotlib.pyplot as plt

#%%

# Get data
t, x = np.loadtxt('Datos.txt')
dt = 1

#%%

# Define general parameters
cn = 8e-15
noise = []
for i in range(50):
    noise.append(np.cos(2*pi*np.random.rand() * (t+dt) / (2*dt)))
noise = np.array(noise)
noise = cn * sum(noise)
coherent_artifact = 1e-15 * np.exp(-t**2 / 2 / 0.07**2)
# Proporcional to cross-correlation

"""
Kind of answer we want:
x = c1.*exp(-b1.*t).*cos(w1*t+p1) + c2.*exp(-b2.*t).*cos(w2.*t+p2) + 
    + c3.*exp(-b3.*t.^2) + noise + coherent_artifact;
"""

#%%

# -----------------------------------------------------------------------------
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
eigenvectors = np.array([l/l[0] for l in eigenvectors.T]).T # Normalize
rank = np.linalg.matrix_rank(np.diag(eigenvalues)) # Size measure

"""
La segunda columna de eV la tengo con signo invertido respecto a Matlab.
==> Obtengo resultados distintos para U.
"""

#%%

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

#%%

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

#%%

# Calculate damping constants 'b' and frequencies 'omega'
damping_constants = (np.log(abs(roots)) / dt)[:rank] # Crop them accordingly
frequencies = (np.angle(roots) / dt)[:rank]

#%%

# Sort them
ordered_index = frequencies.argsort() # From smallest to largest frequency
frequencies = frequencies[ordered_index]
damping_constants = damping_constants[ordered_index]

# Crop them according to number of real roots and rank of diagonalized matrix
Nzeros = len(frequencies) - np.count_nonzero(frequencies)
frequencies = abs(frequencies)[:int(round((rank-Nzeros)/2+Nzeros))]
damping_constants = damping_constants[:int(round((rank-Nzeros)/2+Nzeros))]

# Then crop them according to the number of positive or zero damping constants
Npositives = len(damping_constants[damping_constants>=0])
ordered_index = damping_constants.argsort()[::-1] # From largest to smallest
damping_constants = damping_constants[ordered_index][:Npositives]
frequencies = frequencies[ordered_index][:Npositives]

# Now I have the smallest frequencies and largest damping constants

#%%

# -----------------------------------------------------------------------------
# AMPLITUDES AND PHASES
# -----------------------------------------------------------------------------

# Create modelled data matrix
Nsolutions = len(frequencies)
t2 = np.arange(0, N*dt, dt) # Time starting on zero
X2 = np.zeros((N, 2*Nsolutions))
for i, b, omega in zip(range(Nsolutions), damping_constants, frequencies):
    X2[:, 2*i] = np.exp(-b*t2) * np.cos(omega*t2)
    X2[:, 2*i+1] = -np.exp(-b*t2) * np.sin(omega*t2)

#%%

# Diagonalize square modelled data matrix
[eigenvalues2, eigenvectors2] = np.linalg.eig( np.matmul(X2, X2.T) )
ordered_index = abs(eigenvalues2).argsort() # From smallest to largest absolute
eigenvalues2 = eigenvalues2[ordered_index] # Eigenvalues
eigenvectors2 = eigenvectors2[:, ordered_index] # Eigenvectors on columns
eigenvectors2 = np.array([l/l[0] for l in eigenvectors2.T]).T # Normalize
rank2 = np.linalg.matrix_rank(np.diag(eigenvalues2)) # Size measure

# $&@%!! COMPLEX EIGENVALUES AND EIGENVECTORS >.<
# WE REALLY, REALLY HAVE TO LOOK FOR ANOTHER METHOD!

#%%

# Choose number of significant values
Nsignificant2 = np.linalg.matrix_rank( np.matmul(X2, X2.T) )

# Crop data according to it
F2 = np.zeros((N, N))
F2[-Nsignificant2:,-Nsignificant2:] = np.diag(
        1/np.sqrt(abs(eigenvalues2[-Nsignificant2:])))
auxiliar = np.matmul(eigenvectors2, F2)
U2 = np.matmul(X2.T, auxiliar) # Xmatrix.T * eigenvectors * F

#%%

# Get defining vector
auxiliar = np.matmul(eigenvectors2.T, x)
auxiliar = np.matmul(F2.T, auxiliar)
A2 = np.matmul(U2, auxiliar) # U * F.T * eigenvectors.T * xvector 
# |--> Least-Squares?

#%%

# Calculate phases 'phi' and amplitudes 'C'
amplitudes = []
phases = []
for i in range(Nsolutions):
    
    if A2[2*i]==0 and A2[2*i+1]==0:
        amplitudes[i] == 0
        phases[i] == 0
    elif A2[2*i]==0:
        amplitudes[i] = abs(A2[2*i+1])
        phases[i] = np.sign(A2[2*i+1]) * pi/2
    elif A2[2*i+1]==0:
        amplitudes[i] = abs(A2[2*i])
        phases[i] = (1-np.sign(A2[2*i])) * pi/2
    else:
        amplitudes[i] = np.sqrt(A2[2*i+1]**2 + A2(2*i)**2)
        phases[i] = np.arctan2(A2[2*i+1], A2(2*i))

#%%

# -----------------------------------------------------------------------------
# SOLUTION, PLOTS AND STATISTICS
# -----------------------------------------------------------------------------

# Solution
solutions = np.array([a * np.exp(-b*t) * np.cos(omega*t + fi)
                      for a, b, omega, fi in zip(amplitudes,
                                                 damping_constants,
                                                 frequencies,
                                                 phases)]).T
total_solution = sum(solutions.T)

# Statistics
square_chi = sum( (total_solution - x)**2 ) / N
residue = x - total_solution

# MISSING STATISTICS AND PLOTS BUT THAT'S IT :)