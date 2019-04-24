# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: LEC
"""

from numpy import pi
import numpy as np
import matplotlib.pyplot as plt

#%%

# Get data
t, x = np.loadtxt('Datos.txt')
dt = 1
T = max(t) - min(t)

#%%

# Define general parameters
#cn = 8e-15
#noise = []
#for i in range(50):
#    noise.append(np.cos(2*pi*np.random.rand() * (t+dt) / (2*dt)))
#noise = np.array(noise)
#noise = cn * sum(noise)
#coherent_artifact = 1e-15 * np.exp(-t**2 / 2 / 0.07**2)
# Proporcional to cross-correlation

# Kind of answer we want:
# x = c1.*exp(-b1.*t).*cos(w1*t+p1) + c2.*exp(-b2.*t).*cos(w2.*t+p2) + 
#     + c3.*exp(-b3.*t.^2) + noise + coherent_artifact;

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
Nfit_terms = len(frequencies)
t2 = np.arange(0, N*dt, dt) # Time starting on zero
X2 = np.zeros((N, 2*Nfit_terms))
for i, b, omega in zip(range(Nfit_terms), damping_constants, frequencies):
    X2[:, 2*i] = np.exp(-b*t2) * np.cos(omega*t2)
    X2[:, 2*i+1] = -np.exp(-b*t2) * np.sin(omega*t2)

#%%

# Diagonalize square Hermitian modelled data matrix
[eigenvalues2, eigenvectors2] = np.linalg.eigh( np.matmul(X2, X2.T) )
ordered_index = eigenvalues2.argsort() # From smallest to largest absolute
eigenvalues2 = eigenvalues2[ordered_index] # Eigenvalues
eigenvectors2 = eigenvectors2[:, ordered_index] # Eigenvectors on columns
eigenvectors2 = np.array([l/l[0] for l in eigenvectors2.T]).T # Normalize
rank2 = np.linalg.matrix_rank(np.diag(eigenvalues2)) # Size measure

#%%

# Choose number of significant values
Nsignificant2 = np.linalg.matrix_rank( np.matmul(X2, X2.T) )

# Crop data according to it
F2 = np.zeros((N, N))
F2[-Nsignificant2:,-Nsignificant2:] = np.diag(
        1/np.sqrt(eigenvalues2[-Nsignificant2:]))
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

#%%

# -----------------------------------------------------------------------------
# SOLUTION, PLOTS AND STATISTICS
# -----------------------------------------------------------------------------

# Solution
fit_terms = np.array([a * np.exp(-b*t) * np.cos(omega*t + phi)
                     for a, b, omega, phi in zip(amplitudes,
                                                 damping_constants,
                                                 frequencies,
                                                 phases)]).T
fit = sum(fit_terms.T)
square_chi = sum( (fit - x)**2 ) / N # Best if absolute is smaller

#%%

# Statistics of the residue
residue = x - fit
residue_transform = np.fft.rfft(residue)
residue_frequencies = 1000 * np.fft.rfftfreq(N, d=dt) # in GHz
plt.plot(residue_frequencies, residue_transform)

"""
Y1=fft(residue);
NN=length(Y1);
Pyy= Y1.*conj(Y1) / (NN-1);
frequ = 1/(NN*(t(2)-t(1))) * (0:((NN/2)-1)) * 1000; %in GHz
FFTResidue=Pyy(1:NN/2);
plot(frequ, Pyy(1:NN/2))
"""

# I BELIEVE THE ORIGINAL FFT WAS FAULTY!
# OURS GOES UP TO 500 AND THE OLD ONE GOES TO 250.

#%%

# Raman-like Spectrum parameters
f_damping_constants = 1000 * damping_constants / (2*pi) #in GHz
f_frequencies = 1000 * frequencies / (2*pi) #in GHz
f_max = max(f_frequencies)
f_independent = np.arange(0, 1.5*f_max, f_max/1000)

# Raman-like Spectrum per se
response = np.zeros( (len(f_independent), Nfit_terms) )
for i in range(Nfit_terms):
   if frequencies[i]==0:
      response[:,i] = 0
   else:
      response[:,i] = amplitudes[i] * np.imag( phases[i] / 
         (f_frequencies[i]**2 - f_independent**2 - 
          2j * f_independent * f_damping_constants[i]))
spectrum = np.sum(response, axis=1)

plt.plot(f_independent, spectrum)

"""
b = B * 1000 / (2*pi); %in GHz
f = W * 1000 / (2*pi); %in GHz
Wmax = max(f);
freq = 0 : fMAX/1000 : 1.5*fMAX; %in GHz
freq=freq';

res = zeros(length(freq), length(W));

for i=1:length(W)
   if W(i)==0
      res(:,i)=0;
   else 
      %res(:,i)=C(i)*imag(1i*exp(-1i*fi(i))*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
      res(:,i)=C(i)*imag(1*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
  end
  end
  
  spectrum=sum(res,2);
%figure;
subplot(3,2,6);
plot(freq,spectrum,freq,res);title('Raman-like Data Spectrum');xlabel('freq (GHz)');ylabel('spectrum');
"""

# THIS PART DOES DEFINITELY NOT WORK >.<