# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: LEC
"""

import ivs_plot_module as ivp
import ivs_save_module as ivs
import matplotlib.pyplot as plt
from numpy import pi
import numpy as np

#%%

def roundMatlab(x, round_Matlab_needed=True):
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
    round_Matlab_needed=True : bool
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
    
    if round_Matlab_needed:
        return y
    else:
        return round(x)
    
#%%

def linearPrediction(t, x, dt, autoclose=True, round_Matlab_needed=True):
    
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
    
    # Sort them
    ordered_index = angular_frequencies.argsort() # From smallest to largest freq
    angular_frequencies = angular_frequencies[ordered_index]
    damping_constants = damping_constants[ordered_index]
    
    # Crop them according to number of real roots and rank of diagonalized matrix
    Nzeros = len(angular_frequencies) - np.count_nonzero(angular_frequencies)
    angular_frequencies = abs(angular_frequencies)[:int(roundMatlab(
            (rank-Nzeros)/2+Nzeros,
            round_Matlab_needed))]
    damping_constants = damping_constants[:int(roundMatlab(
            (rank-Nzeros)/2+Nzeros,
            round_Matlab_needed))]
    
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
    
    return results, others