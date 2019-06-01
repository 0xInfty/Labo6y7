# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_plot_module as ivp
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
from numpy import pi
import numpy as np

#%%

def cropData(t0, t, *args, **kwargs):
    
    """Crops data according to a certain logic specifier.
    
    By default, the logic specifier is '>=', so data 't' is cropped from 
    certain value 't0'. It is flexible though. For example, if a parameter 
    'logic="<="' is delivered as argument to this function, then data 't' will 
    be cropped up to certain value 't0'.
    
    Parameters
    ----------
    t0 : int, float
        Value to apply logic specifier to.
    t : np.array
        Data to apply logic specifier to.
    logic='>=' : str
        Logic specifier.
    
    Returns
    -------
    new_args : np.array
        Resultant data.
    
    """
    
    try:
        logic = kwargs['logic']
    except:
        logic = '>='
    
    if t0 not in t:
        raise ValueError("Hey! t0 must be in t")
    
    index = eval("t{}t0".format(logic))
    new_args = []
    for a in args:
        try:
            a = np.array(a)
        except:
            raise TypeError("Extra arguments must be array-like")
        if a.ndim == 1:
            new_args.append(a[index])
        else:
            try:
                new_args.append(a[index, :])
            except:
                raise ValueError("This function takes only 1D or 2D arrays")
    t = t[index]
    new_args = [t, *new_args]
    
    return new_args

#%%

def roundMatlab(x):
    
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
    
    Returns
    -------
    y : int
        Rounded number.
    """
    
    isRoundMatlabNeeded = round(80.5) == 81
    
    if isRoundMatlabNeeded:
        
        xround = int(x)
        even = xround/2 == int(xround/2) # True if multiple of 2
        
        if even:
            y = round(x) + 1
        else:
            y = round(x)
            
        return int(y)
    
    else:
    
        return int(round(x))

#%%

def linearPrediction(t, x, autoclose=True):
    
    """Applies linear prediction fit to data.
    
    Given a set of data :math:`t, x` with independent step :math:`dt`, it looks 
    for the best fit according to the model...
    
    .. math:: f(t) = \sum_i A.cos(\omega_i t + \phi).e^{-\frac{t}{\tau}}
    
    This method does not need initial values for the parameters to fit.
    
    In order for it to work, this function creates a uniform independent variable 
    whose elements are multiples :math:`t_i=i.dt` of a constant step :math:`dt`.
    
    Parameters
    ----------
    t : np.array
        Independent variable :math:`t`.
    x : np.array
        Dependent variable :math:`x`.
    autoclose=True : bool
        Says whether to close the intermediate eigenvalues' plot or not.
    
    Returns
    -------
    results : np.array
        Parameters that best fit the data. On its columns it holds...
        ...frequency :math:`f=2\pi\omega` in :math:`\frac{1}{t}` units.
        ...characteristic time :math:`\tau_i` in :math:`t` units.
        ...quality factors :math:`Q_i=\frac{\omega}{2\gamma}=\pi f \tau`
        ...amplitudes :math:`A_i` in the same units as :math:`x`
        ...phases :math:`\phi_i` written in radians in multiples of :math:`\pi`
    parameters : ivu.InstancesDict
        Other fit parameters...
        ...chi squared :math:`\chi^2`
        ...number of significant values :math:`N`
        ...time range :math:`t_0,t_f`
    
    See also
    --------
    ivp.linearPredictionPlot
    
    """
    
    # -----------------------------------------------------------------------------
    # FREQUENCIES AND DAMPING FACTORS
    # -----------------------------------------------------------------------------
    
    # Create uniform time
    N = len(t)
    dt = np.mean(np.diff(t))
    t = np.arange(t[0], t[0] + N*dt, dt)
    
    # Create data matrix
    M = roundMatlab(0.75 * N)
    X = np.array([x[i+j+1] for j in range(N-M) for i in range(M)]).reshape((N-M,M))
    
    # Diagonalize square data matrix
    [eigenvalues, eigenvectors] = np.linalg.eigh( np.matmul(X, X.T) )
    ordered_index = eigenvalues.argsort() # From smallest to largest value 
    eigenvalues = eigenvalues[ordered_index] # Eigenvalues
    eigenvectors = eigenvectors[:, ordered_index] # Eigenvectors on columns
    #eigenvectors = np.array([l/l[0] for l in eigenvectors.T]).T # Normalize
    rank = np.linalg.matrix_rank(np.diag(eigenvalues)) # Size measure
    
    # Choose number of significant values
    Nsignificant = 4
    fig = plt.figure()
    ax = plt.subplot()
    plt.semilogy(eigenvalues, linestyle='none', marker='o', 
                 fillstyle='none', markersize=10)
    plt.title('¿Número de valores singulares?')
    plt.ylabel("Autovalores")
    Nsignificant = ivp.interactiveIntegerSelector(ax, min_value=0, max_value=8)
    if autoclose:
        plt.close(fig)
    
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
    angular_frequencies = angular_frequencies[ordered_index] # in 2pi/dt units
    damping_constants = damping_constants[ordered_index] # in 1/dt units
    
    # Crop them according to number of real roots and rank of diagonalized matrix
    Nzeros = len(angular_frequencies) - np.count_nonzero(angular_frequencies)
    angular_frequencies = abs(angular_frequencies)[:roundMatlab(
            (rank-Nzeros)/2+Nzeros)]
    damping_constants = damping_constants[:roundMatlab(
            (rank-Nzeros)/2+Nzeros)]
    
    # Then crop them according to the number of positive or zero damping constants
    Npositives = len(damping_constants[damping_constants>=0])
    ordered_index = damping_constants.argsort()[::-1] # From largest to smallest
    damping_constants = damping_constants[ordered_index][:Npositives]
    angular_frequencies = angular_frequencies[ordered_index][:Npositives]
    
    # Now I have the smallest frequencies and largest damping constants
    # Then I calculate the characteristic time tau and the quality factor Q
    quality_factors = angular_frequencies / (2*damping_constants)
    characteristic_times = 1 / damping_constants # in dt units
    
    # -----------------------------------------------------------------------------
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
    frequencies = angular_frequencies / (2*pi) # in 1/dt units
    amplitudes = np.array(amplitudes)
    phases = np.array(phases)
    pi_phases = phases / pi # in radians written as multiples of pi
    if Nfit_terms==0:
        raise ValueError("¡Error! No se encontraron términos de ajuste")
    elif Nfit_terms>1:
        print("¡Listo! Encontramos {} términos".format(Nfit_terms))
    else:
        print("¡Listo! Encontramos {} término".format(Nfit_terms))
    print("Frecuencias: {}".format(frequencies))
    
    # -----------------------------------------------------------------------------
    # FIT AND RESULTS
    # -----------------------------------------------------------------------------
    
    fit_terms = np.array([a * np.exp(-b*(t-t[0])) * np.cos(omega*(t-t[0]) + phi)
                         for a, b, omega, phi in zip(amplitudes,
                                                     damping_constants,
                                                     angular_frequencies,
                                                     phases)]).T
    fit = sum(fit_terms.T)
    chi_squared = sum( (fit - x)**2 ) / N # Best if absolute is smaller
    print("Chi cuadrado \u03C7\u00B2: {:.2e}".format(chi_squared))
                     
    ## Statistics of the residue
    #residue = x - fit
    #residue_transform = abs(np.fft.rfft(residue))
    #residue_frequencies = 1000 * np.fft.rfftfreq(N, d=dt) # in GHz
    #plt.plot(residue_frequencies, residue_transform)
    
    # What I would like this function to return
    results = np.array([frequencies,
                        characteristic_times,
                        quality_factors,
                        amplitudes,
                        pi_phases]).T

    # Some important parameters I would also like to return
    parameters = dict(chi_squared = chi_squared,
                      Nsingular_values = Nsignificant,
                      t0 = t[0], tf = t[-1])
    parameters = ivu.InstancesDict(parameters)
    
    return results, parameters

#%%

def makeRamanLikeSpectrum(frequencies, characteristic_times, amplitudes):
    
    """Makes a Raman-like spectrum for one or more damped sinusoids"""
    
    # Raman-like Spectrum data organization
    try:
        Nfit_terms = len(frequencies)
    except:
        Nfit_terms = 1
        frequencies = np.array([frequencies])
        characteristic_times = np.array([characteristic_times])
        amplitudes = np.array([amplitudes])
    
    # Raman-like Spectrum parameters
    max_frequency = max(frequencies)
    frequencies_damping = 1 / (2*pi*characteristic_times)
    if max_frequency != 0:
        raman_frequencies = np.arange(0, 1.5*max_frequency, max_frequency/1000)
    else:
        raman_frequencies = np.array([0, 12])
        
    # Raman-like Spectrum per se
    raman_spectrum_terms = np.zeros((len(raman_frequencies), Nfit_terms))
    for i in range(Nfit_terms):
       if frequencies[i]==0:
          raman_spectrum_terms[:,i] = 0
       else:
          raman_spectrum_terms[:,i] = amplitudes[i] * np.imag( 
             frequencies[i] / 
             (frequencies[i]**2 - raman_frequencies**2 - 
             2j * raman_frequencies * frequencies_damping[i]))
    raman_spectrum = np.sum(raman_spectrum_terms, axis=1)

#    raman = np.array([raman_frequencies, raman_spectrum,
#                      *list(raman_spectrum_terms.T)]).T
    
    return raman_frequencies, raman_spectrum, raman_spectrum_terms

#%% PMUSIC -----------------------------------------------------------------

## Start by defining general parameters
#PMN = nsize//4 # WHY /4? "cantidad de mediciones"
#PMT = 1200 # size of the window in ps
#PMdt = 20 # time step in ps
#
## Now get PMUSIC's parameters3
#PMdata = detrend(meanV) # BEWARE! THE BEST FIT HORIZONTAL LINE IS FILTERED!
#Mp = [PMN, 200] # This is PMUSIC's most important parameter
## Mp = [components' dimension, harmonics' limit]
## Second variable marks how many harmonics to throw away.
## It can't be greater than the measurement's dimension.
#
## Then define several variables to be filled
#MSn = []
#Mfn = []
#iPMindex=0
#for i in range(PMT+1, 1350, PMdt): # WHY? SHOULD I BE INCLUDING 1350? I'M NOT.
#
#    iPMindex = iPMindex + 1
#
#    # Take a segment of data and apply PMUSIC
#    iPMdata = PMdata[((i-PMT) < t) & (t < i)]
#    [MSn1, Mfn1] = pmusic(iPMdata, Mp, 6000, samplerate, [], 0)
#    # WHY 6000?
#        
#    iPMselection = ((Mfn1 >= 0) & (Mfn1 <= 0.06));
#    MSn[:, iPMindex] = MSn1[iPMselection]
#    Mfn[:, iPMindex] = Mfn1[iPMselection]
#
## Finally... Plot! :)
#plt.figure(1)
#plt.subplot(3,2,1)
#plt.imagesc(np.arange(1,T), Mfn[:,1], MSn)

"""
Problems so far:
    Don't have a pmusic Python equivalent
    Haven't searched an imagesc equivalent
"""

#%%
    
def linearPredictionPlot(t, x, results, showgrid=False):

    """Plots the results of a linear prediction fit.
    
    Parameters
    ----------
    t : np.array
        Independent variable :math:`t`.
    x : np.array
        Dependent variable :math:`x`.
    results : np.array
        Parameters that best fit the data. On its columns it holds...
        ...frequency :math:`f=2\pi\omega` in :math:`\frac{1}{t}` units.
        ...characteristic time :math:`\tau_i` in :math:`t` units.
        ...quality factors :math:`Q_i=\frac{\omega}{2\gamma}=\pi f \tau`
        ...amplitudes :math:`A_i` in the same units as :math:`x`
        ...phases :math:`\phi_i` written in radians in multiples of :math:`\pi`
    showgrid=False : bool
        Says whether to show or not the vertical grid on the time space plot.
    
    Returns
    -------
    fig : plt.Figure instance
        Figure containing the desired plot.
    legend_buttons : wid.CheckButtons
        Interactive legend. Only returned if making an interactive plot.
    
    See also
    --------
    iva.linearPrediction
    
    """
    
    # BIG PROBLEM! I NEED THE TIME!
    
    # First I deglose data
    Nfit_terms = results.shape[0]
    fit_terms = np.array([a * np.exp(-b*(t-t[0])) * np.cos(omega*(t-t[0]) + phi)
                         for a, b, q, omega, phi in results.T])
    fit_total = sum(fit_terms.T)
    fit = np.array([t, fit_total, *fit_terms]).T
    raman = np.array(makeRamanLikeSpectrum).T
    
    # Then, to plot, I first start a figure
    fig = plt.figure()
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
    ldata, = plt.plot(fit[:,0], fit[:,1], 'k', linewidth=0.4)
    ax_data.autoscale(False)
    lfit, = plt.plot(fit[:,0], fit[:,2], linewidth=2)
    lfit_terms = plt.plot(fit[:,0], fit[:,3:], linewidth=2)
    for l in lfit_terms: l.set_visible(False)
    plt.xlabel("Tiempo (ps)")
    plt.ylabel(r"Voltaje ($\mu$V)")
    ax_data.tick_params(labelsize=12)
    if showgrid:
        ax_data.minorticks_on()
        ax_data.tick_params(axis='y', which='minor', left=False)
        ax_data.tick_params(length=5)
        ax_data.grid(axis='x', which='both')
        ldata.set_linewidth(0.6)
        lfit.set_linewidth(2.3)
    
    # Because it's pretty, I make an interactive legend
    ax_legend = plt.axes([0.75, 0.642, 0.155, 0.24])
    legend_buttons = wid.CheckButtons(ax_legend, ('Data', 
                                     'Ajuste', 
                                     *['Término {:.0f}'.format(i+1) 
                                     for i in range(Nfit_terms)]), 
                        (True, True, *[False for i in range(Nfit_terms)]))
    legend_buttons.labels[1].set_color(lfit.get_color())
    for leg, l in zip(legend_buttons.labels[2:], lfit_terms):
        leg.set_color(l.get_color())
    
    # For that, I'll need a callback function  
    def legend_callback(label):
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
    legend_buttons.on_clicked(legend_callback)
    
#    # Since I can, I would also like an interactive 'Save' button
#    save_button = ivp.interactiveSaveButton(filename, sufix='_fit')
#    
#    # Once I have all that, I'll show the plot
#    plt.show()
#    
#    # Like it is shown for the first time, autosave if configured that way
#    if autosave:
#        newpath = os.path.join(path, 'Figuras')
#        if not os.path.isdir(newpath):
#            os.makedirs(newpath)
#        save_button.ax.set_visible(False)
#        plt.savefig(os.path.join(newpath, name+'_fit.png'), bbox_inches='tight')
#        save_button.ax.set_visible(True)
        
    return fig, legend_buttons

#%%

def linearPredictionTables(parameters, results, other_results, units=["Hz","s"]):

    terms_heading = ["F ({})".format(units[0]), "\u03C4 ({})".format(units[1]), 
                     "Q", "A (u.a.)", "Fase (\u03C0rad)"]
    terms_heading = '\t'.join(terms_heading)
    terms_table = ['\t'.join([str(element) for element in row]) for row in results]
    terms_table = '\n'.join(terms_table)
    terms_table = '\n'.join([terms_heading, terms_table])
    
    fit_heading = ["Experimentos utilizados",
                   "Número de valores singulares",
                   "Porcentaje enviado a cero (%)",
                   "Método de corrimiento",
                   "Corrimiento V\u2080 (\u03BCV)",               
                   r"Rango temporal → Inicio (ps)",
                   r"Rango temporal → Final (ps)",
                   "Chi cuadrado \u03C7\u00B2"]
    
    if parameters.use_full_mean:
        used_experiments = 'Todos'
    else:
        used_experiments = ', '.join([str('{:.0f}'.format(i+1)) 
                                      for i in parameters.use_experiments])
        if len(parameters.use_experiments)==1:
            used_experiments = 'Sólo ' + used_experiments
        else:
            used_experiments = 'Sólo ' + used_experiments
    if parameters.send_tail_to_zero:
        tail_percent = parameters.use_fraction*100
    else:
        tail_percent = 0
    if parameters.tail_method=='mean':
        method = 'Promedio'
    elif parameters.tail_method=='min':
        method = 'Mínimo'
    elif parameters.tail_method=='max':
        method = 'Máximo'
    else:
        method = 'Desconocido'
    
    fit = [used_experiments,
           str(other_results['Nsingular_values']),
           '{:.0f}'.format(tail_percent),
           method,
           str(parameters.voltage_zero),
           str(parameters.time_range[0]),
           str(parameters.time_range[1]),
           '{:.2e}'.format(other_results['chi_squared'])]
    fit_table = ['\t'.join([h, f]) for h, f in zip(fit_heading, fit)]
    fit_table = '\n'.join(fit_table)
    
    return terms_table, fit_table

#%%

#def arrayTable(array, heading_list=None, axis=0):
#    
#    if heading_list is not None:
#        heading = '\t'.join(heading_list)
#    if axis==1:
#        array = array.T
#    elif axis!=0:
#        raise ValueError("Axis must be 0 or 1!")
#    items = ['\t'.join([str(element) for element in row]) for row in array]
#    items = '\n'.join(items)
#    if heading_list is not None:
#        table = '\n'.join([heading, items])
#    else:
#        table = items
#    
#    return table