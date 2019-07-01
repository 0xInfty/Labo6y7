# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: Vall
"""

import iv_plot_module as ivp
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
from math import sqrt
from numpy import pi
import numpy as np
from scipy.optimize import curve_fit

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

def linearFit(X, Y, dY=None, showplot=True,
              plot_some_errors=(False, 20), **kwargs):
    
    """Applies linear fit and returns m, b and Rsq. Can also plot it.
    
    By default, it applies minimum-square linear fit 'y = m*x + b'. If 
    dY is specified, it applies weighted minimum-square linear fit.
    
    Parameters
    ----------
    X : np.array, list
        Independent X data to fit.
    Y : np-array, list
        Dependent Y data to fit.
    dY : np-array, list
        Dependent Y data's associated error.
    shoplot : bool
        Says whether to plot or not.
    plot_some_errors : tuple (bool, int)
        Says wehther to plot only some error bars (bool) and specifies 
        the number of error bars to plot.
    
    Returns
    -------
    rsq : float
        Linear fit's R Square Coefficient.
    (m, dm): tuple (float, float)
        Linear fit's slope: value and associated error, both as float.
    (b, db): tuple (float, float)
        Linear fit's origin: value and associated error, both as float.
        
    Other Parameters
    ----------------
    txt_position : tuple (horizontal, vertical), optional
        Indicates the parameters' text position. Each of its values 
        should be a number (distance in points measured on figure). 
        But vertical value can also be 'up' or 'down'.
    mb_units : tuple (m_units, b_units), optional
        Indicates the parameter's units. Each of its values should be a 
        string.
    mb_error_digits : tuple (m_error_digits, b_error_digits), optional
        Indicates the number of digits to print in the parameters' 
        associated error. Default is 3 for slope 'm' and 2 for intercept 
        'b'.
    mb_string_scale : tuple (m_string_scale, b_string_scale), optional
        Indicates whether to apply string prefix's scale to printed 
        parameters. Each of its values should be a bool; i.e.: 'True' 
        means 'm=1230.03 V' with 'dm = 10.32 V' would be printed as 
        'm = (1.230 + 0.010) V'. Default is '(False, False)'.
    rsq_decimal_digits : int, optional.
        Indicates the number of digits to print in the Rsq. Default: 3.
        
    Warnings
    --------
    The returned Rsq doesn't take dY weights into account.
    
    """

    # ¿Cómo hago Rsq si tengo pesos?
    
    if dY is None:
        W = None
    else:
        W = 1/dY**2
        
    fit_data = np.polyfit(X, Y, 1, cov=True, w=W)
    
    m = fit_data[0][0]
    dm = sqrt(fit_data[1][0,0])
    b = fit_data[0][1]
    db = sqrt(fit_data[1][1,1])
    rsq = 1 - sum( (Y - m*X - b)**2 )/sum( (Y - np.mean(Y))**2 )

    try:
        kwargs['text_position']
    except KeyError:
        if m > 1:
            aux = 'up'
        else:
            aux = 'down'
        kwargs['text_position'] = (.02, aux)

    if showplot:

        plt.figure()
        if dY is None:
            plt.plot(X, Y, 'b.', zorder=0)
        else:
            if plot_some_errors[0] == False:
                plt.errorbar(X, Y, yerr=dY, linestyle='', marker='.',
                             ecolor='b', elinewidth=1.5, zorder=0)
            else:
                plt.errorbar(X, Y, yerr=dY, linestyle='', marker='.',
                             color='b', ecolor='b', elinewidth=1.5,
                             errorevery=len(Y)/plot_some_errors[1], 
                             zorder=0)
        plt.plot(X, m*X+b, 'r-', zorder=100)
        plt.legend(["Ajuste lineal ponderado","Datos"])
        
        kwargs_list = ['mb_units', 'mb_string_scale', 
                       'mb_error_digits', 'rsq_decimal_digits']
        kwargs_default = [('', ''), (False, False), (3, 2), 3]
        for key, value in zip(kwargs_list, kwargs_default):
            try:
                kwargs[key]
            except KeyError:
                kwargs[key] = value
        
        if kwargs['text_position'][1] == 'up':
            vertical = [.9, .82, .74]
        elif kwargs['text_position'][1] == 'down':
            vertical = [.05, .13, .21]
        else:
            if kwargs['text_position'][1] <= .08:
                fact = .08
            else:
                fact = -.08
            vertical = [kwargs['text_position'][1]+fact*i for i in range(3)]
        

        plt.annotate('m = {}'.format(ivu.errorValueLatex(
                        m, 
                        dm,
                        error_digits=kwargs['mb_error_digits'][0],
                        units=kwargs['mb_units'][0],
                        string_scale=kwargs['mb_string_scale'][0],
                        one_point_scale=True)),
                    (kwargs['text_position'][0], vertical[0]),
                    xycoords='axes fraction')
        plt.annotate('b = {}'.format(ivu.errorValueLatex(
                        b, 
                        db,
                        error_digits=kwargs['mb_error_digits'][1],
                        units=kwargs['mb_units'][1],
                        string_scale=kwargs['mb_string_scale'][1],
                        one_point_scale=True)),
                    (kwargs['text_position'][0], vertical[1]),
                    xycoords='axes fraction')
        rsqft = r'$R^2$ = {:.' + str(kwargs['rsq_decimal_digits']) + 'f}'
        plt.annotate(rsqft.format(rsq),
                    (kwargs['text_position'][0], vertical[2]),
                    xycoords='axes fraction')

        
        plt.show()

    return rsq, (m, dm), (b, db)

#%%

def nonLinearFit(X, Y, fitfunction, initial_guess=None, dY=None, 
                 showplot=True, plot_some_errors=(False, 20), 
                 **kwargs):
    """Applies nonlinear fit and returns parameters and Rsq. Plots it.
    
    By default, it applies minimum-square fit. If dY is specified, it 
    applies weighted minimum-square fit.
    
    Parameters
    ----------
    X : np.array, list
        Independent X data to fit.
    Y : np-array, list
        Dependent Y data to fit.
    fitfunction : function
        The function you want to apply. Its arguments must be 'X' as 
        np.array followed by the other parameters 'a0', 'a1', etc as 
        float. Must return only 'Y' as np.array.
    initial_guess=None : list, optional
        A list containing a initial guess for each parameter.
    dY : np-array, list, optional
        Dependent Y data's associated error.
    shoplot : bool
        Says whether to plot or not.
    plot_some_errors : tuple (bool, int)
        Says wehther to plot only some error bars (bool) and specifies 
        the number of error bars to plot.
    
    Returns
    -------
    rsq : float
        Fit's R Square Coefficient.
    parameters : list of tuples
        Fit's parameters, each as a tuple containing value and error, 
        both as tuples.
    
    Other Parameters
    -----------------
    txt_position : tuple (horizontal, vertical), optional
        Indicates the parameters' text position. Each of its values 
        should be a number (distance in points measured on figure). 
        But vertical value can also be 'up' or 'down'.
    par_units : list, optional
        Indicates the parameters' units. Each of its values should be a 
        string.
    par_error_digits : list, optional
        Indicates the number of digits to print in the parameters' 
        associated error. Default is 3 for all of them.
    par_string_scale : list, optional
        Indicates whether to apply string prefix's scale to printed 
        parameters. Each of its values should be a bool. Default is 
        False for all of them.
    rsq_decimal_digits : int, optional.
        Indicates the number of digits to print in the Rsq. Default: 3.
        
    Warnings
    --------
    The returned Rsq doesn't take dY weights into account.
    
    """
    
    if not isinstance(X, np.ndarray):
        raise TypeError("X should be a np.array")
    if not isinstance(Y, np.ndarray):
        raise TypeError("Y should be a np.array")
    if not isinstance(dY, np.ndarray) and dY is not None:
        raise TypeError("dY shouuld be a np.array")
    if len(X) != len(Y):
        raise IndexError("X and Y must have same lenght")
    if dY is not None and len(dY) != len(Y):
        raise IndexError("dY and Y must have same lenght")
    
    if dY is None:
        W = None
    else:
        W = 1/dY**2
    
    parameters, covariance = curve_fit(fitfunction, X, Y,
                                       p0=initial_guess, sigma=W)  
    rsq = sum( (Y - fitfunction(X, *parameters))**2 )
    rsq = rsq/sum( (Y - np.mean(Y))**2 )
    rsq = 1 - rsq

    if showplot:
        
        plt.figure()
        if dY is None:
            plt.plot(X, Y, 'b.', zorder=0)
        else:
            if plot_some_errors[0] == False:
                plt.errorbar(X, Y, yerr=dY, linestyle='b', marker='.',
                             ecolor='b', elinewidth=1.5, zorder=0)
            else:
                plt.errorbar(X, Y, yerr=dY, linestyle='-', marker='.',
                             color='b', ecolor='b', elinewidth=1.5,
                             errorevery=len(Y)/plot_some_errors[1], 
                             zorder=0)
        plt.plot(X, fitfunction(X, *parameters), 'r-', zorder=100)        
        plt.legend(["Ajuste lineal ponderado","Datos"])
        
        n = len(parameters)
        kwargs_list = ['text_position', 'par_units', 'par_string_scale', 
                       'par_error_digits', 'rsq_decimal_digits']
        kwargs_default = [(.02,'up'), ['' for i in range(n)], 
                          [False for i in range(n)], 
                          [3 for i in range(n)], 3]
        for key, value in zip(kwargs_list, kwargs_default):
            try:
                kwargs[key]
                if key != 'text_position':
                    try:
                        if len(kwargs[key]) != n:
                            print("Wrong number of parameters",
                                  "on '{}'".format(key))
                            kwargs[key] = value
                    except TypeError:
                        kwargs[key] = [kwargs[key] for i in len(n)]
            except KeyError:
                kwargs[key] = value
        
        if kwargs['text_position'][1] == 'up':
            vertical = [.9-i*.08 for i in range(n+1)]
        elif kwargs['text_position'][1] == 'down':
            vertical = [.05+i*.08 for i in range(n+1)]
        else:
            if kwargs['text_position'][1] <= .08:
                fact = .08
            else:
                fact = -.08
            vertical = [
                kwargs['text_position'][1]+fact*i for i in range(n+1)]
        
        for i in range(n):
            plt.annotate(
                    'a{} = {}'.format(
                        i,
                        ivu.errorValueLatex(
                            parameters[i], 
                            sqrt(covariance[i,i]),
                            error_digits=kwargs['par_error_digits'][i],
                            units=kwargs['par_units'][i],
                            string_scale=kwargs['par_string_scale'][i],
                            one_point_scale=True)),
                    (kwargs['text_position'][0], vertical[i]),
                    xycoords='axes fraction')
        rsqft = r'$R^2$ = {:.' + str(kwargs['rsq_decimal_digits'])+'f}'
        plt.annotate(rsqft.format(rsq),
                    (kwargs['text_position'][0], vertical[-i]),
                    xycoords='axes fraction')
        
        plt.show()
    
    parameters_error = np.array(
            [sqrt(covariance[i,i]) for i in range(n)])
    parameters = list(zip(parameters, parameters_error))
    
    return rsq, parameters

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

def linearPrediction(t, x, dt, svalues=None, max_svalues=8, autoclose=True):
    
    """Applies linear prediction fit to data.
    
    Given a set of data :math:`t, x` with independent step :math:`dt`, it looks 
    for the best fit according to the model...
    
    .. math:: f(t) = \sum_i A.cos(\omega_i t + \phi).e^{-\frac{t}{\tau}}
    
    This method does not need initial values for the parameters to fit.
    
    In order for it to work, it is necesary though to have a uniform 
    independent variable whose elements are multiples :math:`t_i=i.dt` of a 
    constant step :math:`dt`
    
    Parameters
    ----------
    t : np.array
        Independent variable :math:`t` in ps.
    x : np.array
        Dependent variable :math:`x` in any unit.
    dt : float
        Independent variable's step :math:`dt` in ps.
    svalues : None, int
        Number of significant values. If set to None, it's chosen in an 
        interactive way.
    max_svalues : int
        Maximum number of significant values that can be chosen in the 
        interactive mode.
    autoclose=True : bool
        Says whether to close the intermediate eigenvalues' plot or not.
    
    Returns
    -------
    results : np.array
        Parameters that best fit the data. On its columns it holds...
        ...frequency :math:`f=2\pi\omega` in Hz.
        ...characteristic time :math:`\tau_i` in ps.
        ...quality factors :math:`Q_i=\frac{\omega}{2\gamma}=\pi f \tau`
        ...amplitudes :math:`A_i` in the same units as :math:`x`
        ...phases :math:`\phi_i` written in multiples of :math:`\pi`
    other_results : dict
        Other fit parameters...
        ...chi squared :math:`\chi^2`
        ...number of significant values :math:`N`
    plot_results : ivu.InstancesDict
        Fit parameters that allow plotting. In particular, it holds...
        ...'fit' which includes time, data, fit and fit terms.
        ...'raman' which includes frequency, fit spectrum and fit terms spectra.
    
    See also
    --------
    ivp.linearPredictionPlot
    
    """
    
    #%% ---------------------------------------------------------------------------
    # FREQUENCIES AND DAMPING FACTORS
    # -----------------------------------------------------------------------------
    
    # Create data matrix
    N = len(x)
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
    if svalues is None:
        fig = plt.figure()
        ax = plt.subplot()
        plt.semilogy(eigenvalues, linestyle='none', marker='o', 
                     fillstyle='none', markersize=10)
        plt.title('¿Número de valores singulares?')
        plt.ylabel("Autovalores")
        Nsignificant = ivp.interactiveIntegerSelector(ax, 
                                                      min_value=0, 
                                                      max_value=max_svalues)
        if autoclose:
            plt.close(fig)
    else:
        Nsignificant = svalues
    
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
    frequencies = 1000 * angular_frequencies / (2*pi) # in GHz
    amplitudes = np.array(amplitudes)
    phases = np.array(phases)
    pi_phases = phases / pi # in radians written as multiples of pi
    if Nfit_terms==0:
        raise ValueError("¡Error! No se encontraron términos de ajuste")
    elif Nfit_terms>1:
        print("¡Listo! Encontramos {} términos".format(Nfit_terms))
    else:
        print("¡Listo! Encontramos {} término".format(Nfit_terms))
    print("Frecuencias: {} GHz".format(frequencies))
    
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
    print("Chi cuadrado \u03C7\u00B2: {:.2e}".format(chi_squared))
                     
    ## Statistics of the residue
    #residue = x - fit
    #residue_transform = abs(np.fft.rfft(residue))
    #residue_frequencies = 1000 * np.fft.rfftfreq(N, d=dt) # in GHz
    #plt.plot(residue_frequencies, residue_transform)
    
    # Raman-like Spectrum parameters
    max_frequency = max(frequencies)
    frequencies_damping = 1000 * damping_constants / (2*pi) # in GHz
    if max_frequency != 0:
        raman_frequencies = np.arange(0, 1.5*max_frequency, max_frequency/1000)
    else:
        raman_frequencies = np.array([0, 12])
        
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

    # Some other results I need to plot
    other_results = dict(chi_squared = chi_squared,
                         Nsingular_values = Nsignificant)
    # And the data to plot
    plot_results = ivu.InstancesDict(dict(
            fit = np.array([t, x, fit, *list(fit_terms.T)]).T, 
            raman = np.array([raman_frequencies, raman_spectrum,
                              *list(raman_spectrum_terms.T)]).T))
    
    return results, other_results, plot_results
   
#%%

def linearPredictionTables(parameters, results, other_results):

    terms_heading = ["F (GHz)", "\u03C4 (ps)", "Q", "A (u.a.)", "Fase (\u03C0rad)"]
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

def arrayTable(array, heading_list=None, axis=0):
    
    if heading_list is not None:
        heading = '\t'.join(heading_list)
    if axis==1:
        array = array.T
    elif axis!=0:
        raise ValueError("Axis must be 0 or 1!")
    items = ['\t'.join([str(element) for element in row]) for row in array]
    items = '\n'.join(items)
    if heading_list is not None:
        table = '\n'.join([heading, items])
    else:
        table = items
    
    return table