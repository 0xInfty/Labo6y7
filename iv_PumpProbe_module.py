# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 18:00:29 2019

@author: Vall
"""

import iv_save_module as ivs
import iv_plot_module as ivp
import iv_analysis_module as iva
import iv_utilities_module as ivu
import matplotlib.pyplot as plt
#import matplotlib.widgets as wid
import numpy as np
import os

#%%

def loadPumpProbe(filename):
    
    """Retrieves data and details from a PumpProbe measurement's file.
    
    Each PumpProbe file starts with some data heading like this:
            
        '''
        Formula 
        Fecha   10/04/2019  13:49 
        Desde  -40.00 
        Hasta  1320.00 
        Paso  2.00 
        Tiempo de Integracion  100.00 
        Retardo cero  -640.00
        '''
        
    These files contain time in ps and voltage on V.
        
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    
    Returns
    -------
    data : np.array
        Measured data. It has 2*N columns, where N is the number of 
        experiments. Inside, it holds data ordered as [T1, V1, ..., TN, VN] 
        where Ti is time in ps and Vi is voltage in V.
    details : dict
        Details of the measurement, including...
            date : str
                Date and hour, on 'DD/MM/YYYY HH:HH' format.
            time_range : tuple
                time_start : float
                    First time difference, in ps.
                time_end : float
                    Last time difference, in ps.
            time_step : float
                Minimum time step, in ps.
            time_integration : float
                Lock-in's integration time, in ps, that defines how much time 
                will the system retain the same time difference in order to 
                make an average reading using the lock-in.
            time_zero : float
                Time reference, in ps.
    
    Raises
    ------
    ValueError : "Columns have different number of rows :("
        When a numpy array cannot be made because there's a faulty experiment, 
        which doesn't hold as much data as it should.
    
    """
   
    lines = []
    other_lines = []
    
    extras = ['Fecha   ', 'Desde  ',  'Hasta  ', 'Paso  ', 
              'Tiempo de Integracion  ', 'Retardo cero  ']
    names = ['date', 'time_range', 'time_step', 
             'time_integration', 'time_zero']
    
    i = 0
    with open(filename, 'r') as f:
        for line in f:
            if i >= 1 and i < 7: # I must ignore the first line
                lines.append(line)
            elif i >= 7: # From here on there's only data.
                other_lines.append(line)
            i = i+1
    
    details = {}
    details[names[0]] = lines[0].split(extras[0])[-1].split(' \n')[0]
    details[names[1]] = (
            float(lines[1].split(extras[1])[-1].split(' \n')[0]),
            float(lines[2].split(extras[2])[-1].split(' \n')[0]),
                             )
    details[names[2]] = float(lines[3].split(extras[3])[-1].split(' \n')[0])
    details[names[3]] = float(lines[4].split(extras[4])[-1].split(' \n')[0])
    details[names[4]] = float(lines[5].split(extras[5])[-1].split(' \n')[0])

#    other_lines = [[float(number) for number in line.split('\t')] 
#                    for line in other_lines]
#    N = len(other_lines) # Number of measurements each experiment should have.
#
#    data = []
#    for i in range(N):
#        for experiment in range(len(other_lines[0])/2):
#            if other_lines[i][]

    try:
        data = np.array([[float(number) for number in line.split('\t')] 
                        for line in other_lines])   
    except:
        raise ValueError("Columns have different number of rows :(")
        
    return data, details
    
#%%

def loadNicePumpProbe(filename):

    """Loads nice data and details from a PumpProbe measurement's file.
    
    Returns equispaced time in ps, voltage in uV and also calculates mean voltage 
    in uV. Moreover, it adds some parameters to the measurement's details.
    
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    
    Returns
    -------
    t : np.array
        Equispaced time in ps. It has 'nsize' length.
    V : np.array
        Measured voltage in uV. It has 'nsize' rows and 'nrepetitions' columns.
    details : dict
        Details of the measurement, including...
            samplerate : float
                Sampling rate in Hz.
            dt : float
                Time step in ps of the equispaced time.
            nsize : int
                Number of measurements included in each repetition.
            nexperiments : int
                Number of repetitions of the experiment.
    
    Raises
    ------
    ValueError : "Columns have different number of rows :("
        When a numpy array cannot be made because there's a faulty experiment, 
        which doesn't hold as much data as it should.
    
    See also
    --------
    loadPumpProbe
    """
    
    # First get data's name
    [data, details] = loadPumpProbe(filename)
    
    # Define data size
    nrepetitions = int(len(data[0,:]) / 2) # Number of measurements
    nsize = len(data[:,0]) # Length of each measurement
    
    # Get time
    t = data[:,0] # Consider just one time column
    
    # Define time parameters
    T = t[-1] - t[0] # Total time in ps
    samplerate = nsize / (10**12 * T)  # Sampling rate in Hz
    dt = T / nsize # Time step
    
    # Make equispaced time
    t = np.linspace(t[0], t[-1], nsize)
    
    # Add uV voltage
    V = np.array([1e6 * data[:, 2*i+1] for i in range(nrepetitions)]).T
    
    # Add some other relevant details
    details.update(dict(samplerate=samplerate,
                        dt=dt,
                        nsize=nsize,
                        nrepetitions=nrepetitions))
    
    return t, V, details

#%%

def plotPumpProbe(filename, interactive=False, autosave=True, **kwargs):

    """Plots all PumpProbe experiments from a file and its mean.
    
    Can also make an interactive plot, which holds a save button and allows to 
    choose only certain experiments to be shown from the legend.
    
    By default, it also saves a picture on the file's path.
        
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    interactive=True : bool
        Says whether to make an interactive plot or not.
    autosave=True : bool
        Says whether to automatically save or not.
    
    Returns
    -------
    fig : plt.Figure instance
        Figure containing the desired plot.
    legend_buttons : wid.CheckButtons
        Interactive legend. Only returned if making an interactive plot.
    save_button : wid.Button
        Interactive save button. Only returned if making an interactive plot.        
    
    Raises
    ------
    pngfile : .png file
        PNG image file. Only raised if 'autosave=True'.
    
    See also
    --------
    ivs.loadPumpProbe
    ivp.interactiveLegend
    ivp.interactiveSaveButton
    
    """
    
    path = os.path.join(os.path.split(filename)[0], 'Figuras')
    name = os.path.split(os.path.splitext(filename)[0])[-1]
    t, V, details = loadNicePumpProbe(filename)
    meanV = np.mean(V, axis=1)
    Nrepetitions = details['nrepetitions']
    
    fig = plt.figure(figsize=[6.4, 4.4])
    ax = plt.subplot()
    plt.plot(t, V, linewidth=0.8, zorder=0)
    plt.plot(t, meanV, linewidth=1.5, zorder=2)
    labels = ['Experimento {:.0f}'.format(i+1) for i in range(Nrepetitions)]
    labels.append('Promedio')
    plt.ylabel(r'Voltaje ($\mu$V)', fontsize=14)
    plt.xlabel(r'Tiempo (ps)', fontsize=14)
    ax.tick_params(labelsize=12)
    ax.minorticks_on()
    ax.tick_params(axis='y', which='minor', left=False)
    ax.tick_params(length=5)
    ax.grid(axis='x', which='both')
    
    ax = fig.axes[0]
    position = ax.get_position()
    ax.set_position([position.x0*1.2, position.y0*1.3,
                     position.width, position.height])
    
    if interactive:
        show_default = [True for lab in labels]
        legend_buttons = ivp.interactiveLegend(ax, labels, show_default, 
                                           fontsize=12,
                                           x0=(.17, .68), y0=(.06, .84), **kwargs)
        save_button = ivp.interactiveSaveButton(filename)
    else:
        plt.legend(labels, fontsize=12, framealpha=1, **kwargs)
    
    if not os.path.isdir(path):
        os.makedirs(path)
    
    if autosave:
        if interactive:
            save_button.ax.set_visible(False)
        plt.savefig(os.path.join(path,name+'_fig.png'), bbox_inches='tight')
        if interactive:
            save_button.ax.set_visible(True)
    
    if interactive:
        return fig, legend_buttons, save_button
    else:
        return fig

#%%

def plotAllPumpProbe(path, autosave=True, autoclose=False):
    
    """Plots all PumpProbe experiments on the files from a given path.
        
    The data files must be '.txt' files that begin with 'M'.
    
    Parameters
    ----------
    path : str
        Files' folder (must include directory).
    autosave=True : bool
        Says whether to save or not.
    autoclose=False : bool
        Says whether to close the figures or not.
    
    Returns
    -------
    figures : list
        A list containing plt.Figure instances -only returned if autoclose is 
        deactivated.
    
    Raises
    ------
    pngfiles : .png files
        PNG image files. Only raised if 'autosave=True'.
    
    See also
    --------
    ivp.plotPumpProbe
    
    """
    
    files = []
    for file in os.listdir(path):
        if file.endswith(".txt") and file.startswith("M"):
            files.append(os.path.join(path,file))
    
    figures = []
    for f in files:
        fig = plotPumpProbe(f, interactive=False, autosave=autosave)
        if autoclose:
            plt.close(fig)
        else:
            figures.append(fig)
    
    if not autoclose:
        return figures

#%%
        
class FitConfigPumpProbe(ivu.InstancesDict):
    
    def __init__(self, **kwargs):
        
        default_properties = dict(
                use_full_mean = True,
                use_experiments = [0], # First is 0, not 1!
                send_tail_to_zero = False,
                send_tail_method = 'mean', # Could also be 'min' or 'max' or any numpy function
                send_tail_fraction = .2,
                choose_t0 = True,
                choose_tf = False,
                choose_V0 = False,
                )
        
        dic = {}
        for k, v in default_properties.items():
            try:
                dic[k] = kwargs[k]
            except:
                dic.update({k: v})
        
        super().__init__(dic)

#%%
                
class FitParamsPumpProbe(ivu.InstancesDict):
    
    def __init__(self, **kwargs):
        
        default_properties = dict(
                t0 = None,
                tf = None,
                V0 = None
                )
        
        dic = {}
        for k, v in default_properties.items():
            try:
                dic[k] = kwargs[k]
            except:
                dic.update({k: v})
        
        super().__init__(dic)

#%%

def linearPredictionPumpProbe(filename, autosave=True, autoclose=True, **kwargs):
    
    # Load data
    t, V, details = ivs.loadNicePumpProbe(filename) # t in ps, V in uV
    
    # Set configuration
    fit_config = FitConfigPumpProbe(**kwargs)
    fit_params = FitParamsPumpProbe(**kwargs)
    
    # Choose specific data to fit
    if fit_params.use_full_mean:
        V = np.mean(V, axis=1)
    else:
        V = np.mean(V[:, fit_params.use_experiments], axis=1)

    # Choose initial time t0
    if fit_config.choose_t0 and fit_params.t0 is None:
        fit_params.t0 = ivp.interactiveTimeSelector(filename, 
                                                    autoclose=autoclose)
    elif fit_params.t0 is None:
        fit_params.t0 = t[0]
    
    # Choose final time tf
    if fit_config.choose_tf and fit_params.tf is None:
        fit_params.tf = ivp.interactiveTimeSelector(filename, 
                                                    autoclose=autoclose)
    elif fit_params.tf is None:
        fit_params.tf = t[-1]
    
    # Crop data for the chosen time interval
    t, V = iva.cropData(fit_params.t0, t, V)
    t, V = iva.cropData(fit_params.tf, t, V, logic='<=')
    
    # If needed, make a vertical shift
    if fit_config.send_tail_to_zero and fit_config.choose_V0:
        print("¡Ojo! Corrimiento vertical automático")
        fit_config.choose_V0 = False
    if fit_config.send_tail_to_zero and fit_params.V0 is None:
        function = eval('np.{}'.format(fit_config.send_tail_method))
        fit_params.V0 = function(V[int( (1-fit_config.send_tail_fraction) * len(V)):])
        del function
    elif fit_params.V0 is None:
        fit_params.V0 = 0
    V = V - fit_params.V0
    
    # Use linear prediction
    results, parameters = iva.linearPrediction(t, V, dt=details['dt'], 
                                               autoclose=autoclose)
    results[:,0] = results[:,0] / 1000 # t is ps, so f was THz but now is GHz
    
#    if autosave:
#        ivs.linearPredictionSave(filename, results, other_results, fit_params)
#    
#    # Plot linear prediction
#    ivp.linearPredictionPlot(filename, plot_results, autosave=autosave)
#    
#    # Generate fit tables
#    tables = iva.linearPredictionTables(fit_params, results, other_results)
#    ivu.copy(tables[0]) 

#%%
 
def interactiveTimeSelector(filename, autoclose=True):
    
    """Allows to select a particular time instant on a Pump Probe file.
    
    Parameters
    ----------
    filename : str
        Filename, which must include full path and extension.
    autoclose=True : bool
        Says whether to automatically close this picture or not.
    
    Returns
    -------
    ti : float
        Selected value.
    
    See also
    --------
    ivs.loadNicePumpProbe
    
    """
    
    t, V, details = loadNicePumpProbe(filename)
    fig = plotPumpProbe(filename, autosave=False)
    ax = fig.axes[0]
    ti = ivp.interactiveValueSelector(ax, y_value=False)
    ti = t[np.argmin(abs(t-ti))]
    
    if autoclose:
        plt.close(fig)

    return ti
   
#%%
"""
def linearPredictionPumpProbeTables(results, other_results, units=["Hz","s"]):

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
"""