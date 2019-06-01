# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 18:00:29 2019

@author: Usuario
"""

import iv_plot_module as ivp
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