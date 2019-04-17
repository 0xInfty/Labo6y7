# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:08:08 2019

@author: LEC
"""

from loadPumpProbe import loadPumpProbe
import matplotlib.pyplot as plt
import numpy as np
import os

def plotPumpProbe(file, save=True):

    """Plots all PumpProbe experiments from a file and its mean.
        
    Parameters
    ----------
    file : str
        File's root (must include directory and termination).
    save=True : bool
        Says whether to save or not.
    
    Returns
    -------
    fig : matplotlib.pyplot.Figure instance
        Figure containing the desired plot.
    
    Raises
    ------
    pngfile : .png file
        PNG image file (only if 'save=True').
    
    See also
    --------
    loadPumpProbe
    
    """
    
    path = os.path.join(os.path.split(file)[0], 'Figuras')
    name = os.path.split(os.path.splitext(file)[0])[-1]
    data = loadPumpProbe(file)[0]
    
    N = len(data[0,:]) // 2 # Number of experiments
    
    t = np.array([data[:,2*i] for i in range(N)]).T
    V = np.array([data[:,2*i+1] for i in range(N)]).T * 1e6
    meanV = np.mean(V, axis=1)
    meant = t[:,0]
    
    fig = plt.figure()
    plt.plot(t, V, linewidth=0.8)
    plt.plot(meant, meanV, linewidth=1.5)
    legends = ['Experimento {:.0f}'.format(i+1) for i in range(N)]
    legends.append('Promedio')
    plt.legend(legends)
    plt.ylabel(r'Voltaje ($\mu$V)')
    plt.xlabel(r'Tiempo (ps)')
    
    if not os.path.isdir(path):
        os.makedirs(path)
    
    if save:
        plt.savefig(os.path.join(path,name+'_fig.png'), bbox_inches='tight')
    
    return fig
    
def fullplotPumpProbe(file, save=True):
    
    """Plots all PumpProbe experiments from a file on a set of subplots.
        
    Parameters
    ----------
    file : str
        File's root (must include directory and termination).
    save=True : bool
        Says whether to save or not.
    
    Returns
    -------
    fig : matplotlib.pyplot.Figure instance
        Figure containing the desired plot.
    
    Raises
    ------
    pngfile : .png file
        PNG image file (only if 'save=True').
    
    See also
    --------
    loadPumpProbe
    
    """
    
    path = os.path.join(os.path.split(file)[0], 'Figuras')
    name = os.path.split(os.path.splitext(file)[0])[-1]
    data = loadPumpProbe(file)[0]
    
    N = len(data[0,:]) // 2 # Number of experiments
    
    t = np.array([data[:,2*i] for i in range(N)]).T
    V = np.array([data[:,2*i+1] for i in range(N)]).T * 1e6
    meanV = np.mean(V, axis=1)
    meant = t[:,0] #np.mean(t, axis=1)
    
    fig = plt.figure()
    grid = plt.GridSpec(N, 2, wspace=0.4, hspace=0.3)
    
    for i in range(N):
        plt.subplot(N, 2, 2*i+1)
        plt.plot(t[:,i], V[:,i], linewidth=0.8)
    
    plt.subplot(grid[0,1])
    plt.plot(meant, meanV, 'r', linewidth=0.8)
    
    plt.subplot(grid[1:,1])
    plt.plot(t, V, linewidth=0.8)
    plt.plot(meant, meanV, linewidth=1.5)

    if not os.path.isdir(path):
        os.makedirs(path)
    
    if save:
        plt.savefig(os.path.join(path,name+'_full.png'), bbox_inches='tight')
    
    return fig
   
def plotallPumpProbe(path, full=False, save=True):
    
    """Plots all PumpProbe experiments on the files from a given path.
        
    Parameters
    ----------
    file : str
        File's root (must include directory and termination).
    save=True : bool
        Says whether to save or not.
    
    Returns
    -------
    matplotlib.pyplot figure
    png image files
    
    See also
    --------
    plotPumpProbe
    fullplotPumpProbe
    
    """
    
    files = []
    for file in os.listdir(path):
        if file.endswith(".txt"):
            files.append(os.path.join(path,file))
    
    for f in files:
        if full:
            fullplotPumpProbe(f, save=save);
        else:
            plotPumpProbe(f, save=save);
