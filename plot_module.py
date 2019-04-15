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
    
    path = os.path.join(os.path.split(file)[0], 'Figuras')
    name = os.path.split(os.path.splitext(file)[0])[-1]
    data = loadPumpProbe(file)[0]
    
    N = len(data[0,:]) // 2 # Number of experiments
    
    t = np.array([data[:,2*i] for i in range(N)]).T
    V = np.array([data[:,2*i+1] for i in range(N)]).T * 1e6
    meanV = np.mean(V, axis=1)
    meant = t[:,0] #np.mean(t, axis=1)
    
    plt.figure()
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
    
def fullplotPumpProbe(file, save=True):
    
    path = os.path.join(os.path.split(file)[0], 'Figuras')
    name = os.path.split(os.path.splitext(file)[0])[-1]
    data = loadPumpProbe(file)[0]
    
    N = len(data[0,:]) // 2 # Number of experiments
    
    t = np.array([data[:,2*i] for i in range(N)]).T
    V = np.array([data[:,2*i+1] for i in range(N)]).T * 1e6
    meanV = np.mean(V, axis=1)
    meant = t[:,0] #np.mean(t, axis=1)
    
    plt.figure()
    grid = plt.GridSpec(N, 2, wspace=0.4, hspace=0.3)
    
    for i in range(N):
        plt.subplot(N, 2, 2*i+1)
        plt.plot(t[:,i], V[:,i], linewidth=0.8)
    
    plt.subplot(grid[0,1])
    plt.plot(meant, meanV, 'r', linewidth=0.8)
    
    plt.subplot(grid[1:,1])
    plt.plot(t, V, linewidth=0.8)
    plt.plot(meant, meanV, linewidth=1.5)
   
def plotallPumpProbe(path, save=True):
    
    files = []
    for file in os.listdir(path):
        if file.endswith(".txt"):
            files.append(os.path.join(path,file))
    
    for f in files:
        plotPumpProbe(f, save=save)
    
    if save:
        plt.savefig(os.path.join(path,name+'_full.png'), bbox_inches='tight')
