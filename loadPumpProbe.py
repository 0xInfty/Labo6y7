# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 10:10:42 2019

@author: quimica
"""

"""
En los archivos que guarda el PumpProbe
El tiempo está en ps
Pero el voltaje está en V
`\(T.T)/´
"""

import numpy as np

def loadPumpProbe(file):
    
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
        
    Parameters
    ----------
    file : str
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
    with open(file, 'r') as f:
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
