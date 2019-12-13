# -*- coding: utf-8 -*-
"""The 'save' module loads and saves data, dealing with overwriting.

It could be divided into 3 sections:
    (1) making new directories and free files to avoid overwriting 
    ('newDir', 'freeFile')
    (2) loading data from PumpProbe experiments ('loadPumpProbe', 
    'loadNicePumpProbe', 'filenameTo...')
    (2) saving data into files with the option of not overwriting 
    ('saveFig', 'saveTxt')
    (4) loading data from files ('retrieveHeader', 'retrieveFooter')

@author: Vall
"""

import iv_utilities_module as ivu
import numpy as np
import os
import matplotlib.pyplot as plt

#%%

def newDir(my_dir, newformat='{}_{}'):
    
    """Makes and returns a new directory to avoid overwriting.
    
    Takes a directory name 'my_dir' and checks whether it already 
    exists. If it doesn't, it returns 'dirname'. If it does, it 
    returns a related unoccupied directory name. In both cases, 
    the returned directory is initialized.
    
    Parameters
    ----------
    my_dir : str
        Desired directory (should also contain full path).
    
    Returns
    -------
    new_dir : str
        New directory (contains full path)
    
    Yields
    ------
    new_dir : directory
    
    """
    
    sepformat = newformat.split('{}')
    base = os.path.split(my_dir)[0]
    
    new_dir = my_dir
    while os.path.isdir(new_dir):
        new_dir = os.path.basename(new_dir)
        new_dir = new_dir.split(sepformat[-2])[-1]
        try:
            new_dir = new_dir.split(sepformat[-1])[0]
        except ValueError:
            new_dir = new_dir
        try:
            new_dir = newformat.format(my_dir, str(int(new_dir)+1))
        except ValueError:
            new_dir = newformat.format(my_dir, 2)
        new_dir = os.path.join(base, new_dir)
    os.makedirs(new_dir)
        
    return new_dir

#%%

def filenameCreator(my_file, folder='', sufix=''):
    
    """Returns same filename in a different folder and with a sufix if given.
    
    Parameters
    ----------
    my_file : str
        Tentative filename (must contain full path and extension).
    folder='' : str
        A new folder or set of folders to be added to the folder.
    sufix='' : str
        A sufix to be added to the given filename.
    
    Returns
    -------
    new_file : str
        New filename (also contains full path and same extension)
    """

    base = os.path.split(my_file)[0]
    if folder!='':
        base = os.path.join(base, folder)
    extension = os.path.splitext(my_file)[-1]
    new_file = os.path.join(
        base, 
        os.path.splitext( os.path.split(my_file)[1] )[0] + sufix + extension)
    
    if not os.path.isdir(base):
        os.makedirs(base)
    
    return new_file

#%%

def freeFile(my_file, newformat='{}_v{}'):
    
    """Returns a name for a new file to avoid overwriting.
        
    Takes a file name 'my_file'. It returns a related unnocupied 
    file name 'free_file'. If necessary, it makes a new 
    directory to agree with 'my_file' path.
        
    Parameters
    ----------
    my_file : str
        Tentative filename (must contain full path and extension).
    newformat='{}_v{}' : str
        A formatter that allows to make new filenames in order to avoid 
        overwriting. If 'F:\Hola.png' does already exist, new file is saved as 
        'F:\Hola_v2.png'.
    
    Returns
    -------
    free_file : str
        Unoccupied filename (also contains full path and extension).
    
    """
    
    base = os.path.split(my_file)[0]
    extension = os.path.splitext(my_file)[-1]
    
    if not os.path.isdir(base):
        os.makedirs(base)
        free_file = my_file
    
    else:
        sepformat = newformat.split('{}')[-2]
        free_file = my_file
        while os.path.isfile(free_file):
            free_file = os.path.splitext(free_file)[0]
            free_file = free_file.split(sepformat)
            number = free_file[-1]
            free_file = free_file[0]
            try:
                free_file = newformat.format(
                        free_file,
                        str(int(number)+1),
                        )
            except ValueError:
                free_file = newformat.format(
                        os.path.splitext(my_file)[0], 
                        2)
            free_file = os.path.join(base, free_file+extension)
    
    return free_file

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

def filenameToDate(filename):
    
    """Given a filename 'M_20190610_01', returns date on '2019-06-10' format"""
    
    date = filename.split('_')[1] # From 'M_20190610_01' take '20190610'
    date = '-'.join([date[:4], date[4:6], date[6:]]) # Transfrom to '2019-06-10'
    
    return date
 
#%%

def filenameToMeasureFilename(filename, 
                              home=os.getcwd()):
    
    """Given a filename 'M_20190610_01', returns path to fits' data"""
    
    date = filenameToDate(filename) # Transfrom to '2019-06-10'
    fits_filename = os.path.join(home, 'Mediciones', date, filename+'.txt')
    
    return fits_filename

#%%

def filenameToFitsFilename(filename,
                           home=os.getcwd()):
    
    """Given a filename 'M_20190610_01', returns path to fits' data"""
    
    date = filenameToDate(filename) # Transfrom to '2019-06-10'
    fits_filename = os.path.join(home, 'Mediciones', date, 
                                 'Ajustes', filename+'.txt')
    
    return fits_filename

#%%

def linearPredictionSave(filename, results, other_results, fit_parameters, 
                         overwrite=False):
    
    """Saves the data from a linear prediction fit on '.txt' file.
    
    Parameters
    ----------
    filename : str
        The name you wish (must include full path and extension)
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
    fit_parameters : ivu.InstancesDict
        Several fit configuration parameters, including...
            use_full_mean=True : bool
                Whether to use full mean or not.
            send_tail_to_zero=False : bool
                Whether to apply a vertical shift to send the last data to zero 
                or not.
            voltage_zero : float, int
                Vertical shift.
            time_range : tuple
                Initial and final time to fit.
    overwrite=False
        Whether to allow overwriting or not.
    
    Returns
    -------
    None
    
    Yields
    ------
    .txt file
    
    See also
    --------
    saveTxt
    
    """
        
    fit_params = fit_parameters.__dict__ # Because it's an ivu.InstancesDict
    
    footer = {}
    footer.update(other_results)
    footer.update(fit_params)
    
    saveTxt(filename, results,
            header=["F (GHz)", "Tau (ps)", "Q", "A (u.a.)", "Phi (pi rad)"],
            footer=footer,
            overwrite=overwrite,
            folder='Ajustes')
    
    return

#%%
    
def saveFig(filename, extension='.png', overwrite=False, 
            newformat='{}_v{}', **kwargs):
    
    """Saves current matplotlib figure in a compact format.
    
    This function takes per default the current plot and saves it on file. If 
    'overwrite=False', it checks whether 'filename' exists or not; if it already 
    exists, it defines a new file in order to not allow overwritting. If 
    overwrite=True, it saves on 'filename' even if it already exists.
    
    Variables
    ---------
    filename : str
        The name you wish (must include full path).
    extension='.png' : str
        An image extension; i.e.: '.pdf', '.jpg'.
    overwrite=False : bool, optional
        Indicates whether to overwrite or not.
    newformat='{}_v{}' : str
        A formatter that allows to make new filenames in order to avoid 
        overwriting. If 'F:\Hola.png' does already exist, new file is saved as 
        'F:\Hola_v2.png'.
    
    Other parameters
    ----------------
    folder='' : str
        A new folder or set of folders to be added to the folder.
    sufix='' : str
        A sufix to be added to the given filename.
    
    Return
    ------
    nothing
    
    Yield
    -----
    image file
    
    See Also
    --------
    freeFile
    
    """
    
    filename = os.path.splitext(filename)[0] + extension
    filename = filenameCreator(filename, **kwargs)
    
    if not overwrite:
        filename = freeFile(filename, newformat=newformat)
    
    plt.savefig(filename, bbox_inches='tight')
    
    print('Imagen guardada en {}'.format(filename))
    
    return

#%%

def saveTxt(filename, datanumpylike, header='', footer='', 
            overwrite=False, newformat='{}_v{}', **kwargs):
    
    """Takes some array-like data and saves it on a '.txt' file.
    
    This function takes some data and saves it on a '.txt' file.
    If 'overwrite=False', it checks whether 'filename' exists or not; if it 
    already exists, it defines a new file in order to not allow 
    overwritting. If overwrite=True, it saves on 'filename' even if 
    it already exists.
    
    Variables
    ---------
    filename : string
        The name you wish (must include full path and extension)
    datanumpylike : array, list
        The data to be saved.
    overwrite=False : bool, optional
        Indicates whether to overwrite or not.
    header='' : list, str, optional
        Data's descriptor. Its elements should be str, one per column.
        But header could also be a single string.
    footer='' : dict, str, optional
        Data's specifications. Its elements and keys should be str. 
        But footer could also be a single string. Otherwise, an element 
        could be a tuple containing value and units; i.e.: (100, 'Hz').
    newformat='{}_v{}' : str
        A formatter that allows to make new filenames in order to avoid 
        overwriting. If 'F:\Hola.png' does already exist, new file is saved as 
        'F:\Hola_v2.png'.
    
    Other parameters
    ----------------
    folder='' : str
        A new folder or set of folders to be added to the folder.
    sufix='' : str
        A sufix to be added to the given filename.
    
    Return
    ------
    nothing
    
    Yield
    -----
    '.txt' file
    
    See Also
    --------
    freeFile
    
    """
        
    if header != '':
        if not isinstance(header, str):
            try:
                header = '\t'.join(header)
            except:
                TypeError('Header should be a list or a string')

    if footer != '':
        if not isinstance(footer, str):
            try:
                aux = []
                for key, value in footer.items():
                    if isinstance(value, tuple) and len(value) == 2:
                        condition = isinstance(value[0], str)
                        if not condition and isinstance(value[1], str):
                            value = '"{} {}"'.format(*value)
                    elif isinstance(value, str):
                        value = '"{}"'.format(value)
                    aux.append('{}={}'.format(key, value) + ', ')
                footer = ''.join(aux)
            except:
                TypeError('Header should be a dict or a string')
   
    filename = filenameCreator(filename, **kwargs)
    if not overwrite:
        filename = freeFile(filename, newformat=newformat)
        
    np.savetxt(filename, np.array(datanumpylike), 
               delimiter='\t', newline='\n', header=header, footer=footer)
    
    print('Archivo guardado en {}'.format(filename))
    
    return

#%%

def retrieveFooter(filename, comment_marker='#'):
    
    """Retrieves the footer of a .txt file saved with np.savetxt or saveTxt.
    
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    comment_marker='#' : str, optional
        Sign that indicates a line is a comment on np.savetxt.
    
    Returns
    -------
    last_line : str, dict
        File's footer
    
    Raises
    ------
    ValueError : "Footer not found. Sorry!"
        When the last line doesn't begin with 'comment_marker'.
        
    See Also
    --------
    saveTxt
    
    """
    
    
    with open(filename, 'r') as f:
        for line in f:
            last_line = line
    
    if last_line[0] == comment_marker:
        try:
            last_line = last_line.split(comment_marker + ' ')[-1]
            last_line = last_line.split('\n')[0]
            footer = eval('dict({})'.format(last_line))
            for key, value in footer.items():
                try:
                    number = ivu.findNumbers(value)
                    if len(number) == 1:
                        number = number[0]
                        if len(value.split(' ')) == 2:
                            footer[key] = (
                                number, 
                                value.split(' ')[-1]
                                )
                        else:
                            footer[key] = number
                except TypeError:
                    value = value
        except:
            footer = last_line
        return footer
        
    else:
        raise ValueError("No footer found. Sorry!")

#%%

def retrieveHeader(filename, comment_marker='#'):
    
    """Retrieves the header of a .txt file saved with np.savetxt or saveTxt.
    
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    comment_marker='#' : str, optional
        Sign that indicates a line is a comment on np.savetxt.
    
    Returns
    -------
    last_line : str, list
        File's header
    
    Raises
    ------
    ValueError : "Header not found. Sorry!"
        When the first line doesn't begin with 'comment_marker'.
    
    See Also
    --------
    saveTxt
    
    """
    
    
    with open(filename, 'r') as f:
        for line in f:
            first_line = line
            break
    
    if first_line[0] == comment_marker:
        header = first_line.split(comment_marker + ' ')[-1]
        header = header.split('\n')[0]
        header = header.split('\t')
        if len(header) > 1:
            return header
        else:
            return header[0]
        
    else:
        raise ValueError("No header found. Sorry!")

#%%

def loadTxt(filename, comment_marker='#', **kwargs):
    
    """Loads data of a .txt file saved with np.savetxt or saveTxt.

    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    comment_marker='#' : str, optional
        Sign that indicates a line is a comment on np.savetxt.

    Returns
    -------
    data : np.array
        File's data.
    header : str, list or None
        File's header.
    footer : str, dict or None
        File's footer.

    See also
    --------
    saveTxt
    retrieveHeader
    retrieveFooter
    
    """
    
    data = np.loadtxt(filename, comments=comment_marker, **kwargs)
    try:
        header = retrieveHeader(filename, comment_marker)
    except ValueError:
        header = None
    try:
        footer = retrieveFooter(filename, comment_marker)
    except ValueError:
        footer = None
    
    return data, header, footer