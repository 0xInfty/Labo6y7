# -*- coding: utf-8 -*-
"""The 'save' module loads and saves data, dealing with overwriting.

It could be divided into 3 sections:
    (1) making new directories and free files to avoid overwriting 
    ('newDir', 'freeFile')
    (2) loading data from PumpProbe experiments ('loadPumpProbe', 
    'loadNicePumpProbe')
    (2) saving data into files with the option of not overwriting 
    ('saveTxt')
    (4) loading data from files ('retrieveHeader', 'retrieveFooter')

new_dir : function
    Makes and returns a new related directory to avoid overwriting.
free_file : function
    Returns a name for a new file to avoid overwriting.

savetxt : function
    Saves some np.array like data on a '.txt' file.

@author: Vall
"""

import iv_utilities_module as ivu
import numpy as np
import os

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

def freeFile(my_file, newformat='{}_{}'):
    
    """Returns a name for a new file to avoid overwriting.
        
    Takes a file name 'my_file'. It returns a related unnocupied 
    file name 'free_file'. If necessary, it makes a new 
    directory to agree with 'my_file' path.
        
    Parameters
    ----------
    my_file : str
        Tentative filename (must contain full path and extension).
    newformat='{}_{}' : str
        Format string that indicates how to make new names.
    
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
    
    path, name = os.path.split(filename)
    filename = os.path.join(path, 'Ajustes', name)
    
    fit_params = fit_parameters.__dict__ # Because it's an ivu.InstancesDict
    
    footer = {}
    footer.update(other_results)
    footer.update(fit_params)
    
    saveTxt(filename, results,
            header=["F (GHz)", "Tau (ps)", "Q", "A (u.a.)", "Phi (pi rad)"],
            footer=footer,
            overwrite=overwrite, newformat='{}_v{}')
    
    return

#%%

def saveTxt(filename, datanumpylike, header='', footer='', 
            overwrite=False, newformat='{}_{}'):
    
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
    
    base = os.path.split(filename)[0]
    if not os.path.isdir(base):
        os.makedirs(base)
    
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

    filename = os.path.join(
            base,
            (os.path.splitext(os.path.basename(filename))[0] + '.txt'),
            )
    
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