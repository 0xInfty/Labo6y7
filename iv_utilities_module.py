# -*- coding: utf-8 -*-
"""The 'utilities' module contains a few generic tools.

@author: Vall
"""

from tkinter import Tk
import re

#%%

def copy(string):
    
    """Copies a string to the clipboard.
    
    Parameters
    ----------
    string : str
        The string to be copied.
    
    Returns
    -------
    nothing
    
    """
    
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(string)
    r.update() # now it stays on the clipboard
    r.destroy()
    
    print("Copied")

#%%

def findNumbers(string):
    
    """Returns a list of numbers found on a given string
    
    Parameters
    ----------
    string: str
        The string where you search.
    
    Returns
    -------
    list
        A list of numbers (each an int or float).
    
    Raises
    ------
    "There's no number in this string" : TypeError
        If no number is found.
    """
    
    numbers = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", string)
    
    if not numbers:
        raise TypeError("There's no number in this string")
    
    for i, n in enumerate(numbers):
        if '.' in n:
            numbers[i] = float(n)
        else:
            numbers[i] = int(n) 
    
    return numbers

#%%

def countingSufix(number):
    
    """Returns a number's suffix string to use for counting.
    
    Parameters
    ----------
    number: int, float
        Any number, though it is designed to work with integers.
    
    Returns
    -------
    ans: str
        A string representing the integer number plus a suffix.
    
    Examples
    --------
    >> counting_sufix(1)
    '1st'
    >> counting_sufix(22)
    '22nd'
    >> counting_sufix(1.56)
    '2nd'
    
    """
    
    number = round(number)
    unit = int(str(number)[-1])
    
    if unit == 1:
        ans = 'st'
    if unit == 2:
        ans = 'nd'
    if unit == 3:
        ans = 'rd'
    else:
        ans = 'th'
    
    return ans

#%%

class InstancesDict:
    
    """Example of a class that holds a callable dictionary of instances.
    Examples
    --------
    >> class MyClass:
        def __init__(self, value=10):
            self.sub_prop = value
    >> instance_a, instance_b = MyClass(), MyClass(12)
    >> instance_c = ClassWithInstances(dict(a=instance_a,
                                            b=instance_b))
    >> Z = InstancesDict({1: instance_a,
                          2: instance_b,
                          3: instance_c})
    >> Z(1)
    <__main__.MyClass at 0x2e5572a2dd8>
    >> Z(1).sub_prop
    10
    >> Z(1).sub_prop = 30
    >> Z(1).sub_prop
    >> Z(3).b.sub_prop
    12
    >> Z(1,2)
    [<__main__.MyClass at 0x2e5573cfb00>, 
    <__main__.MyClass at 0x2e5573cf160>]
    
    Warnings
    --------
    'Z(1,2).prop' can't be done.
    
    """
    
    def __init__(self, dic):#, methods):
        
        self.__dict__.update(dic)
    
    def __call__(self, *key):

        if len(key) == 1:
            return self.__dict__[key[0]]
        
        else:
            return [self.__dict__[k] for k in key]
        
    def __repr__(self):
        
        return str(self.__dict__)
    
    def __str__(self):
        
        return str(self.__dict__)
                
    def update(self, dic):
        
        self.__dict__.update(dic)
    
    def is_empty(self, key):
        
        if key in self.__dict__.keys():
            return False
        else:
            return True

#%%

#def errorValue(X, dX, error_digits=2, units='',
#               string_scale=True, one_point_scale=False, 
#               legend=False, unpack=False):
#    
#    """Rounds up value and error of a measure. Also makes a latex string.
#    
#    This function takes a measure and its error as input. Then, it 
#    rounds up both of them in order to share the same amount of decimal 
#    places.
#    
#    After that, it generates a latex string containing the rounded up 
#    measure. For that, it can rewrite both value and error so that the 
#    classical prefix scale of units can be applied.
#    
#    Parameters
#    ----------
#    X : float
#        Measurement's value.
#    dX : float
#        Measurement's associated error.
#    error_digits=2 : int, optional.
#        Desired number of error digits.
#    units='' : str, optional.
#        Measurement's units.
#    string_scale=True : bool, optional.
#        Whether to apply the classical prefix scale or not.        
#    one_point_scale=False : bool, optional.
#        Applies prefix with one order less.
#    legend=False : bool, optional.
#        Says whether it is for the legend of a plot or not.
#    
#    Returns
#    -------
#    latex_str : str
#        Latex string containing value and error.
#    
#    Examples
#    --------
#    >> errorValue(1.325412, 0.2343413)
#    '(1.33$\\pm$0.23)'
#    >> errorValue(1.325412, 0.2343413, error_digits=3)
#    '(1.325$\\pm$0.234)'
#    >> errorValue(.133432, .00332, units='V')
#    '\\mbox{(133.4$\\pm$3.3) mV}'
#    >> errorValue(.133432, .00332, one_point_scale=True, units='V')
#    '\\mbox{(0.1334$\\pm$0.0033) V}'
#    >> errorValue(.133432, .00332, string_scale=False, units='V')
#    '\\mbox{(1.334$\\pm$0.033)$10^{-1}$ V}'
#    
#    See Also
#    --------
#    copy
#    
#    """
#    
#    # First, I string-format the error to scientific notation with a 
#    # certain number of digits
#    if error_digits >= 1:
#        aux = '{:.' + str(error_digits) + 'E}'
#    else:
#        print("Unvalid 'number_of_digits'! Changed to 1 digit")
#        aux = '{:.0E}'
#    error = aux.format(dX)
#    error = error.split("E") # full error (string)
#    
#    error_order = int(error[1]) # error's order (int)
#    error_value = error[0] # error's value (string)
#
#    # Now I string-format the measure to scientific notation
#    measure = '{:E}'.format(X)
#    measure = measure.split("E") # full measure (string)
#    measure_order = int(measure[1]) # measure's order (int)
#    measure_value = float(measure[0]) # measure's value (string)
#    
#    # Second, I choose the scale I will put both measure and error on
#    # If I want to use the string scale...
#    if -12 <= measure_order < 12 and string_scale:
#        prefix = ['p', 'n', r'$\mu$', 'm', '', 'k', 'M', 'G']
#        scale = [-12, -9, -6, -3, 0, 3, 6, 9, 12]
#        for i in range(8):
#            if not one_point_scale:
#                if scale[i] <= measure_order < scale[i+1]:
#                    prefix = prefix[i] # prefix to the unit
#                    scale = scale[i] # order of both measure and error
#                    break
#            else:
#                if scale[i]-1 <= measure_order < scale[i+1]-1:
#                    prefix = prefix[i]
#                    scale = scale[i]
#                    break
#        used_string_scale = True
#    # ...else, if I don't or can't...
#    else:
#        scale = measure_order
#        prefix = ''
#        used_string_scale = False
#    
#    # Third, I actually scale measure and error according to 'scale'
#    # If error_order is smaller than scale...
#    if error_order < scale:
#        if error_digits - error_order + scale - 1 >= 0:
#            aux = '{:.' + str(error_digits - error_order + scale - 1)
#            aux = aux + 'f}'
#            error_value = aux.format(
#                    float(error_value) * 10**(error_order - scale))
#            measure_value = aux.format(
#                    measure_value * 10**(measure_order - scale))
#        else:
#            error_value = float(error_value) * 10**(error_order - scale)
#            measure_value = float(measure_value)
#            measure_value = measure_value * 10**(measure_order - scale)
#    # Else, if error_order is equal or bigger than scale...
#    else:
#        aux = '{:.' + str(error_digits - 1) + 'f}'
#        error_value = aux.format(
#                float(error_value) * 10**(error_order - scale))
#        measure_value = aux.format(
#                float(measure_value) * 10**(measure_order - scale))
#    
#    # Forth, I make a latex string. Ex.: '(1.34$pm$0.32) kV'
#    if not unpack:
#        latex_str = r'({}$\pm${})'.format(measure_value, error_value)
#        if not used_string_scale and measure_order != 0:
#            latex_str = latex_str + r'$10^{' + '{:.0f}'.format(scale) + '}$'     
#        elif used_string_scale and prefix!='':
#            latex_str = latex_str + ' ' + prefix
#        if prefix!='':
#            latex_str = latex_str + units
#        else:
#            latex_str = latex_str + ' ' + units
#        if units != '' or prefix:
#            if not legend:
#                latex_str = r'\mbox{' + latex_str + '}'
#        if latex_str[-1]==' ':
#            latex_str = latex_str[:len(latex_str)-1]
#                
#    else:
#        latex_str = [r'{}'.format(measure_value), r'{}'.format(error_value)]
#        if not used_string_scale and measure_order != 0:
#            latex_str = [st + r'$10^{' + '{:.0f}'.format(scale) + '}$' 
#                         for st in latex_str]
#        elif used_string_scale and prefix!='':
#            latex_str = [st + ' ' + prefix for st in latex_str]
#        if prefix!='':
#            latex_str = [st + units for st in latex_str]
#        else:
#            latex_str = [st + ' ' + units for st in latex_str]
#        aux = []
#        for st in latex_str:
#            if st[-1]==' ':
#                aux.append(st[:len(st)-1])
#        latex_str = aux
#    
#    return latex_str

"""
>>> errorValue(.133432, .00332, string_scale=False, units='V', unpack=True)
[]
"""