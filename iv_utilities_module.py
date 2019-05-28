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