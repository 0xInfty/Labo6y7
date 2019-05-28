# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:58:29 2019

@author: 0xInfty
"""

from tkinter import Tk

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