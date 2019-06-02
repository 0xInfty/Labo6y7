# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:08:08 2019

@author: Vall
"""

from iv_save_module import freeFile, loadNicePumpProbe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import os
from tkinter import Tk, messagebox

#%%

def interactiveLegend(ax, labels=False, show_default=True, 
                      location='best', **kwargs):

    """Adds an interactive save button to a given figure.
    
    Parameters
    ----------
    ax : plt.Axes
        The axes to which the interactive legend should be added.
    labels=False : bool, list
        If not false, the list of string names for the different lines that 
        are plotted.
    show_default=True : bool, list
        If not bool, the list of boolean values that say whether to show at 
        first or not the different lines that are plotted.
    location='best' : str
        A string that indicates where to add the legend on the plot area. 
        Can be 'best', 'upper right', 'upper left', 'lower right', 
        'lower left'.
    
    Returns
    -------
    buttons : wid.Button
        Interactive legend instance.
    """
    
    # First, get the lines that are currently plotted
    lines = ax.lines
    if labels is False:
        labels = [l.get_label() for l in lines]

    # Now, if needed, correct labels and default show parameters
    try:
        N = len(labels)
    except:
        N = 1
    if N == 1:
        labels = list(labels)
    try:
        M = len(show_default)
    except:
        M = 1
    if M != N and M == 1:
        show_default = [show_default for l in labels]
    
    # Choose legend location
    number = len(labels)
    height = .05 * number
    extra_y = .05 * (number - 1)
    try:
        fsize = kwargs['fontsize']
    except:
        fsize = 10
    if fsize == 10:
        width = .23
        extra_x = 0
    else:
        width = .23 * fsize / 10
        extra_x = .23 * (fsize/10 - 1)
    try:
        x0 = kwargs.pop('x0')
    except:
        x0 = (.14, .65)
    try:
        y0 = kwargs.pop('y0')
    except:
        y0 = (.03, .81)
    if location=='best':
        xmin = min([min(l.get_data()[0]) for l in lines])
        xmax = max([max(l.get_data()[0]) for l in lines])
        ymin = min([min(l.get_data()[1]) for l in lines])
        ymax = max([max(l.get_data()[1]) for l in lines])
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        if abs(ymin-ylim[0]) > abs(ymax-ylim[1]):
            location = 'lower '
        else:
            location = 'upper '
        if abs(xmin-xlim[0]) > abs(xmax-xlim[1]):
            location = location + 'left'
        else:
            location = location + 'right'
    if location=='upper right':
        position = [x0[1] + extra_x, y0[1] + extra_y, width, height]
    elif location=='upper left':
        position = [x0[0] - extra_x, y0[1] + extra_y, width, height]
    elif location=='lower right':
        position = [x0[0] + extra_x, y0[0] - extra_y, width, height]
    elif location=='lower left':
        position = [x0[1] - extra_x, y0[0] - extra_y, width, height]
    else:
        raise ValueError("Unvalid legend location")
 
    # Create legend buttons
    ax_buttons = plt.axes(position)
    buttons = wid.CheckButtons(ax_buttons, labels, show_default)
    legends = buttons.labels
    for l, leg in zip(lines, legends):
        leg.set_color(l.get_color())
        leg.set(**kwargs)
    for l, sd in zip(lines, show_default):
        l.set_visible(sd)
    
    # Add callback function and run
    def buttons_callback(label):
        for l, leg in zip(lines, legends):
            if label == leg.get_text():
                l.set_visible(not l.get_visible())
        plt.draw()
        return
    buttons.on_clicked(buttons_callback)
    
    plt.show()
    
    return buttons

#%%

def interactiveSaveButton(filename, sufix='', new_format='{}_v{}'):

    """Adds an interactive save button to a given figure.
    
    Parameters
    ----------
    filename : str
        A model filename, which must include full path.
    sufix='' : str
        A sufix to be always added to the given filename.
    new_format='{}_v{}' : str
        A formatter that allows to make new filenames in order to avoid 
        overwriting. If 'F:\Hola.png' does already exist, new file is saved as 
        'F:\Hola_v2.png'.
    
    Returns
    -------
    save_button : wid.Button
        Interactive save button instance.
    """
    
    path, name = os.path.split(filename)
    name = os.path.splitext(name)[0]
    
    # Since I can, I would also like an interactive 'Save' button
    ax_save = plt.axes([0.8, 0.01, 0.1, 0.04])
    save_button = wid.Button(ax_save, 'Guardar')
    
    # For that, I'll need another callback function
    def check_save_callback(event):
        Tk().withdraw()
    #   tk.newfilename = askopenfilename()
        newpath = os.path.join(path, 'Figuras')
        if not os.path.isdir(newpath):
            os.makedirs(newpath)
        newfilename = freeFile(os.path.join(newpath, name+sufix+'.png'),
                               newformat=new_format)
        ax_save.set_visible(False)
        plt.savefig(newfilename, bbox_inches='tight')
        ax_save.set_visible(True)
        messagebox.showinfo('¡Listo!',
                            'Imagen guardada como {}.png'.format(
                    os.path.split(os.path.splitext(newfilename)[0])[-1]))
    save_button.on_clicked(check_save_callback)
    plt.show()
    
    return save_button

#%%

def interactiveValueSelector(ax, x_value=True, y_value=True):
    
    """Allows to choose values from the axes of a plot.
    
    Parameters
    ----------
    ax : plt.Axes
        The axes instance of the plot from where you want to choose.
    x_value=True : bool
        Whether to return the x value.
    y_value=True : bool
        Whether to return the y value.
    
    Returns
    -------
    value : float
        If only one value is required. This is the x value if 'x_value=True' 
        and 'y_value=False'. Otherwise, it is the y value.
    values : tuple
        If both values are required. Then it returns (x value, y value).
    
    See also
    --------
    plt.Axes
    wid.Cursor
    """
    
    ax.autoscale(False)
    cursor = wid.Cursor(ax, useblit=True, 
                        linestyle='--', color='red', linewidth=2)
    if not y_value:
        cursor.horizOn = False
    if not x_value:
        cursor.vertOn = False
    plt.show()
    
    values = plt.ginput()[0]
    if x_value:
        plt.vlines(values[0], ax.get_ylim()[0], ax.get_ylim()[1], 
                   linestyle='--', linewidth=2, color='red')
    if y_value:
        plt.hlines(values[1], ax.get_xlim()[0], ax.get_xlim()[1], 
                   linestyle='--', linewidth=2, color='red')
    cursor.visible = False
    cursor.active = False

    if x_value and y_value:
        return values
    elif x_value:
        return values[0]
    else:
        return values[1]

#%%

def interactiveIntegerSelector(ax, min_value=0, max_value=5):
    
    """Adds an integer selector bar to a single-plot figure.
    
    Allows to choose an integer value looking at a plot.
    
    Parameters
    ----------
    ax : plt.Axes
        The axis instance from the single-plot figure.
    min_value=0 : int
        Minimum integer value that can be chosen.
    max_value=5 : int
        Maximum integer value that can be chosen.
    
    Returns
    -------
    integer : int
        Selected integer value.
    
    See also
    --------
    ivp.IntFillingCursor
    plt.Axes
    """
    
    position = ax.get_position()   
    if ax.xaxis.label.get_text() == '':
        ax.set_position([position.x0,
                         position.y0 + position.height*0.16,
                         position.width,
                         position.height*0.84])
    else:
        ax.set_position([position.x0,
                         position.y0 + position.height*0.18,
                         position.width,
                         position.height*0.82])
    
    ax_selector = plt.axes([0.18, 0.1, 0.65, 0.03])        
    ax_selector.yaxis.set_visible(False)
    ax_selector.set_xlim(min_value, max_value+1)
    selector = IntFillingCursor(ax_selector, color='r', linewidth=2)
    selector.horizOn = False
    plt.show()
    plt.annotate("¿Cantidad?", (0.01, 1.3), xycoords='axes fraction');
    plt.annotate(
            "Elija un número desde {:.0f} hasta {:.0f}.".format(
                    min_value, 
                    max_value), 
            (0.45, 1.3), xycoords='axes fraction');        
    
    integer = int(plt.ginput()[0][0])
    ax_selector.autoscale(False)
    plt.fill([ax_selector.get_xlim()[0], integer,
              integer, ax_selector.get_xlim()[0]],
              [ax_selector.get_ylim()[0], ax_selector.get_ylim()[0],
               ax_selector.get_ylim()[1], ax_selector.get_ylim()[1]],
              'r')
    selector.visible = False
    
    return integer

#%%

class FillingCursor(wid.Cursor):
    
    """Subclass that fills one side of the cursor"""
    
    def __init__(self, ax, horizOn=True, vertOn=True, **lineprops):
        self.fill, = ax.fill([ax.get_xbound()[0], ax.get_xbound()[0],
                             ax.get_xbound()[0], ax.get_xbound()[0]],
                             [ax.get_xbound()[0], ax.get_xbound()[0],
                             ax.get_xbound()[0], ax.get_xbound()[0]],
                             **lineprops)
#        self.fill.set_visible(False)
        self.myax = ax
        super().__init__(ax, horizOn=horizOn, vertOn=vertOn, 
                         useblit=False, **lineprops)
    
    def clear(self, event):
        """Internal event handler to clear the cursor."""
        self.fill.set_visible(False)
        super().clear(event)
    
    def onmove(self, event):
        """Internal event handler to draw the cursor when the mouse moves."""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)
            self.fill.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        self.needclear = True
        if not self.visible:
            return
        self.linev.set_xdata((event.xdata, event.xdata))
        self.lineh.set_ydata((event.ydata, event.ydata))
        if self.vertOn and self.horizOn:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[0],
                                        event.xdata,
                                        event.xdata],
                                        [self.myax.get_ybound()[0],
                                         event.ydata,
                                         event.ydata,
                                         self.myax.get_xbound()[0]]]).T)
        elif self.horizOn:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[1],
                                        self.myax.get_xbound()[1]],
                                        [self.myax.get_ybound()[0],
                                         event.ydata,
                                         event.ydata,
                                         self.myax.get_ybound()[0]]]).T)
        else:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        event.xdata,
                                        event.xdata,
                                        self.myax.get_xbound()[0]],
                                       [self.myax.get_ybound()[0],
                                        self.myax.get_ybound()[0],
                                        self.myax.get_ybound()[1],
                                        self.myax.get_ybound()[1]]]).T)           
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)
        self.fill.set_visible(self.visible and (self.horizOn or self.vertOn))

        self._update()

#%%

class IntFillingCursor(FillingCursor):
    
    """Subclass that only allows integer values on the filling cursor"""
    
    def __init__(self, ax, horizOn=True, vertOn=True,
                 **lineprops):
        super().__init__(ax, horizOn=horizOn, vertOn=vertOn, **lineprops)
        
    def onmove(self, event):
        """Internal event handler to draw the cursor when the mouse moves."""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)
            self.fill.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        self.needclear = True
        if not self.visible:
            return
        self.linev.set_xdata((int(event.xdata), int(event.xdata)))
        self.lineh.set_ydata((int(event.ydata), int(event.ydata)))
        if self.vertOn and self.horizOn:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[0],
                                        int(event.xdata),
                                        int(event.xdata)],
                                        [self.myax.get_ybound()[0],
                                         int(event.ydata),
                                         int(event.ydata),
                                         self.myax.get_xbound()[0]]]).T)
        elif self.horizOn:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[0],
                                        self.myax.get_xbound()[1],
                                        self.myax.get_xbound()[1]],
                                        [self.myax.get_ybound()[0],
                                         int(event.ydata),
                                         int(event.ydata),
                                         self.myax.get_ybound()[0]]]).T)
        else:
            self.fill.set_xy(np.array([[self.myax.get_xbound()[0],
                                        int(event.xdata),
                                        int(event.xdata),
                                        self.myax.get_xbound()[0]],
                                       [self.myax.get_ybound()[0],
                                        self.myax.get_ybound()[0],
                                        self.myax.get_ybound()[1],
                                        self.myax.get_ybound()[1]]]).T)           
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)
        self.fill.set_visible(self.visible and (self.horizOn or self.vertOn))

        self._update()