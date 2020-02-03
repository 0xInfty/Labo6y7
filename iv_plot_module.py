# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:08:08 2019

@author: Vall
"""

import iv_save_module as ivs
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import os
from tkinter import Tk, messagebox

#%%

def interactiveLegend(ax, labels=False, show_default=True, 
                      loc='best', **kwargs):

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
    loc='best' : str
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
    if loc=='best':
        xmin = min([min(l.get_data()[0]) for l in lines])
        xmax = max([max(l.get_data()[0]) for l in lines])
        ymin = min([min(l.get_data()[1]) for l in lines])
        ymax = max([max(l.get_data()[1]) for l in lines])
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        if abs(ymin-ylim[0]) > abs(ymax-ylim[1]):
            loc = 'lower '
        else:
            loc = 'upper '
        if abs(xmin-xlim[0]) > abs(xmax-xlim[1]):
            loc = loc + 'left'
        else:
            loc = loc + 'right'
    if loc=='upper right':
        position = [x0[1] - extra_x, y0[1] - extra_y, width, height]
    elif loc=='upper left':
        position = [x0[0] + extra_x, y0[1] - extra_y, width, height]
    elif loc=='lower right':
        position = [x0[1] - extra_x, y0[0] + extra_y, width, height]
    elif loc=='lower left':
        position = [x0[0] + extra_x, y0[0] + extra_y, width, height]
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

def interactiveSaveButton(filename, **kwargs):

    """Adds an interactive save button to a given figure.
    
    Parameters
    ----------
    filename : str
        A model filename, which must include full path.
    
    Other parameters
    ----------------
    overwrite=False : bool
        Says whether to overwrite files or not.
    sufix='' : str
        A sufix to be always added to the given filename.
    newformat='{}_v{}' : str
        A formatter that allows to make new filenames in order to avoid 
        overwriting. If 'F:\Hola.png' does already exist, new file is saved as 
        'F:\Hola_v2.png'.
    
    Returns
    -------
    save_button : wid.Button
        Interactive save button instance.
    """
    
    # Since I can, I would also like an interactive 'Save' button
    ax_save = plt.axes([0.8, 0.01, 0.1, 0.04])
    save_button = wid.Button(ax_save, 'Guardar')
    
    # For that, I'll need another callback function
    def check_save_callback(event):
        Tk().withdraw()
    #   tk.newfilename = askopenfilename()
        ax_save.set_visible(False)
        ivs.saveFig(filename, **kwargs)
        ax_save.set_visible(True)
        messagebox.showinfo('¡Listo!', 'Imagen guardada')
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
 
def interactiveTimeSelector(filename, autoclose=True):
    
    """Allows to select a particular time instant on a Pump Probe file.
    
    Parameters
    ----------
    filename : str
        Filename, which must include full path and extension.
    autoclose=True : bool
        Says whether to automatically close this picture or not.
    
    Returns
    -------
    ti : float
        Selected value.
    
    See also
    --------
    ivs.loadNicePumpProbe
    """
    
    t, V, details = ivs.loadNicePumpProbe(filename)
    fig = plotPumpProbe(filename, autosave=False)[0]
    ax = fig.axes[0]
    ti = interactiveValueSelector(ax, y_value=False)
    ti = t[np.argmin(abs(t-ti))]
    
    if autoclose:
        plt.close(fig)

    return ti

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

#%%

def plotPumpProbe(filename, extension='.png', interactive=False, autosave=True, 
                  overwrite=False, **kwargs):

    """Plots a PumpProbe experiment from a file and its mean.
    
    Can also make an interactive plot, which holds a save button and allows to 
    choose only certain experiments to be shown from the legend.
    
    By default, it also saves a picture on the file's path.
        
    Parameters
    ----------
    filename : str
        File's root (must include directory and termination).
    extension='.png' : str
        Image file's format.
    interactive=True : bool
        Says whether to make an interactive plot or not.
    autosave=True : bool
        Says whether to automatically save or not.
    overwrite=False : bool
        Says whether to allow overwriting or not.
    
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
    
    t, V, details = ivs.loadNicePumpProbe(filename)
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
        legend_buttons = interactiveLegend(ax, labels, show_default, 
                                           fontsize=12,
                                           x0=(.17, .68), y0=(.06, .84), **kwargs)
        save_button = interactiveSaveButton(filename, extension=extension, 
                                            overwrite=overwrite,
                                            folder='Figuras',
                                            sufix='_fig')
    else:
        plt.legend(labels, fontsize=12, framealpha=1, **kwargs)
       
    if autosave:
        if interactive:
            save_button.ax.set_visible(False)
        save_kwargs = dict()
        if 'newformat' in kwargs.keys():
            save_kwargs.add('newformat', kwargs['newformat'])        
        ivs.saveFig(filename, extension=extension, overwrite=overwrite, 
                    folder='Figuras', sufix='_fig')
        if interactive:
            save_button.ax.set_visible(True)
    
    if interactive:
        return fig, legend_buttons, save_button
    else:
        return fig, None, None

#%%

def plotAllPumpProbe(path, extension='.png', autosave=True, autoclose=False, 
                     **kwargs):
    
    """Plots all PumpProbe experiments on the files from a given path.
        
    The data files must be '.txt' files that begin with 'M'.
    
    Parameters
    ----------
    path : str
        Files' folder (must include directory).
    extension='.png' : str
        Image file's format.
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
        fig = plotPumpProbe(f, extension=extension, interactive=False, 
                            autosave=autosave, **kwargs)[0]
        if autoclose:
            plt.close(fig)
        else:
            figures.append(fig)
    
    if not autoclose:
        return figures

#%%
    
def linearPredictionPlot(filename, plot_results, extension='.png', 
                         folder='Figuras', autosave=True, overwrite=False,
                         showgrid=False):

    """Plots the results of a linear prediction plot.
    
    Parameters
    ----------
    filename : str
        File's root (must include directory and extension).
    plot_results : ivu.InstancesDict
        Fit results that allow to plot. Must include...
        ...numpy array 'fit', that holds time, data, fit and fit terms
        ...numpy.array 'raman', that holds frequencies, fit spectrum and fit 
        terms' spectrum.
    extension='.png' : str
        Image file's format.
    folder='Figuras' : str
        Folder to include in figure's filename.
    autosave=True : bool
        Says whether to save or not.
    overwrite=False : bool
        Says whether to allow overwriting or not.
    showgrid=False : bool
        Says whether to show or not the vertical grid on the time space plot.
    
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
    Image file. Only raised if 'autosave=True'.
    
    See also
    --------
    iva.linearPrediction
    
    """
    
    # First I deglose data
    fit = plot_results.fit
    raman = plot_results.raman
    Nfit_terms = fit.shape[1] - 3
    
    # In order to save, if needed, I will need...
    filename = os.path.splitext(filename)[0] + extension
    
    # Then, to plot, I first start a figure
    fig = plt.figure()
    grid = plt.GridSpec(3, 5, hspace=0.1)
    
    # In the upper subplot, I put the Raman-like spectrum
    ax_spectrum = plt.subplot(grid[0,:4])
    plt.plot(raman[:,0], raman[:,1], linewidth=2)
    lspectrum_terms = plt.plot(raman[:,0], raman[:,2:], 
                               linewidth=2)
    for l in lspectrum_terms: l.set_visible(False)
    plt.xlabel("Frecuencia (GHz)")
    plt.ylabel("Amplitud (u.a.)")
    ax_spectrum.xaxis.tick_top()
    ax_spectrum.xaxis.set_label_position('top')
    
    # In the lower subplot, I put the data and fit
    ax_data = plt.subplot(grid[1:,:])
    ldata, = plt.plot(fit[:,0], fit[:,1], 'k', linewidth=0.4)
#    ax_data.autoscale(False)
    lfit, = plt.plot(fit[:,0], fit[:,2], linewidth=2)
    lfit_terms = plt.plot(fit[:,0], fit[:,3:], linewidth=1)
    for l in lfit_terms: l.set_visible(False)
    plt.xlabel("Tiempo (ps)")
    plt.ylabel(r"Voltaje ($\mu$V)")
    ax_data.tick_params(labelsize=12)
    if showgrid:
        ax_data.minorticks_on()
        ax_data.tick_params(axis='y', which='minor', left=False)
        ax_data.tick_params(length=5)
        ax_data.grid(axis='x', which='both')
        ldata.set_linewidth(0.6)
        lfit.set_linewidth(2.3)
    
    # Because it's pretty, I make an interactive legend
    ax_legend = plt.axes([0.75, 0.642, 0.155, 0.24])
    legend_buttons = wid.CheckButtons(ax_legend, ('Data', 
                                     'Ajuste', 
                                     *['Término {:.0f}'.format(i+1) 
                                     for i in range(Nfit_terms)]), 
                        (True, True, *[False for i in range(Nfit_terms)]))
    legend_buttons.labels[1].set_color(lfit.get_color())
    for leg, l in zip(legend_buttons.labels[2:], lfit_terms):
        leg.set_color(l.get_color())
    
    # For that, I'll need a callback function  
    def legend_callback(label):
        if label == 'Data':
            ldata.set_visible(not ldata.get_visible())
        elif label == 'Ajuste':
            lfit.set_visible(not lfit.get_visible())
        else:
            for i in range(Nfit_terms):
                if label == 'Término {:.0f}'.format(i+1):
                    lfit_terms[i].set_visible(not lfit_terms[i].get_visible())
                    lspectrum_terms[i].set_visible(
                            not lspectrum_terms[i].get_visible())
        plt.draw()
    legend_buttons.on_clicked(legend_callback)
    
    # Since I can, I would also like an interactive 'Save' button
    save_button = interactiveSaveButton(filename, overwrite=overwrite, 
                                        folder=folder, sufix='_fit')
    
    # Once I have all that, I'll show the plot
    plt.show()
    
    # Like it is shown for the first time, autosave if configured that way
    if autosave:
        save_button.ax.set_visible(False)
        ivs.saveFig(filename, overwrite=overwrite, folder=folder, 
                    sufix='_fit')
        save_button.ax.set_visible(True)
        
    return fig, legend_buttons, save_button