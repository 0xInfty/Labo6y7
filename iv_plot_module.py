# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:08:08 2019

@author: LEC
"""

from iv_save_module import freeFile, loadNicePumpProbe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as wid
import os
from tkinter import Tk, messagebox

#%%

def interactiveValueSelector(ax, x_value=True, y_value=True):
    
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

def plotPumpProbe(filename, save=True):

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
    
    path = os.path.join(os.path.split(filename)[0], 'Figuras')
    name = os.path.split(os.path.splitext(filename)[0])[-1]
    t, V, meanV, details = loadNicePumpProbe(filename)
    Nrepetitions = details['Nrepetitions']
    
    fig = plt.figure()
    plt.plot(t, V, linewidth=0.8)
    plt.plot(t, meanV, linewidth=1.5)
    legends = ['Experimento {:.0f}'.format(i+1) for i in range(Nrepetitions)]
    legends.append('Promedio')
    plt.legend(legends)
    plt.ylabel(r'Voltaje ($\mu$V)')
    plt.xlabel(r'Tiempo (ps)')
    
    if not os.path.isdir(path):
        os.makedirs(path)
    
    if save:
        plt.savefig(os.path.join(path,name+'_fig.png'), bbox_inches='tight')
    
    return fig

#%%

def fullplotPumpProbe(filename, save=True):
    
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
    
    path = os.path.join(os.path.split(filename)[0], 'Figuras')
    name = os.path.split(os.path.splitext(filename)[0])[-1]
    t, V, meanV, details = loadNicePumpProbe(filename)
    Nrepetitions = details['Nrepetitions']
    
    fig = plt.figure()
    grid = plt.GridSpec(Nrepetitions, 2, wspace=0.4, hspace=0.3)
    
    for i, v in enumerate(V.T):
        plt.subplot(Nrepetitions, 2, 2*i+1)
        plt.plot(t, V, linewidth=0.8)
    
    plt.subplot(grid[0,1])
    plt.plot(t, meanV, 'r', linewidth=0.8)
    
    plt.subplot(grid[1:,1])
    plt.plot(t, V, linewidth=0.8)
    plt.plot(t, meanV, linewidth=1.5)

    if not os.path.isdir(path):
        os.makedirs(path)
    
    if save:
        plt.savefig(os.path.join(path,name+'_full.png'), bbox_inches='tight')
    
    return fig

#%%

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

#%%
    
def linearPredictionPlot(filename, others, autosave=True):

    # First I deglose data
    fit = others['fit']
    raman = others['raman']
    Nfit_terms = fit.shape[1] - 3
    
    # In order to save, if needed, I will need...
    name = os.path.split(os.path.splitext(filename)[0])[-1]
    path = os.path.split(filename)[0]
    
    # Then, to plot, I first start a figure
    plt.figure()
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
    ldata, = plt.plot(fit[:,0], fit[:,1], 'k', linewidth=0.5)
    ax_data.autoscale(False)
    lfit, = plt.plot(fit[:,0], fit[:,2], linewidth=2)
    lfit_terms = plt.plot(fit[:,0], fit[:,3:], linewidth=2)
    for l in lfit_terms: l.set_visible(False)
    plt.xlabel("Tiempo (ps)")
    plt.ylabel(r"Voltaje ($\mu$V)")
    
    # Because it's pretty, I make an interactive legend
    ax_legend = plt.axes([0.75, 0.642, 0.155, 0.24])
    check_legend = wid.CheckButtons(ax_legend, ('Data', 
                                     'Ajuste', 
                                     *['Término {:.0f}'.format(i+1) 
                                     for i in range(Nfit_terms)]), 
                        (True, True, *[False for i in range(Nfit_terms)]))
    check_legend.labels[1].set_color(lfit.get_color())
    for leg, l in zip(check_legend.labels[2:], lfit_terms):
        leg.set_color(l.get_color())
    
    # For that, I'll need a callback function  
    def check_legend_callback(label):
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
    check_legend.on_clicked(check_legend_callback)
    
    # Since I can, I would also like an interactive 'Save' button
    ax_save = plt.axes([0.8, 0.01, 0.1, 0.04])
    check_save = wid.Button(ax_save, 'Guardar')
    
    # For that, I'll need another callback function
    def check_save_callback(event):
        Tk().withdraw()
    #   tk.newfilename = askopenfilename()
        newpath = os.path.join(path, 'Figuras')
        if not os.path.isdir(newpath):
            os.makedirs(newpath)
        newfilename = freeFile(os.path.join(newpath, name+'_fit.png'),
                               newformat='{}__{}')
        ax_save.set_visible(False)
        plt.savefig(newfilename, bbox_inches='tight')
        ax_save.set_visible(True)
        messagebox.showinfo('¡Listo!',
                            'Imagen guardada como {}.png'.format(
                    os.path.split(os.path.splitext(newfilename)[0])[-1]))
    check_save.on_clicked(check_save_callback)
    
    # Once I have all that, I'll show the plot
    plt.show()
    
    # Like it is shown for the first time, autosave if configured that way
    if autosave:
        newpath = os.path.join(path, 'Figuras')
        if not os.path.isdir(newpath):
            os.makedirs(newpath)
        ax_save.set_visible(False)
        plt.savefig(os.path.join(newpath, name+'_fit.png'), bbox_inches='tight')
        ax_save.set_visible(True)
        
    return