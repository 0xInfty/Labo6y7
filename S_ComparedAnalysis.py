# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 09:58:43 2019

@author: Vall
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import scipy.stats as st
import iv_save_module as ivs
import iv_utilities_module as ivu
import iv_analysis_module as iva

#%% PARAMETERS ----------------------------------------------------------------

# Saving parameters
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
figs_extension = '.png'
overwrite = True

# Analysis Parameters
load_sem = True
filter_not_in_common_rods = True
filter_air_outliers = False

# For each data to compare, we need one value on each list
desired_frequency = [12, 16]#, 8] # in GHz
full_series = ['Fused Silica + Aire', 
               'Fused Silica + Ta2O5']#,
#               'Ta2O5 + Aire']
series = ['LIGO1', 'LIGO1_PostUSA']#, 'LIGO5bis']
sem_series = ['LIGO1_1', 'LIGO1_1']#, 'LIGO5bis_1'] # ['nan']
sem_full_series = ['L1 1', 'L1 1']#, 'L5bis 1']
name = 'FusedSilica' # 'All'

#%% OTHERS PARAMETERS ---------------------------------------------------------

# Some more strings
filter_conditions = ''
if filter_not_in_common_rods:
    filter_conditions += '_CommonRods'
elif filter_air_outliers:
    filter_conditions += '_NoOutliers'

# Some functions and variables to manege filenames
def semFilename(sem_series, home=home):
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    filename = 'Resultados_SEM_{}.txt'.format(sem_series)
    sem_series = sem_series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *sem_series, filename)
    return sem_filename
rodsFilename = lambda series : os.path.join(home, 
                                           r'Análisis\Rods_{}.txt'.format(
                                                   series))
paramsFilename = lambda series : os.path.join(home, 
                                              r'Análisis\Params_{}.txt'.format(
                                                   series))
figsFilename = lambda fig_name : os.path.join(home, fig_name+'.png')
figs_folder = r'Análisis\ComparedAnalysis{}_{}'.format(filter_conditions,
                                                       name)

#%% PHYSICS -------------------------------------------------------------------

# Physics' Parameters
physics = {}
physics['density_gold'] = 19.3e3 # kg/m3 for gold
physics['shear_fsilica'] = np.mean([30.8e9, 32.3e9]) # Pa for fused silica
physics['viscosity_gold'] = 2e-3 # Pa/s for gold
physics['young_gold'] = np.mean([71.2e9, 74.8e9])  # Pa for fused silica
physics['density_fsilica'] = np.mean([2.17e3, 2.22e3]) # kg/m3 for fused silica
physics['density_ta2o5'] = 8.18e3 # kg/m3 for Ta2O5

#%% LOAD DATA -----------------------------------------------------------------

filenames = []
rods = []
params = []
fits_data = []
fits_footer = []
frequency = []
damping_time = []
quality_factor = []
if load_sem:
    sem_data = []
    sem_footer = []
    length = []
    width = []

for s, ss, f in zip(series, sem_series, desired_frequency):

    # Look for the list of rods and filenames
    sfilenames = [] # Will contain filenames like 'M_20190610_01'
    srods = [] # Will contain rods' positions like '1,2'
    with open(rodsFilename(s), 'r') as file:
        for line in file:
            if line[0]!='#':
                sfilenames.append(line.split('\t')[0]) # Save filenames
                aux = line.split('\t')[1:]
                aux = r' '.join(aux)
                srods.append(aux.split('\n')[0]) # Save rods
                del aux
        del line
    
    # Then load parameters
    params_filenames = [] # Will contain filenames like 'M_20190610_01'
    amplitude = []
    power = []
    wavelength = []
    spectral_width = []
    with open(paramsFilename(s), 'r') as file:
        for line in file:
            if line[0]!='#':
                params_filenames.append(line.split('\t')[0])
                amplitude.append(float(line.split('\t')[1]))
                power.append(float(line.split('\t')[2]))
                wavelength.append(float(line.split('\t')[3]))
                spectral_width.append(float(line.split('\t')[-1]))
        del line
    sparams = np.array([amplitude, power, wavelength, spectral_width]).T
    index = [params_filenames.index(f) for f in sfilenames]
    sparams = sparams[index,:]
    params_header = ['Amplitud (mVpp)', 'Potencia Pump post-MOA (muW)', 
                     'Longitud de onda (nm)', 'Ancho medio de la campana (nm)']
    del params_filenames, index, amplitude, power, wavelength, spectral_width
    
    # Now create a list of folders for each filename    
    fits_filenames = [ivs.filenameToFitsFilename(file, home) for file in sfilenames]
     
    # Load data from each fit
    sfits_data = []
    sfits_footer = []
    for file in fits_filenames:
        data, fits_header, footer = ivs.loadTxt(file)
        sfits_data.append(data)
        sfits_footer.append(footer)
    del file, data, footer, fits_filenames
    
    # Keep only the fit term that has the closest frequency to the desired one
    fits_new_data = []
    for rod, fit in zip(srods, sfits_data):
        try:
            i = np.argmin(abs(fit[:,0] - f*np.ones(fit.shape[0])))
            fits_new_data.append([*fit[i,:]])
        except IndexError:
            fits_new_data.append([*fit])
    sfits_data = np.array(fits_new_data)
    sfrequency = sfits_data[:,0]*1e9 # Hz
    sdamping_time = sfits_data[:,1]*1e-12 # s
    squality_factor = sfits_data[:,2]
    del rod, fit, i, fits_new_data
    
    if load_sem:
    
        # Also load data from SEM dimension analysis
        ssem_data, sem_header, ssem_footer = ivs.loadTxt(semFilename(ss))
          
        # Now lets put every rod in the same order for SEM and fits
        index = [ssem_footer['rods'].index(r) for r in srods]
        ssem_data = ssem_data[index,:]
        slength = ssem_data[:,2] * 1e-9 # m
        swidth = ssem_data[:,0] * 1e-9 # m
        del index
    
    # Now we can filter the results
    if s == 'LIGO1' and filter_air_outliers:
        index = np.argsort(sfrequency)[2:] # Remove the two lowest frequencies
        sfilenames = [sfilenames[k] for k in index]
        srods = [srods[k] for k in index]
        if load_sem:
            ssem_data = ssem_data[index,:]
            slength = slength[index]
            swidth = swidth[index]
        sparams = sparams[index]
        sfits_data = sfits_data[index,:]
        sfits_footer = [sfits_footer[k] for k in index]
        sfrequency = sfrequency[index]
        sdamping_time = sdamping_time[index]
        squality_factor = squality_factor[index]
        del index
    
    # Since I'll be analysing frequency vs length mostly...
    if load_sem:
        index = list(np.argsort(slength))
        sfilenames = [sfilenames[k] for k in index]
        srods = [srods[k] for k in index]
        ssem_data = ssem_data[index,:]
        slength = slength[index]
        swidth = swidth[index]
        sparams = sparams[index]
        sfits_data = sfits_data[index,:]
        sfits_footer = [sfits_footer[k] for k in index]
        sfrequency = sfrequency[index]
        sdamping_time = sdamping_time[index]
        squality_factor = squality_factor[index]
        del index
    
    # Now add all that data to a list outside the loop
    filenames.append(sfilenames)
    rods.append(srods)
    params.append(sparams)
    fits_data.append(sfits_data)
    fits_footer.append(sfits_footer)
    frequency.append(sfrequency)
    damping_time.append(sdamping_time)
    quality_factor.append(squality_factor)
    if load_sem:
        sem_data.append(ssem_data)
        length.append(slength)
        width.append(swidth)
    
del s, ss, f
del sfilenames, srods, sfits_data, sfits_footer, sfrequency, sdamping_time
del squality_factor, sparams
if load_sem:
    del ssem_data, ssem_footer, slength, swidth

# Now lets discard rods that aren't in all of the samples
if filter_not_in_common_rods:
    
    remove_rods = []
    for j in range(len(rods)):
        for j2 in range(len(rods)):
            if j2 != j:
                for r in rods[j]:
                    if r not in rods[j2]:
                        remove_rods.append(r)
    del j, j2

    nfilenames = []
    nrods = []
    if load_sem:
        nsem_data = []
        nlength = []
        nwidth = []
    nparams = []
    nfits_data = []
    nfits_footer = []
    nfrequency = []
    ndamping_time = []
    nquality_factor = []
    for j in range(len(rods)):
        index = []
        for r in rods[j]:
            if r not in remove_rods:
                index.append(rods[j].index(r))
        nfilenames.append([filenames[j][k] for k in index])
        nrods.append([rods[j][k] for k in index])
        if load_sem:
            nsem_data.append(sem_data[j][index,:])
            nlength.append(length[j][index])
            nwidth.append(width[j][index])
        nparams.append(params[j][index])
        nfits_data.append(fits_data[j][index,:])
        nfits_footer.append([fits_footer[j][k] for k in index])
        nfrequency.append(frequency[j][index])
        ndamping_time.append(damping_time[j][index])
        nquality_factor.append(quality_factor[j][index])
        del index
    del j
    
    filenames = nfilenames
    rods = nrods
    if load_sem:
        sem_data = nsem_data
        length = nlength
        width = nwidth
    params = nparams
    fits_data = nfits_data
    fits_footer = nfits_footer
    frequency = nfrequency
    damping_time = ndamping_time
    quality_factor = nquality_factor
    del nfilenames, nrods, nparams, nfits_data, nfits_footer, nfrequency
    del ndamping_time, nquality_factor
    if load_sem:
        del nsem_data, nlength, nwidth

#%% SAVE DATA -----------------------------------------------------------------

if filter_not_in_common_rods:

    # Make OneNote table
    heading = '\t'.join(["Rod", "Longitud (nm)", "Error (nm)", 
                         *["Frecuencia {} (GHz)".format(fs) for fs in full_series], 
                         *["Factor de calidad {} (GHz)".format(fs) 
                            for fs in full_series]])
    items = []
    for r in range(len(rods[0])):
        
        h = '\t'.join(ivu.errorValue(sem_data[0][r,2], sem_data[0][r,3]))
        auxf = []
        for j in range(len(full_series)):
            auxf.append('{:.2f}'.format(fits_data[j][r,0]))
        auxf = '\t'.join(auxf)
        auxt = []
        for j in range(len(full_series)):
            auxt.append('{:.2f}'.format(fits_data[j][r,1]))
        auxt = '\t'.join(auxt)
        items.append('\t'.join([h, auxf, auxt]))
    del h, auxf, auxt
    items = ['\t'.join([ri, i]) for ri, i in zip(rods[j], items)]
    items = '\n'.join(items)
    table = '\n'.join([heading, items])
    ivu.copy(table)
    
    # Save all important data to a single file
    whole_filename = os.path.join(home, figs_folder, 'Resultados_Comparados.txt')
    whole_data = np.array([*sem_data[0][:,:6].T, 
                           *[fits_data[j][:,0] for j in range(len(full_series))],
                           *[fits_data[j][:,1] for j in range(len(full_series))], 
                           *[fits_data[j][:,2] for j in range(len(full_series))]
                           ])
    ivs.saveTxt(whole_filename, whole_data.T, 
                header=["Ancho (nm)", "Error (nm)",
                        "Longitud (nm)", "Error (nm)", 
                        "Relación de aspecto", "Error",
                        *["Frecuencia {} (GHz)".format(fs) for fs in full_series], 
                        *["Tiempo de decaimiento {} (GHz)".format(fs) 
                            for fs in full_series],
                        *["Factor de calidad {} (GHz)".format(fs) 
                            for fs in full_series]],
                footer=dict(rods=rods, filenames=filenames),
                overwrite=True)
    del whole_data, whole_filename

else:

    # Make OneNote tables
    tables = []
    for j, fs, ss in zip(range(len(series)), full_series, sem_full_series):
        heading = '\t'.join(["Rod {}".format(ss), 
                             "Longitud {} (nm)".format(ss), 
                             "Error {} (nm)".format(ss), 
                             "Frecuencia {} (GHz)".format(fs), 
                             "Tiempo de decaimiento {} (ps)".format(fs),
                             "Factor de calidad {} (GHz)".format(fs)])
        items = []
        for r in range(len(rods[j])):
            
            h = '\t'.join(ivu.errorValue(sem_data[j][r,2], sem_data[j][r,3]))
            f = '{:.2f}'.format(fits_data[j][r,0])
            t = '{:.2f}'.format(fits_data[j][r,1])
            q = '{:.2f}'.format(fits_data[j][r,2])
            items.append('\t'.join([h, f, t, q]))
        del h, f, t, q
        items = ['\t'.join([ri, i]) for ri, i in zip(rods[j], items)]
        items = '\n'.join(items)
        tables.append('\n'.join([heading, items]))
    del items, heading, r, j
        
    # Save all important data to a single file
    for j, s, fs, ss in zip(range(len(series)), series, 
                            full_series, sem_full_series):
        whole_filename = os.path.join(home, figs_folder, 
                                      'Resultados_Comparados_{}.txt'.format(s))
        whole_data = np.array([*sem_data[j][:,:6].T, 
                               *fits_data[j][:,:3].T])
        ivs.saveTxt(whole_filename, whole_data.T, 
                    header=["Ancho {} (nm)".format(ss), "Error (nm)",
                            "Longitud {} (nm)".format(ss), "Error (nm)", 
                            "Relación de aspecto {}".format(ss), "Error",
                            "Frecuencia {} (GHz)".format(fs), 
                            "Tiempo de decaimiento {} (GHz)".format(fs),
                            "Factor de calidad {} (GHz)".format(fs)
                            ],
                    footer=dict(rods=rods[j], filenames=filenames[j]),
                    overwrite=True)
    del whole_data, whole_filename, j, s, fs, ss
    
#%% *) FREQUENCY PER ROD

if filter_not_in_common_rods:
    
    # Plot results for the different rods
    fig, ax1 = plt.subplots()
    
    # Frequency plot, right axis
    ax1.set_xlabel('Antena')
    ax1.set_ylabel('Frecuencia (GHz)')
    for j in range(len(series)):
        ax1.plot(fits_data[j][:,0], 'o')
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    ax1.legend(full_series)
    plt.show()
    
    # Format graph
    plt.xticks(np.arange(len(rods[0])), rods[0], rotation='vertical')
    plt.grid(which='both', axis='x')
    ax1.tick_params(length=2)
    ax1.grid(axis='x', which='minor')
    ax1.tick_params(axis='x')#, labelrotation=90)
    plt.show()
    
    # Save plot
    ivs.saveFig(figsFilename('FvsRod'), extension=figs_extension, 
                folder=figs_folder, overwrite=overwrite)

else:

    for j, s, fs in zip(range(len(series)), series, full_series):
    
        # Plot results for the different rods
        fig, ax1 = plt.subplots()
        
        # Frequency plot, right axis
        ax1.set_xlabel('Antena')
        ax1.set_ylabel('Frecuencia (GHz)')
        ax1.plot(fits_data[j][:,0], 'o')
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        ax1.legend([fs])
        plt.show()
        
        # Format graph
        plt.xticks(np.arange(len(rods[j])), rods[j], rotation='vertical')
        plt.grid(which='both', axis='x')
        ax1.tick_params(length=2)
        ax1.grid(axis='x', which='minor')
        ax1.tick_params(axis='x')#, labelrotation=90)
        plt.show()
        
        # Save plot
        ivs.saveFig(figsFilename('FvsRod_{}'.format(s)), 
                    extension=figs_extension, 
                    folder=figs_folder, overwrite=overwrite)

#%% *) BOXPLOTS

fig = plt.figure()
grid = plt.GridSpec(1, 2, wspace=0.5, hspace=0)

ax = plt.subplot(grid[0,0])
ax.boxplot([fits_data[j][:,0] for j in range(len(series))], 
           showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2},
           labels=full_series)
plt.ylabel("Frecuencia (GHz)")
ax.tick_params(axis='y', direction='in')
    
ax = plt.subplot(grid[0,1])
ax.boxplot([fits_data[j][:,2] for j in range(len(series))], 
           showmeans=True, meanline=True, 
           meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
           medianprops={'color':'r', 'linewidth':2},
           labels=full_series)
plt.ylabel("Factor de calidad")
ax.tick_params(axis='y', direction='in')

ivs.saveFig(figsFilename('Boxplots'), extension=figs_extension, 
            folder=figs_folder, overwrite=overwrite)

#%% 

"""UP TO THERE IT'S ALL GENERIC. FROM HERE ON... JUST FUSED SILICA + AIR (0)
VS FUSED SILICA + TA2O5 (1)"""

#%% *) FREQUENCY ON SAME RODS

nbins = 10
set_bigger_percent = 20

if filter_not_in_common_rods:
    
    fig = plt.figure()
    grid = plt.GridSpec(4, 5, wspace=0, hspace=0)
#    fig.set_figheight(fig.get_figwidth())
    
    ax = plt.subplot(grid[1:,:-1])
    ax.plot(frequency[0]*1e-9, frequency[1]*1e-9, 'xk', markersize=7)
    plt.xlabel("Frecuencia {} (GHz)".format(full_series[0]))
    plt.ylabel("Frecuencia {} (GHz)".format(full_series[1]))
    
    # Grid's format
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which='major', axis='both')
    ax.grid(which='minor', axis='both', linestyle=':')
    
    # Limits' format
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    delta_xlims = xlims[1] - xlims[0]
    delta_ylims = ylims[1] - ylims[0]
    new_xlims = (xlims[0] - delta_xlims*set_bigger_percent/100,
                 xlims[1] + delta_xlims*set_bigger_percent/100)
    new_ylims = (ylims[0] - delta_ylims*set_bigger_percent/100,
                 ylims[1] + delta_ylims*set_bigger_percent/100)
    ax.set_xlim(new_xlims)
    ax.set_ylim(new_ylims)
    new_lims = [new_xlims, new_ylims]
    new_lims_T = [new_ylims, new_xlims]

    # Mean values
    line_functions = [plt.vlines, plt.hlines]
    colors = ['blue', 'red']
    for i in range(2):
        line_functions[i](np.mean(frequency[i])*1e-9, 
                          *new_lims_T[i], colors=colors[i], linestyle='--')
    del i
    
    # Standard deviation
    fill_function = [ax.fill_betweenx, ax.fill_between]
    for i in range(2):
        fill_function[i](new_lims_T[i], 
                         (np.mean(frequency[i])-np.std(frequency[i]))*1e-9,
                         (np.mean(frequency[i])+np.std(frequency[i]))*1e-9,
                         color=colors[i],
                         alpha=0.1)
    
    # Histograms
    axh = []
    limsh = []
    grid_places = [grid[0,:-1], grid[1:,-1]]
    orientations = ['vertical', 'horizontal']
    function_lims = [plt.xlim, plt.ylim]
    function_lims_T = [plt.ylim, plt.xlim]
    frequency_linspaces = [np.linspace(nl[0], nl[1], 50) for nl in new_lims]
    normal_distributions = [st.norm.pdf(flins, np.mean(f)*1e-9, np.std(f)*1e-9)
                            for f, flins in zip(frequency, frequency_linspaces)]
    normal_pairs = [[flins, ndist] for flins, ndist 
                    in zip(frequency_linspaces, normal_distributions)]
    normal_pairs[1].reverse()
    for i in range(2):
        axh.append(plt.subplot(grid_places[i]))        
        # Histogram
        n, b, p = axh[i].hist(frequency[i]*1e-9, nbins, density=True,
                              alpha=0.4, facecolor=colors[i],                               orientation=orientations[i])
        # Curve over histogram
        axh[i].plot(*normal_pairs[i], color=colors[i])    
        # Format
        axh[i].axis('off')
        function_lims[i](new_lims[i])
        limsh.append(function_lims_T[i]())
        # Mean values
        line_functions[i](np.mean(frequency[i])*1e-9, *limsh[i], 
                          colors=colors[i], linestyle='--')     
    del i

    ivs.saveFig(figsFilename('FvsF0'), extension=figs_extension, 
            folder=figs_folder, overwrite=overwrite)

#%% *) FREQUENCY ON SAME RODS - ANDREA'S

nbins = 10
set_bigger_percent = 20

if filter_not_in_common_rods:
    
    fig = plt.figure()
    grid = plt.GridSpec(4, 5, wspace=0, hspace=0)
#    fig.set_figheight(fig.get_figwidth())
    
    ax = plt.subplot(grid[1:,:-1])
    ax.plot(frequency[0]*1e-9, frequency[1]*1e-9, 'ko', markersize=8, mfc='w')
    plt.xlabel("Frecuencia {} (GHz)".format(full_series[0]))
    plt.ylabel("Frecuencia {} (GHz)".format(full_series[1]))
    
    # Grid's format
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which='major', axis='both')
    ax.grid(which='minor', axis='both', linestyle=':')
    
    # Limits' format
    xlims = ax.get_xlim() # GHz
    ylims = ax.get_ylim() # GHz
    delta_xlims = xlims[1] - xlims[0]
    delta_ylims = ylims[1] - ylims[0]
    new_xlims = (xlims[0] - delta_xlims*set_bigger_percent/100,
                 xlims[1] + delta_xlims*set_bigger_percent/100)
    new_ylims = (ylims[0] - delta_ylims*set_bigger_percent/100,
                 ylims[1] + delta_ylims*set_bigger_percent/100)
    ax.set_xlim(new_xlims)
    ax.set_ylim(new_ylims)
    new_lims = [new_xlims, new_ylims]
    new_lims_T = [new_ylims, new_xlims]
    
    # Identity
    frequency_linspaces = [np.linspace(nl[0], nl[1], 50) for nl in new_lims]
    ax.plot(frequency_linspaces[0], frequency_linspaces[0], '-k')
    
    # Identity with mean difference vertical shift
    delta_frequency = frequency[1] - frequency[0] # Hz
    ax.plot(frequency_linspaces[0], 
            frequency_linspaces[0] + np.mean(delta_frequency)*1e-9,
            '--k')
    
    # Identity with mean difference standard deviation vertical shift
    ax.fill_between(
        frequency_linspaces[0], 
        frequency_linspaces[0] + (np.mean(delta_frequency) 
                - np.std(delta_frequency))*1e-9,
        frequency_linspaces[0] + (np.mean(delta_frequency) 
                + np.std(delta_frequency))*1e-9,
        color='m',
        alpha=0.2)
    
    # Mean values
    line_functions = [plt.vlines, plt.hlines]
    colors = ['blue', 'red']
    for i in range(2):
        line_functions[i](np.mean(frequency[i])*1e-9, 
                          *new_lims_T[i], colors=colors[i], linestyle='--')
    del i
    
    # Histograms
    axh = []
    limsh = []
    grid_places = [grid[0,:-1], grid[1:,-1]]
    orientations = ['vertical', 'horizontal']
    function_lims = [plt.xlim, plt.ylim]
    function_lims_T = [plt.ylim, plt.xlim]
    normal_distributions = [st.norm.pdf(flins, np.mean(f)*1e-9, np.std(f)*1e-9)
                            for f, flins in zip(frequency, frequency_linspaces)]
    normal_pairs = [[flins, ndist] for flins, ndist 
                    in zip(frequency_linspaces, normal_distributions)]
    normal_pairs[1].reverse()
    for i in range(2):
        axh.append(plt.subplot(grid_places[i]))        
        # Histogram
        n, b, p = axh[i].hist(frequency[i]*1e-9, nbins, density=True,
                              alpha=0.4, facecolor=colors[i],                               orientation=orientations[i])
        # Curve over histogram
        axh[i].plot(*normal_pairs[i], color=colors[i])    
        # Format
        axh[i].axis('off')
        function_lims[i](new_lims[i])
        limsh.append(function_lims_T[i]())
        # Mean values
        line_functions[i](np.mean(frequency[i])*1e-9, *limsh[i], 
                          colors=colors[i], linestyle='--')     
    del i

    ivs.saveFig(figsFilename('FvsF0Identity'), extension=figs_extension, 
                folder=figs_folder, overwrite=overwrite)

#%% *) FREQUENCY AND LENGTH FIT ON AIR (Eef) - DO NOT USE AGAIN LIGHTLY!

# Experimental values
physics['diameter'] = np.mean(width[0])
physics['area'] = np.pi * physics['diameter']**2 / 4

physics['K1_fsilica'] = physics['shear_fsilica'] * 2.75 # Pa
physics['K2_fsilica'] = np.pi * physics['diameter'] 
aux = np.sqrt(physics['density_fsilica'] * physics['shear_fsilica']) 
physics['K2_fsilica'] = physics['K2_fsilica'] * aux # Pa.s (viscosity's units)
del aux

# Theory model
def f_simple(length, young):
    f_0 = (np.sqrt(young/physics['density_gold']) / (2 * length))
    return f_0

if load_sem:
    young_gold_effective = iva.nonLinearFit(length[0], frequency[0], 
                                            f_simple, showplot=False)[-1][0]
    print(r"Módulo de Young: {}".format(ivu.errorValueLatex(
            *young_gold_effective, 
            units="Pa")))
    chi_squared = sum( (f_simple(length[0],
                                 young_gold_effective[0]) - frequency[0])**2 ) 
    chi_squared = chi_squared / len(length[0])

#%% *) FREQUENCY AND LENGTH FIT ON AIR (Eef PLOT) - DO NOT USE AGAIN LIGHTLY!
    
if load_sem:

    # Make a plot with fit and predictions
    plt.figure()
    ax = plt.subplot()
    plt.title('Análisis de resultados')
    plt.ylabel('Frecuencia (GHz)')
    plt.xlabel('Longitud (nm)')
    plt.plot(length[0]*1e9, frequency[0]*1e-9,'o', label=full_series[0])
    plt.plot(length[0]*1e9, 
             1e-9*f_simple(length[0], young_gold_effective[0]), 
             '-r', 
             label=r"Ajuste en vacío con {}".format(ivu.errorValueLatex(
                         *young_gold_effective, 
                         units="Pa")))
    ax.minorticks_on()
    ax.tick_params(axis='y')
    ax.tick_params(axis='y', which='minor', length=0)
    ax.grid(axis='both', which='both')
    ax.legend()
    plt.show()
    
    # Save plot
    ivs.saveFig(figsFilename('YoungEfectivo'), extension=figs_extension, 
                folder=figs_folder, overwrite=overwrite)

#%% FREQUENCY AND LENGTH MODEL ON TA2O5 (G) - DO NOT USE AGAIN LIGHTLY!

if load_sem and filter_not_in_common_rods:
    
    def G_iv_ta2o5(ta2o5_frequency, air_frequency, width):
        area = np.pi * width**2 / 4
        num = (ta2o5_frequency**2 - air_frequency**2) * 4 * np.pi**2
        den = 2.75 / (physics['density_gold'] * area)
        den = den - ( ( (np.pi * width) / (2 * physics['density_gold'] * 
                     area))**2 )* physics['density_ta2o5']
        return num/den
    
    G = {r'No Mean [NM] ($F_i$, $d_i$)' : G_iv_ta2o5(frequency[1], 
                                                     frequency[0], 
                                                     width[0]),
         r'All Mean [AM] ($\bar{F}$, $\bar{d}$)' : G_iv_ta2o5(np.mean(frequency[1]),
                                                              np.mean(frequency[0]),
                                                              np.mean(width))
        }    
    full_name_key = lambda k : k[0 : k.find('[')-1]
    short_name_key = lambda k : k[k.find('[')+1 : k.find(']')]
    symbols_key = lambda k : k[k.find('(') : k.find(')')+1]
    
    for k in G.keys():
        if 'AM' not in k:
            print("G = {} +- {} for {} method".format(
                *ivu.errorValue(np.mean(G[k]),np.std(G[k]),units='Pa'),
                short_name_key(k)))
        else:
            print("G = {:.1f} GPa for AM method".format(G[k]/1e9))

#%% FREQUENCY AND LENGTH MODEL ON TA2O5 (BETA TERM) - DO NOT USE AGAIN LIGHTLY!

if load_sem and filter_not_in_common_rods:

    # Let's study if the approximation was valid
    for k in G.keys():
        beta_term = np.pi**2 * physics['viscosity_gold'] / (
                2 * physics['length']**2 * physics['density_gold'])
        K2_term = np.pi * physics['length'] * np.sqrt(
                physics['density_ta2o5'] * np.mean(G[k])) / (
                2 * physics['density_gold'] * physics['area'])
        print("K2 term is {:.2f} times bigger than Beta term for {} method".format(
                K2_term/beta_term,
                short_name_key(k)))

#%% FREQUENCY AND LENGTH MODEL ON TA2O5 (PLOT) - DO NOT USE AGAIN LIGHTLY!

if load_sem and filter_not_in_common_rods:

    def f_simple(length, young):
        f_0 = (np.sqrt(young/physics['density_gold']) / (2 * length))
        return f_0
    
    def f_iv_ta2o5(frequency_air, G_ta2o5):
        K1_ta2o5 = G_ta2o5 * 2.75 # Pa
        K2_ta2o5 = np.pi * physics['diameter']
        aux = np.sqrt(physics['density_ta2o5'] * G_ta2o5) 
        K2_ta2o5 = K2_ta2o5 * aux # Pa.s (viscosity's units)
        K1_term = ( K1_ta2o5 / 
                  ( np.pi**2 * physics['density_gold'] * physics['area'] ) )
        K2_subterm = ( K2_ta2o5 / 
                     ( 2 * np.pi * physics['density_gold'] * physics['area'] ) )
        f = np.sqrt(frequency_air**2 + K1_term/4 - (K2_subterm)**2/4 )
        return f
    
    plt.figure()
    plt.plot(length[0]*1e9, frequency[1]/1e9, 'ok', label='Datos')
    plt.ylabel('Frecuencia f (GHz)')
    plt.xlabel('Longitud L (nm)')
    full_lines = []
    dotted_lines = []
    for k, c in zip(G.keys(), ['r', 'g', 'b']):
        plt.plot(length[0]*1e9, f_iv_ta2o5(f_simple(length[0],
                                                    young_gold_effective[0]),
                                           np.mean(G[k]))/1e9,
                 '-'+c,
                 label=r'Método $E_{ef}$ ' + symbols_key(k))[0]
        plt.plot(length[0]*1e9, f_iv_ta2o5(frequency[0],
                                           np.mean(G[k]))/1e9,
                 ':'+c,
                 label=r'Método $F$ ' + symbols_key(k))[0]
    plt.legend()
    
#%% *) FREQUENCY AND LENGTH - DO NOT USE AGAIN LIGHTLY!

# --> Final plots of F vs L

if load_sem and filter_not_in_common_rods:

    # Make a plot with fit and predictions
    plt.figure()
    ax = plt.subplot()
    plt.title('Análisis de resultados')
    plt.ylabel('Frecuencia (GHz)')
    plt.xlabel('Longitud (nm)')
    for j in range(len(series)):
        plt.plot(length[j]*1e9, frequency[j]*1e-9,'o', label=full_series[j])
    plt.plot(length[0]*1e9, 
             1e-9*f_simple(length[0], young_gold_effective[0]), 
             '-r', 
             label=r"Ajuste en vacío con {}".format(ivu.errorValueLatex(
                         *young_gold_effective, 
                         units="Pa")))
    plt.plot(length[0]*1e9, f_iv_ta2o5(f_simple(length[0],
                                                young_gold_effective[0]),
                                       np.mean(G[list(G.keys())[-1]]))/1e9,
             '-k',
             label=r'Método $E_{ef}$ ' + symbols_key(list(G.keys())[-1]))[0]
    ax.minorticks_on()
    ax.tick_params(axis='y')
    ax.tick_params(axis='y', which='minor', length=0)
    ax.grid(axis='both', which='both')
    ax.legend()
    plt.show()
    
    # Save plot
    ivs.saveFig(figsFilename('PlotSuperFinal'), extension=figs_extension, 
                folder=figs_folder, overwrite=overwrite)

"""
Multiple legends on the same Axes
https://matplotlib.org/users/legend_guide.html
"""