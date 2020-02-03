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

name = 'All' # 'FS', 'Ta2O5'

# Analysis Parameters
load_sem = True
filter_not_in_common_rods = False
filter_air_outliers = False # two largest f
filter_Ta2O5_outliers = False # largest L
filter_FS_not_common = False # rod '9,10'
symbols=['% ', '$\pm$'] #['', '±']

# Saving parameters
home = r'C:\Users\Valeria\OneDrive\Labo 6 y 7'
figs_extension = '.png'
overwrite = True

# Plotting Parameters
show_air_outliers = False
fvsf0_nbins = 10
fvsf0_set_bigger_percent = 20

#%% OTHERS PARAMETERS ---------------------------------------------------------

# Data Parameters -for each data to compare, we need one value on each list-
desired_frequency = [12, 16, 8] # in GHz
full_series = [r'FS + Aire', 
               r"FS + Ta$_2$O$_5$",
               r'Ta$_2$O$_5$ + Aire']
series = ['LIGO1', 'LIGO1_PostUSA', 'LIGO5bis']
sem_series = ['LIGO1_1', 'LIGO1_1', 'LIGO5bis_1'] # ['nan']
sem_full_series = ['L1 1', 'L1 1', 'L5bis 1']

# Data selection
if name not in ['FS', 'Ta2O5', 'All']:
    raise ValueError("Accepted routine names: 'All', 'FS', 'Ta2O5'")
if name != 'All':
    if name == 'FS':
        desired_frequency = desired_frequency[:2]
        full_series = full_series[:2]
        series = series[:2]
        sem_series = sem_series[:2]
        sem_full_series = sem_full_series[:2]
    else:
        desired_frequency = desired_frequency[2]
        full_series = full_series[2]
        series = series[2]
        sem_series = sem_series[2]
        sem_full_series = sem_full_series[2]

# Parameters logic
if 'LIGO5bis' in series and filter_not_in_common_rods:
    filter_not_in_common_rods = False
    print("Watch it! Since loading L5bis, I'm not using only common rods")
if filter_not_in_common_rods:
    filter_FS_not_common = False
if filter_air_outliers:
    show_air_outliers = False

# Some more strings
filter_conditions = ''
if filter_not_in_common_rods or filter_FS_not_common:
    filter_conditions = 'CommonRods'
if filter_air_outliers and filter_Ta2O5_outliers:
    filter_conditions = 'NoOutliers'
elif filter_air_outliers:
    filter_conditions = 'NoAirOutliers'
elif filter_Ta2O5_outliers:
    filter_conditions = 'NoTa2O5Outliers'    

# Some functions and variables to manege filenames
def semFilename(sem_series, home=home):
    """Given a series 'M135_7B_1D', returns path to SEM data"""
    filename = 'Resultados_SEM_{}.txt'.format(sem_series)
    sem_series = sem_series.split('_') # From 'M_20190610_01' take '20190610'
    sem_filename = os.path.join(home, 'Muestras\SEM', *sem_series, filename)
    return sem_filename
rodsFilename = lambda series : os.path.join(
        home, r'Análisis\Rods_{}.txt'.format(series))
paramsFilename = lambda series : os.path.join(
        home, r'Análisis\Params_{}.txt'.format(series))
figsFilename = lambda fig_name : os.path.join(home, fig_name+'.png')
if filter_conditions != '':
    full_name = 'ComparedAnalysis_{}_{}'.format(name, filter_conditions)
else:
    full_name = 'ComparedAnalysis_{}'.format(name)
figs_folder = r'Análisis\{}'.format(full_name)

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
    
    ############ FILTERING ALGORITHMS #########################################
    
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
    
    if s == 'LIGO5bis' and filter_Ta2O5_outliers and load_sem:
        index = np.argsort(slength)[:-1] # Remove the largest length
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
    
    if (s=='LIGO1' or s=='LIGO1postUSA') and filter_FS_not_common:
        i = srods.index('9,10')
        index = list(range(len(srods)))
        index.pop(i) # Remove one rod
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
        
    ############ ORDER DATA AND SAVE TO OUT-LOOP VARIABLES ####################
        
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
    whole_filename = os.path.join(
            home, figs_folder,
            'Resultados_{}.txt'.format(full_name))
    whole_data = np.array([*sem_data[0][:,:8].T, 
                           *[fits_data[j][:,0] for j in range(len(full_series))],
                           *[fits_data[j][:,1] for j in range(len(full_series))], 
                           *[fits_data[j][:,2] for j in range(len(full_series))]
                           ])
    ivs.saveTxt(whole_filename, whole_data.T, 
                header=["Ancho (nm)", "Error (nm)",
                        "Longitud (nm)", "Error (nm)", 
                        "Relación de aspecto", "Error",
                        "Ángulo (º)", "Error (º)",
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
        whole_filename = os.path.join(
                home, figs_folder,
                'Resultados_{}_{}.txt'.format(full_name, s))
        whole_data = np.array([*sem_data[j][:,:8].T, 
                               *fits_data[j][:,:3].T])
        ivs.saveTxt(whole_filename, whole_data.T, 
                    header=["Ancho {} (nm)".format(ss), "Error (nm)",
                            "Longitud {} (nm)".format(ss), "Error (nm)", 
                            "Relación de aspecto {}".format(ss), "Error",
                            "Ángulo (º)", "Error (º)",
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
    ax1.set_xlabel('Nanoantena (Total {} NPs)'.format(len(rods[0])))
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
        ax1.set_xlabel('Nanoantena (Total {} NPs)'.format(len(rods[j])))
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
    
    del j, s, fs

#%% *) BOXPLOTS FQ and Ld

# Series to plot
make_boxplot_of = list(range(len(full_series)))
population = [fits_data[s].shape[0] for s in range(len(full_series))]
labels = [('{}\n'+r'$\hookrightarrow${} NPS').format(full_series[s], population[s]) 
          for s in make_boxplot_of]
del population

# Data inside each series to plot
choose_index_from_data_header = [[2,0], [0,2]]
boxplot_data_series = [sem_data, fits_data]
ax_labels = [[r"Longitud $L$ (nm)", r"Diámetro $d$ (nm)"],
              [r"Frecuencia $F$ (GHz)", r"Factor de calidad $Q$"]]
name_mask = ['BoxplotsLd{}', 'BoxplotsFQ{}']

for k, c in enumerate(choose_index_from_data_header):
    
    # Define boxplot data
    boxplot_data = [[boxplot_data_series[k][j][:,i] for j in make_boxplot_of] 
                    for i in c]
    
    # Format
    base_height = .1
    base_width = .6
    label_right_space = .1
    label_left_space = .08
    if len(make_boxplot_of)==1:
        alpha = [.25, .5]
    elif len(make_boxplot_of)==2:
        alpha = [0.12, .28]
    elif len(make_boxplot_of)==3:
        alpha = [0, .05]
    
    # Begin Figure
    fig = plt.figure()
    grid = plt.GridSpec(len(choose_index_from_data_header[k]), 1, hspace=0.1)
    ax = [plt.subplot(g) for g in grid]
    
    index = 0
    for a, dat, lab in zip(ax, boxplot_data, ax_labels[k]):
        
        # Boxplot
        bplot = a.boxplot(
            dat, 
            showmeans=True, meanline=True, 
            meanprops={'color':'k', 'linewidth':2, 'linestyle':':'},
            medianprops={'color':'r', 'linewidth':2},
            flierprops={'markersize':7},
            patch_artist=True,
            widths=base_width,
            labels=labels,
            vert=False)
        for p in bplot['boxes']:
            p.set_facecolor('w') # paint white boxes
        del p, bplot    
    
        # Labels' format
        a.xaxis.set_label_text(lab, va='center')
    
        # Grid's format
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.grid(which='major', axis='x')
        a.grid(which='minor', axis='x', linestyle=':')
        a.grid(which='major', axis='y')
        a.yaxis.tick_right()
        a.yaxis.set_label_position('right')
        
        # Axes size
        box = a.get_position()
        w = box.width
        box.x0 = box.x0 - w * label_left_space
        box.x1 = box.x1 - w * label_right_space
        box.y1 = box.y0 + base_height * len(make_boxplot_of)
        box.y0 = box.y0 + alpha[index]
        box.y1 = box.y1 + alpha[index]
        a.set_position(box)
        
        index += 1
    
    del index, box, w, a
    
    ax[0].xaxis.tick_top()
    ax[0].xaxis.set_label_position('top')
    
    ivs.saveFig(figsFilename(name_mask[k].format(make_boxplot_of)), 
                extension=figs_extension, 
                folder=figs_folder, overwrite=overwrite)
    
del make_boxplot_of, labels
del base_height, base_width, label_right_space, label_left_space
del choose_index_from_data_header, boxplot_data, ax_labels
del ax, grid

#%% PRINT AND COPY RESULTS

variables = ['Longitud L', 'Diámetro d', 'Relación de aspecto',
             'Ángulo', 'Frecuencia F', 'Factor de calidad Q',
             'Tiempo de decaimiento']
variables_units = ['nm', 'nm', '', 'º', 'GHz', '', 'ps']
variables_data = lambda i : [
        iva.getValueError(sem_data[i][:,2], sem_data[i][:,3]),
        iva.getValueError(sem_data[i][:,0], sem_data[i][:,1]),
        iva.getValueError(sem_data[i][:,4], sem_data[i][:,5]),
        iva.getValueError(sem_data[i][:,6], sem_data[i][:,7]),
        iva.getValueError(fits_data[i][:,0]),
        iva.getValueError(fits_data[i][:,2]),
        iva.getValueError(fits_data[i][:,1]),]

values_string = '{}Serie {} {}\n\n'.format(symbols[0], name, filter_conditions)
for i, fs in enumerate(full_series):
    values_string += "{}Resultados de {}\n".format(symbols[0], fs)
    values_string += "{}Cantidad de NPs: {:.0f}\n".format(
            symbols[0],
            len(fits_data[i][:,0]))
    for v, vd, vu in zip(variables, variables_data(i), variables_units):
        values_string += '{}{} = {}\n'.format(symbols[0], 
                                              v, 
                                              ivu.errorValueLatex(
                                                      *vd, 
                                                      units=vu, 
                                                      symbol=symbols[1]))
    values_string += '\n'

print(values_string)
ivu.copy(values_string)

variables = ['Amplitud (mVpp)', r'Potencia ($\mu$W)', 'Longitud de onda (nm)',
             'Ancho medio de la campana (nm)']
variables_units = ['mVpp', r'$\mu$W', 'nm', 'nm']
variables_data = lambda i : [
        iva.getValueError(params[i][:,0]),
        iva.getValueError(params[i][:,1]),
        iva.getValueError(params[i][:,2]),
        iva.getValueError(params[i][:,3])]

values_string = '{}Serie {} {}\n\n'.format(symbols[0], name, filter_conditions)
for i, fs in enumerate(full_series):
    values_string += "{}Parámetros de {}\n".format(symbols[0], fs)
    values_string += "{}Cantidad de NPs: {:.0f}\n".format(
            symbols[0],
            len(fits_data[i][:,0]))
    for v, vd, vu in zip(variables, variables_data(i), variables_units):
        values_string += '{}{} = {}\n'.format(symbols[0], 
                                              v, 
                                              ivu.errorValueLatex(
                                                      *vd, 
                                                      units=vu, 
                                                      symbol=symbols[1]))
    values_string += '\n'

print(values_string)
ivu.copy(values_string)

#%% *) FREQUENCY ON SAME RODS - ANDREA'S

if filter_not_in_common_rods:
    
    plot_data = [[f*1e-9 for f in frequency[:2]],
                 [(f*1e-9)**2 for f in frequency[:2]]]
    plot_symbols = [[r'f', r'f_0'],
                    [r'f^2', r'f_0^2']]
    plot_name = ['FvsF0',
                 'FvsF0Sqrd']
    
    for pd, ps, pn in zip(plot_data, plot_symbols, plot_name):
    
        fig = plt.figure(figsize=[4.8, 4.8]) # size in inches
        grid = plt.GridSpec(5, 5, wspace=0, hspace=0)
        
        ax = plt.subplot(grid[1:,:-1])
        ax.plot(pd[0], pd[1], 'ko', markersize=8, mfc='w')
        axis_functions = [plt.xlabel, plt.ylabel]
        for f, s, fs in zip(axis_functions, ps, full_series):
            f(r"Frecuencia ${}$ (GHz) $\rightarrow$ {}".format(s, fs))
        
        # Grid's format
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(which='major', axis='both')
        ax.grid(which='minor', axis='both', linestyle=':')
        
        # Limits' format
        delta = max(pd[1]) - min(pd[0])
        lims = [(min(pd[0]) - fvsf0_set_bigger_percent*delta/100), 
                (max(pd[1]) + fvsf0_set_bigger_percent*delta/100)]
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        new_lims = [lims, lims]
        new_lims_T = [lims, lims]
    
        il = []
        delta_data = pd[1] - pd[0] # Hz
        
        # Identity
        data_linspaces = [np.linspace(nl[0], nl[1], 50) for nl in new_lims]
        il.append(ax.plot(
                 data_linspaces[0], data_linspaces[0], '-k',
                 label=r'${b} = {a}$'.format(a=ps[0], b=ps[1]))[0])
        
        # Identity with mean difference vertical shift
        il.append(ax.plot(
            data_linspaces[0], 
            data_linspaces[0] + np.mean(delta_data),
            '--k',
            label=r'${b} = {a} + \langle {b} - {a} \rangle$'.format(
                    a=ps[0], b=ps[1]))[0])
        
        # Identity with mean difference standard deviation vertical shift
        ax.fill_between(
            data_linspaces[0], 
            data_linspaces[0] + (np.mean(delta_data) 
                    - np.std(delta_data)),
            data_linspaces[0] + (np.mean(delta_data) 
                    + np.std(delta_data)),
            color='m',
            alpha=0.2)
        
        # Mean values
        line_functions = [plt.vlines, plt.hlines]
        colors = ['blue', 'red']
        ml = []
        for i in range(2):
            ml.append(line_functions[i](
                     np.mean(pd[i]), 
                     *new_lims_T[i], colors=colors[i], linestyle='--',
                     label=r'$\langle {} \rangle$'.format(ps[i])))
        del i
        
        # Standard deviation
        fill_function = [ax.fill_betweenx, ax.fill_between]
        for i in range(2):
            fill_function[i](new_lims_T[i], 
                             (np.mean(pd[i])-np.std(pd[i])),
                             (np.mean(pd[i])+np.std(pd[i])),
                             color=colors[i],
                             alpha=0.1)
        
        # Histograms
        axh = []
        limsh = []
        grid_places = [grid[0,:-1], grid[1:,-1]]
        orientations = ['vertical', 'horizontal']
        function_lims = [plt.xlim, plt.ylim]
        function_lims_T = [plt.ylim, plt.xlim]
        normal_distributions = [st.norm.pdf(dlins, np.mean(d), np.std(d))
                                for d, dlins in zip(pd, data_linspaces)]
        normal_pairs = [[flins, ndist] for flins, ndist 
                        in zip(data_linspaces, normal_distributions)]
        normal_pairs[1].reverse()
        for i in range(2):
            axh.append(plt.subplot(grid_places[i]))        
            # Histogram
            n, b, p = axh[i].hist(pd[i], fvsf0_nbins, density=True,
                                  alpha=0.4, facecolor=colors[i], 
                                  orientation=orientations[i])
            # Curve over histogram
            axh[i].plot(*normal_pairs[i], color=colors[i])    
            # Format
            axh[i].axis('off')
            function_lims[i](new_lims[i])
            limsh.append(function_lims_T[i]())
            # Mean values
            line_functions[i](np.mean(pd[i]), *limsh[i], 
                              colors=colors[i], linestyle='--')     
        del i
    
        leg = plt.legend(handles=il, loc=(-.8, 1.05), frameon=False)
        ax_leg = plt.gca().add_artist(leg)
        ax.legend(handles=ml, loc=(1.03, 0), frameon=False)
    
        ivs.saveFig(figsFilename(pn), extension=figs_extension, 
                    folder=figs_folder, overwrite=overwrite)
        
        print(r"Varianza $\sigma^2$:")
        for s, d in zip(ps, pd):
            print(r"- ${}$ --> {}".format(s, np.var(d)))
        print(r"- ${} - {}$ --> {}".format(*ps, np.var(delta_data)))
    
    del pd, ps, pn