import iv_analysis_module as iva
import iv_plot_module as ivp
import iv_save_module as ivs
import os

#%% PARAMETERS --------------------------------------------------------------------

name = 'M_20190408_03.txt'
path = os.getcwd()
autoclose = True
autosave = True
round_Matlab_needed = True # Pyhon 3.6.2 needs it

#max_nrepetitions = 3 # BEWARE OF FALSE COLUMNS!

#%% LOAD AND CROP DATA ------------------------------------------------------------

# Load data
filename = os.path.join(path, name)
t, V, meanV, details = ivs.loadZeroPumpProbe(filename)
dt = details['dt']

#%% LINEAR PREDICTION -------------------------------------------------------------

# Use linear prediction
results, others = iva.linearPrediction(t, meanV, dt, 
                                      autoclose=autoclose,
                                      round_Matlab_needed=round_Matlab_needed)

#%% Plot linear prediction
ivp.linearPredictionPlot(filename, others, autosave=autosave)

#%% PMUSIC -----------------------------------------------------------------

## Start by defining general parameters
#PMN = nsize//4 # WHY /4? "cantidad de mediciones"
#PMT = 1200 # size of the window in ps
#PMdt = 20 # time step in ps
#
## Now get PMUSIC's parameters3
#PMdata = detrend(meanV) # BEWARE! THE BEST FIT HORIZONTAL LINE IS FILTERED!
#Mp = [PMN, 200] # This is PMUSIC's most important parameter
## Mp = [components' dimension, harmonics' limit]
## Second variable marks how many harmonics to throw away.
## It can't be greater than the measurement's dimension.
#
## Then define several variables to be filled
#MSn = []
#Mfn = []
#iPMindex=0
#for i in range(PMT+1, 1350, PMdt): # WHY? SHOULD I BE INCLUDING 1350? I'M NOT.
#
#    iPMindex = iPMindex + 1
#
#    # Take a segment of data and apply PMUSIC
#    iPMdata = PMdata[((i-PMT) < t) & (t < i)]
#    [MSn1, Mfn1] = pmusic(iPMdata, Mp, 6000, samplerate, [], 0)
#    # WHY 6000?
#        
#    iPMselection = ((Mfn1 >= 0) & (Mfn1 <= 0.06));
#    MSn[:, iPMindex] = MSn1[iPMselection]
#    Mfn[:, iPMindex] = Mfn1[iPMselection]
#
## Finally... Plot! :)
#plt.figure(1)
#plt.subplot(3,2,1)
#plt.imagesc(np.arange(1,T), Mfn[:,1], MSn)

"""
Problems so far:
    Don't have a pmusic Python equivalent
    Haven't searched an imagesc equivalent
"""

