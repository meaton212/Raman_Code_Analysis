# load packages
import pandas as pd
import glob, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit

from tqdm import tqdm

from scipy.signal import find_peaks

plt.close('all')

#%% Define functions used in this code
def findpeaks(arr, h, w=1, d=1):
    # Adjust baseline to emulate 'DoubleSided'
    adjusted = arr  - arr.mean()
   # adjusted = arr
    indexes, _ = find_peaks(adjusted, height=h)
    # Retrieve extrema heights from original signal rather
    # than from the properties dict with adjusted heights
    heights = arr[indexes]

    return heights, indexes

#Converts nested list into single list
def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

#Rounds Number to Specified number of significant digits
def round_sig(x, sig=2):
        if x!=0:
            rnd=round(x, sig-int(np.floor(np.log10(abs(x))))-1)
        else:
            rnd=0
        return rnd

#%%Input data:
#The program expects multiple (up to 10) files in .txt or .dat format. In each file, the Raman shift
#should be included in the first column, followed by the Intensity data in
#consecutive columns (as many spectra available, and not neccesarily the same number of spectra per file).

#The location and names of the files should be indicated in path and file_name
#variables, respectively. The total number of files should be included in
#total.
path=r"C:\Users\matte\Documents\IMDEA\Data\Raman\2023\02-Feb\07-Test Lorentz\\"
os.chdir(path) #Changes folder to the selected path


file_name=['6,5-SWCNTs(Porf)_10mW_25x1000,9-15cm-1,2co,1s_2',
           'CuMINT_10mW_25x1000,9-15cm-1,2co,1s_0',
           ]

total=len(file_name) #Total Number of files


delim='\t'  #Define delimeter of data in data file
type_='.txt'
name=['6,5 SWNT','CuMINT']#Indicate here the names of the data for legend and titles

imgtype='.svg' #What type of images do you want to save


#Histogram Normalization?
#Set to True if you want your histograms to be normalized
#Else, set to False
den=True


#%%Normalization: 
#Choose normalization spectral rng. The spectra will be normalised to
#the maximun within the specified spectral rng: 
normLow=1500; #Lower limit in cm-1
normHigh=1700; #Upper limit in cm-1

#%%Peak identification: Intensity and Shifts
#Intensity ratio Id/Ig will be calculated by taking: max Int value-min Int
#value within the spectral rng selected. Select appropriate rng where 
#the full peaks are resolved. 
#Raman shifts for D and G modes are caculated as the position of the
#maximun intensity 
#RBM modes

#G band
band1Low=1450; #Lower limit in cm-1
band1High=1700; #Upper limit in cm-1
#D band
band2Low=1200; #Lower limit in cm-1
band2High=1400; #Upper limit in cm-1
#2D band
band3Low=2450; #Lower limit in cm-1
band3High=2800; #Upper limit in cm-1


#RBM modes
rbm=1# Set to 1 if RBM analysis desired
RBMregion_Low=200; #Lower limit in cm-1
RBMregion_High=360; #Upper limit in cm-1
Prom=[0.01]; #This value sets the max limit at which peaks will be considered. 
    #This can be a single number, or a list of values that different for each file
    #Local maxima with a maximum height lower than Prom will be discarded.
if len(Prom)==1 and total!=1:
    Prom=Prom*np.ones(total)
if rbm==1:
    meanRBM=dict.fromkeys(name)
    stdRBM=dict.fromkeys(name)    
    
#%% Model: Choose the method to analyse spectral features (simple peak
# identification or lorentzianpeak fitting)

lorentz=1; # lorentz fits for G, D and 2D (=0 for simple identification; =1 for lorentzian peak fitting)
nt=1; # Splitting G band in G+ and G-?

# If lorentzian fits are selected, the user is required to indicate
# initial guess for the peak position, FWHM and intensity

Gmin_init=[[1530, 20, 0.4], [1530, 25, 0.31]] #[center FWHM intensity_max] for G-, only if nt=1
Gplus_init=[[1580, 30, 1], [1580, 25, 1]] #[center FWHM intensity_max] fof G+
D_init=[[1300, 40, 0.08], [1325, 10, 0.04]] #[center FWHM intensity_max] for D band
init_2D=[[2600, 50, 0.4], [2640, 50, 0.3]] #[center FWHM intensity_max] for 2D band



#%%Ploting options:
#Choose the desired output plots (1 is yes and 0 is no). 
raw=1; #Make figure af all raw spectra.
norm=1; #Plot all normalised spectra.
rng=1; #Plot spectral regions chosen for G, D and 2D bands in intensity calculation
peaks=1; #plot found peaks in the RBM region, Note that if the
#number of spectra is very high, the computing is going to slow down
#significantly.

correlations=1 #Plot correlations between peaks shifts and Id/Ig


#%%Fontsize for Plots
fs=16 # This is the fontsize you want to use for all the plots
width=3 #width of histograms




#%%Is it a 2D map? In case it is and you want the 3D plots in form of maps
maps=0; #set to 0 if not a maps and to 1 if you want a 2D Map
map_var='I'  # What do you want to map?
                                         #'I' for Intensity Ratio, 
                                         #'G' for Gband
                                         #'D' for Dband
                                         #'2D' for 2DBand
rows=32 ; #Number of rows in 2D Map
col=32 ; #Number of columns in 2D Map
use_leng=1; #If 1, then map will use x and y dimensions defined below. Else if 0, will use pixel number
xlen=42; #x len in um
ylen=38;  #y len in um


#%%Code Starts here:
for z in range(0,total):

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#Import data #%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%

    data = pd.read_csv(path+file_name[z]+type_,delimiter=delim,decimal='.',header=None);
    data=data.to_numpy()
    Shift = np.array(data[:,0]); #Raman shift is the first column of the file
    Intensity=data[:,1:]; #Raman intensities are all data from column 2 til end
    
    rng2=Shift>0
    Intensity=Intensity[rng2,:]
    Shift=Shift[rng2]
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#   
    #%%#%%#%%#Set color scales for plots #%%#%%#%%#%%
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
    gray=cm.get_cmap('gist_yarg',len(Intensity[0,:])+20)
    winter=cm.get_cmap('winter',len(Intensity[0,:])+5)
    summer=cm.get_cmap('summer',len(Intensity[0,:])+5)
    autumn=cm.get_cmap('autumn',len(Intensity[0,:])+5)
    bone= cm.get_cmap('bone',len(Intensity[0,:])+5)
    cool=cm.get_cmap('cool',len(Intensity[0,:])+5)
    spring=cm.get_cmap('spring',len(Intensity[0,:])+5)
    copper=cm.get_cmap('copper',len(Intensity[0,:])+5)
    pink=cm.get_cmap('pink',len(Intensity[0,:])+5)
    hot=cm.get_cmap('hot',len(Intensity[0,:])+5)
   
    cmap=[gray,winter,summer,autumn,
          cool,spring,copper,pink,hot,bone]
    #Sets color maps for different files
    
    color=cmap[z];
    clr=color(np.linspace(0,1,len(Intensity.transpose())))
    
    if total==1:
        div=[1, 1];
    elif total==2:
        div=[1, 2];
    elif total==3:
        div=[1, 3];
    elif total==4:
        div=[2, 2];
    elif total==5:
        div=[2, 3];
    elif total==6:
        div=[2, 3];
    elif total==10:
        div=[2, 5];
    else:
        div=[3, 3];

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%Normalization #%%#%%#%%#%%#%%#%%#%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%

    indLow=np.where(abs(Shift-normLow)==min(abs(Shift-normLow)))[0][0]
    indHigh=np.where(abs(Shift-normHigh)==min(abs(Shift-normHigh)))[0][0];
    
    Intensity_norm=np.empty(np.shape(Intensity));
    for n in range(0,len(Intensity[1,:])):
        Intensity_norm[:,n]=(Intensity[:,n]-min(Intensity[:,n])); #substract min
        Intensity_norm[:,n]=Intensity_norm[:,n]/max(Intensity_norm[indLow:indHigh,n]); #divide by norm
    
    
    Intensity_av=np.mean(Intensity_norm,axis=1); #Calculate average spectra
    
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
    #%%#%%#%%#%%G and D modes: Intensity and Shift #%%#%%#%%#%%#%%
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
    
    
    ind1Low=np.where(abs(Shift-band1Low)==min(abs(Shift-band1Low)))[0][0]
    ind1High=np.where(abs(Shift-band1High)==min(abs(Shift-band1High)))[0][0];
    
    ind2Low=np.where(abs(Shift-band2Low)==min(abs(Shift-band2Low)))[0][0]
    ind2High=np.where(abs(Shift-band2High)==min(abs(Shift-band2High)))[0][0];
    
    ind3Low=np.where(abs(Shift-band3Low)==min(abs(Shift-band3Low)))[0][0]
    ind3High=np.where(abs(Shift-band3High)==min(abs(Shift-band3High)))[0][0];
    
    Intensity_1rng=Intensity_norm[ind1Low:ind1High,:];
    Intensity_2rng=Intensity_norm[ind2Low:ind2High,:];
    Intensity_3rng=Intensity_norm[ind3Low:ind3High,:];
    
    Shift_rng1=Shift[ind1Low:ind1High];
    Shift_rng2=Shift[ind2Low:ind2High];
    Shift_rng3=Shift[ind3Low:ind3High];
    
    Int_G=np.array([]);
    Int_D=np.array([]);
    Int_2D=np.array([]);
    
    center_G=np.array([]);
    center_D=np.array([]);
    center_2D=np.array([]);
    
    if nt==1:
        I_G_fit=[]
        center_Gmin=[]
        FWHM_Gmin=[]
        Int_Gmin=[]
        center_Gplus=[]
        FWHM_Gplus=[]
        Int_Gplus=[]
    if lorentz==1:
        I_D_fit=[]
        FWHM_D=[]
        I_2D_fit=[]
        FWHM_2D=[]
        if nt==0:
            I_G_fit=[]
            FWHM_G=[]
            
        
        
    for n in tqdm(range(0,len(Intensity_norm[0,:]))):
        
       if np.logical_and(lorentz==0,nt==0)==1:
           #peak 1/G: 
           Int_G=np.append(Int_G,max(Intensity_1rng[:,n])-min(Intensity_1rng[:,n]));
           indMax1=np.where(Intensity_1rng[:,n]==max(Intensity_1rng[:,n]))[0][0];
           center_G=np.append(center_G,Shift_rng1[indMax1]);
       
       if lorentz==0:
           #peak 2/D: 
           Int_D=np.append(Int_D,max(Intensity_2rng[:,n])-min(Intensity_2rng[:,n]));   
           indMax2=np.where(Intensity_2rng[:,n]==max(Intensity_2rng[:,n]))[0][0];
           center_D=np.append(center_D,Shift_rng2[indMax2]);
           
           #peak 3/2D: 
           Int_2D=np.append(Int_2D,max(Intensity_3rng[:,n])-min(Intensity_3rng[:,n]));
           indMax3=np.where(Intensity_3rng[:,n]==max(Intensity_3rng[:,n]))[0][0];
           center_2D=np.append(center_2D,Shift_rng3[indMax3]);
    
    
       
       #%%G+ and G- fitting   
       if nt==1:
           
           #Define initial guesses
           InitGuess_G=[Gmin_init[z][2]*Gmin_init[z][1],
                        Gmin_init[z][0],
                        (Gmin_init[z][1]/2)**2,
                        Gplus_init[z][2]*Gplus_init[z][1],
                        Gplus_init[z][0],
                        (Gplus_init[z][1]/2)**2,
                        0.1]
           def fit_func(x, g1,g2,g3,g4,g5,g6,g7):
               return (g1/((x-g2)**2+g3)+g4/((x-g5)**2+g6)+g7)# lorentz

           gamma_G ,pcov2= curve_fit(fit_func, Shift_rng1, Intensity_1rng[...,n], InitGuess_G, 
                                     bounds=(0,np.inf),maxfev=10000)
                               
           I_G_fit.append(fit_func(Shift_rng1, *gamma_G))
           
           if gamma_G[1]<gamma_G[4]:
                  center_Gmin.append(gamma_G[1])
                  FWHM_Gmin.append(2*np.sqrt(gamma_G[2]))
                  Int_Gmin.append(gamma_G[0]/gamma_G[2])
                  center_Gplus.append(gamma_G[4])
                  FWHM_Gplus.append(2*np.sqrt(gamma_G[5]))
                  Int_Gplus.append(gamma_G[3]/gamma_G[5])
           else:
                 center_Gplus.append(gamma_G[1])
                 FWHM_Gplus.append(2*np.sqrt(gamma_G[2]))
                 Int_Gplus.append(gamma_G[0]/gamma_G[2])
                 center_Gmin.append(gamma_G[4])
                 FWHM_Gmin.append(2*np.sqrt(gamma_G[5]))
                 Int_Gmin.append(gamma_G[3]/gamma_G[5])
                 
        #%%Lorentizan Fits                
       if lorentz==1:
           
           #D bands
           #Define Initial guesses
           InitGuess_D=[D_init[z][2]*D_init[z][1],
                        D_init[z][0],
                        (D_init[z][1]/2)**2,
                        0];
           
           def fit_func2(x, g1,g2,g3,g4):
               return (g1/((x-g2)**2+g3)+g4)# lorentz
           
           gamma_D ,pcov2= curve_fit(fit_func2, Shift_rng2, Intensity_2rng[...,n], 
                                      InitGuess_D,bounds=(0,np.inf),maxfev=10000)
           
           I_D_fit.append(fit_func2(Shift_rng2, *gamma_D))
           center_D=np.append(center_D,gamma_D[1]);
           FWHM_D.append(2*np.sqrt(gamma_D[2]));
           Int_D=np.append(Int_D,gamma_D[0]/gamma_D[2]); 
 
    
           #2D bands
           #Define Initial guesses
           InitGuess_2D=[init_2D[z][2]*init_2D[z][1],
                        init_2D[z][0],
                        (init_2D[z][1]/2)**2,
                        0];    
      
           
           gamma_2D ,pcov2= curve_fit(fit_func2, Shift_rng3, Intensity_3rng[...,n], 
                                      InitGuess_2D, bounds=(0,np.inf),maxfev=10000)
           
           I_2D_fit.append(fit_func2(Shift_rng3, *gamma_2D))
           center_2D=np.append(center_2D,gamma_2D[1]);
           FWHM_2D.append(2*np.sqrt(gamma_2D[2]));
           Int_2D=np.append(Int_2D,gamma_2D[0]/gamma_2D[2]);  
    
           if nt==0:
               #Lorentzian G fit if no G+/G- fit
               InitGuess_G=[Gplus_init[z][2]*Gplus_init[z][1],
                            Gplus_init[z][0],
                            (Gplus_init[z][1]/2)**2,
                            0];
              
               
               gamma_G ,pcov2= curve_fit(fit_func2, Shift_rng1, Intensity_1rng[...,n], 
                                InitGuess_G, bounds=(0,np.inf),maxfev=10000)
               
               I_G_fit.append(fit_func2(Shift_rng1, *gamma_G))
               center_G=np.append(center_G,gamma_G[1]);
               FWHM_G.append(2*np.sqrt(gamma_G[2]));
               Int_G=np.append(Int_G,gamma_G[0]/gamma_G[2]); 
               
 
    if nt==1:
        I21=Int_D/Int_Gplus;
        G_av=np.mean(center_Gplus)
        G_std=np.std(center_Gplus)    
        Gmin_av=np.mean(center_Gmin)
        Gmin_std=np.std(center_Gmin)  
    else:
        I21=Int_D/Int_G;
        G_av=np.mean(center_G)
        G_std=np.std(center_G)
        
    I21_av=np.mean(I21);
    I21_error=np.std(I21);
    D_av=np.mean(center_D)
    D_std=np.std(center_D)
    D2_av=np.mean(center_2D)
    D2_std=np.std(center_2D)
    

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
#%%#%%#%%#%%#RBM modes Shifts  #%%#%%#%%#%%#%%#%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
    
    if rbm == 1:
        indRBMLow=np.where(Shift-RBMregion_Low==min(abs(Shift-RBMregion_Low)))[0][0];
        indRBMHigh=np.where(Shift-RBMregion_High==min(abs(Shift-RBMregion_High)))[0][0];
        PeaksInt=[];
        PeaksLoc1=[];
    
    
        if peaks == 1:
            if z==0:
                fig_peaks,ax_peaks= plt.subplots(div[0],div[1],figsize=(div[1]*5,div[0]*5), constrained_layout=True)  #Create the figure only once, on the first loop
                fig_peaks.canvas.manager.set_window_title('Region Peaks')
                if total==1:
                    ax_peaks=[ax_peaks]
                else:
                    ax_peaks=ax_peaks.flatten()
                    
        for n in range(0,len(Intensity_norm[0,:])):    
            [pksRBM,locsRBM]=findpeaks(Intensity_norm[indRBMLow:indRBMHigh,n],Prom[z]);
            PeaksInt.append(pksRBM);
            PeaksLoc1.append(Shift[locsRBM+indRBMLow]);        
            
            if peaks == 1:
                ax_peaks[z].plot(Shift[indRBMLow:indRBMHigh],Intensity_norm[indRBMLow:indRBMHigh,n],color='C'+str(n))
                ax_peaks[z].plot(Shift[locsRBM+indRBMLow],pksRBM,ls='',marker='v',markersize=10,color='C'+str(n))
    
    
                ax_peaks[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
                ax_peaks[z].set_ylabel('Intensity / a.u',fontsize=fs)
                ax_peaks[z].tick_params(axis="both", labelsize=fs)   
                ax_peaks[z].set_title(name[z]+' :Region Peaks');
                ax_peaks[z].set_xlim((Shift[indRBMLow],Shift[indRBMHigh]))
    



#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
#%%#%%#%%Figures #%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%

    #Raw spectra
    if raw==1:
       if z==0:
           fig_raw,ax_raw= plt.subplots(div[0],div[1],figsize=(div[1]*5,div[0]*5), constrained_layout=True)  #Create the figure only once, on the first loop
           fig_raw.canvas.manager.set_window_title('Raw Data')
           if total==1:
               ax_raw=[ax_raw]
           else:
               ax_raw=ax_raw.flatten()

           
    
       for n,clr2 in enumerate(clr):
           ax_raw[z].plot(Shift,Intensity[:,n],color=clr2,label=str(n));
    

     #  ax_raw[z].legend()
       ax_raw[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
       ax_raw[z].set_ylabel('Intensity / a.u',fontsize=fs)
       ax_raw[z].tick_params(axis="both", labelsize=fs)   
       ax_raw[z].set_title(name[z]+': Raw data');
       ax_raw[z].set_xlim((0,np.max(Shift)))
#%%Normalised spectra
    if norm==1:
        if z==0:
           fig_norm,ax_norm= plt.subplots(div[0],div[1], figsize=(div[1]*5,div[0]*5), constrained_layout=True)
           fig_norm.canvas.manager.set_window_title('Normalized Spectra')
           if total==1:
               ax_norm=[ax_norm]
           else:
               ax_norm=ax_norm.flatten()    
        for n,clr2 in enumerate(clr):
           ax_norm[z].plot(Shift,Intensity_norm[:,n],color=clr2,label=str(n));
        #ax_norm[z].legend()  
        ax_norm[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_norm[z].set_ylabel('Intensity / a.u',fontsize=fs)
        ax_norm[z].tick_params(axis="both", labelsize=fs)   
        ax_norm[z].set_title(name[z]+': Normalized spectra');    
        ax_norm[z].set_ylim((0,np.max(Intensity_norm)))
        ax_norm[z].set_xlim((0,np.max(Shift)))
           

#%%Average spectra
    if z==0:
       fig_avg,ax_avg= plt.subplots(1,1, figsize=(6,6), constrained_layout=True)
       fig_avg.canvas.manager.set_window_title('Average Spectra')
           
    
    ax_avg.plot(Shift,Intensity_av,color=clr[round(len(clr)/2)],label=name[z]);
    ax_avg.legend(fontsize=fs)  
    ax_avg.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
    ax_avg.set_ylabel('Intensity / a.u',fontsize=fs)
    ax_avg.tick_params(axis="both", labelsize=fs)   
    ax_avg.set_title('Average spectra');



#%% #Spectral rng for Intensity ratio
    if rng ==1:
        if z==0:
           fig_rng,ax_rng= plt.subplots(div[0],div[1], figsize=(div[1]*5,div[0]*5), constrained_layout=True)
           fig_rng.canvas.manager.set_window_title('Range Used for Intensity Ratio Calculation')
           if total==1:
               ax_rng=[ax_rng]
           else:
               ax_rng=ax_rng.flatten()
        for n,clr2 in enumerate(clr):
            ax_rng[z].plot(Shift[ind1Low:ind1High],Intensity_1rng[:,n],color=clr2);
            ax_rng[z].plot(Shift[ind2Low:ind2High],Intensity_2rng[:,n],color=clr2);
            ax_rng[z].plot(Shift[ind3Low:ind3High],Intensity_3rng[:,n],color=clr2);
            
            
            
        ax_rng[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_rng[z].set_ylabel('Intensity / a.u',fontsize=fs)
        ax_rng[z].tick_params(axis="both", labelsize=fs)   
        ax_rng[z].set_title(name[z]+': Range for $I_{d}/I_{g}$ Calculation');    
        ax_rng[z].set_ylim((0,1))


#%% Histogram intensity ratio
    

    if z==0:
        fig_IR,ax_IR= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
        fig_IR.canvas.manager.set_window_title('Intensity Ratio')


    #Define Bin Size
    binsInt=np.arange(min(I21), max(I21) + width/250, width/250)
    if len(binsInt)==1:
        binsInt=np.arange(min(I21)- width/250, max(I21) + width/250, width/250)

         
    ax_IR.hist(I21,binsInt,color=clr[round(len(clr)/2)],
               label=name[z]+': $I_{d}/I_{g}$='+str(round(I21_av,4))+'$\pm$'
               +str(round(I21_error,4)),alpha=0.5,ec='k',align='left',density=den)

    ax_IR.legend(fontsize=fs-2)  
    ax_IR.set_xlabel('$I_{d}/I_{g}$',fontsize=fs)
    ax_IR.set_ylabel('counts',fontsize=fs)
    ax_IR.tick_params(axis="both", labelsize=fs)   
    ax_IR.set_title('Intensity ratio: $I_{d}/I_{g}$',fontsize=fs+2);


#%% Raman shift G mode
    
    if z==0:
        fig_G,ax_G= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
        fig_G.canvas.manager.set_window_title('Raman shift G band')
        if nt==1:
            fig_G_I,ax_G_I= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_G_I.canvas.manager.set_window_title('Intensity ratio G modes')

    if nt==0:
        

        binsShift=np.arange(min(center_G)-width, max(center_G) + width, width)
        
        
        # if len(binsShift)==1:
        #     ax_G.hist(center_G,width, color=clr[round(len(clr)/2)],
        #               label=name[z]+': $G$ band='+str(round(G_av,2))+'$\pm$'
        #               +str(round(G_std,2))+' $cm^{-1}$',
        #               alpha=0.5,ec='k',align='left')
            
        #else:
        ax_G.hist(center_G,binsShift,color=clr[round(len(clr)/2)],
                      label=name[z]+': $G$ band='+str(round(G_av,2))+'$\pm$'
                      +str(round(G_std,2))+' $cm^{-1}$',
                      alpha=0.5,ec='k',align='left',density=den)
        
        
 
        ax_G.set_title('Raman shift $G$ band',fontsize=fs+2);
    else:
      #  binsShift=round(max(center_Gplus)-min(center_Gmin));
        binsShift=np.arange(min(center_Gmin), max(center_Gplus) + width, width)
        ax_G.hist(center_Gplus,binsShift, color=clr[round(len(clr)/3)],
                  label=name[z]+': $G^{+}$ band='+str(round(G_av,2))+'$\pm$'
                  +str(round(G_std,2))+' $cm^{-1}$',
                  alpha=0.5,ec='k',align='left',density=den)
        ax_G.hist(center_Gmin,binsShift,color=clr[round(len(clr)/2)],
                  label=name[z]+': $G^{-}$ band='+str(round(Gmin_av,2))+'$\pm$'
                  +str(round(Gmin_std,2))+' $cm^{-1}$',
                  alpha=0.5,ec='k',align='left',density=den)
        ax_G.set_title('Raman shift $G$ modes',fontsize=fs+2);
        
        

        Iplus_minus=np.array(Int_Gplus)/np.array(Int_Gmin);
        binsInt=np.arange(min(Iplus_minus), max(Iplus_minus) + width/10, width/10)
        ax_G_I.hist(Iplus_minus,binsInt,color=clr[round(len(clr)/2)],
                  label=name[z]+': $I_{G^+}/I_{G^-}$ ='+str(round(np.mean(Iplus_minus),2))+'$\pm$'
                  +str(round(np.std(Iplus_minus),2)),
                  alpha=0.5,ec='k',align='left',density=den)
        ax_G_I.legend(fontsize=fs-2)  
        ax_G_I.set_xlabel('$I_{G^+}/I_{G^-}$',fontsize=fs)
        ax_G_I.set_ylabel('counts',fontsize=fs)
        ax_G_I.tick_params(axis="both", labelsize=fs)
        ax_G_I.set_title('Intensity ratio G modes',fontsize=fs+2);



    ax_G.legend(fontsize=fs-2)  
    ax_G.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
    ax_G.set_ylabel('counts',fontsize=fs)
    ax_G.tick_params(axis="both", labelsize=fs)  
#%% #Raman shift D mode

    if z==0:
        fig_D,ax_D= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
        fig_D.canvas.manager.set_window_title('Raman shift D band')



    binsShift=np.arange(min(center_D)-width, max(center_D) + width, width)

    ax_D.hist(center_D,binsShift,color=clr[round(len(clr)/2)],
              label=name[z]+': $D$ band='+str(round(D_av,2))+'$\pm$'
              +str(round(D_std,2))+' $cm^{-1}$',
              alpha=0.5,ec='k',align='left',density=den)
    
    
    ax_D.legend(fontsize=fs-2)  
    ax_D.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
    ax_D.set_ylabel('counts',fontsize=fs)
    ax_D.tick_params(axis="both", labelsize=fs)   
    ax_D.set_title('Raman shift $D$ band',fontsize=fs+2);
    ax_D.xaxis.set_major_locator(MaxNLocator(integer=True))
#%% #Raman shift 2D mode
    if z==0:
        fig_2D,ax_2D= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
        fig_2D.canvas.manager.set_window_title('Raman shift D 2band')



    binsShift=np.arange(min(center_2D)-width, max(center_2D) + width, width)
    
    ax_2D.hist(center_2D,binsShift,color=clr[round(len(clr)/2)],
               label=name[z]+': $2D$ band='+str(round(D2_av,2))+'$\pm$'
               +str(round(D2_std,2))+' $cm^{-1}$',
               alpha=0.5,ec='k',align='left',density=den)
    
    
    ax_2D.legend(fontsize=fs-2)  
    ax_2D.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
    ax_2D.set_ylabel('counts',fontsize=fs)
    ax_2D.tick_params(axis="both", labelsize=fs)   
    ax_2D.set_title('Raman shift $2D$ band',fontsize=fs+2);
    ax_2D.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    
#%%%

# Lorentz fit results

    if nt==1 or lorentz==1:
        if nt==1 and lorentz==0:
                
                fig_lor_fit,ax_lor_fit= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_lor_fit.canvas.manager.set_window_title('Lorenztian Fitting results')
                ax_lor_fit=[ax_lor_fit]
                for n,clr2 in enumerate(clr):
                    ax_lor_fit[0].plot(Shift_rng1,Intensity_1rng[:,n],color=clr2)
                    ax_lor_fit[0].plot(Shift_rng1,np.transpose(I_G_fit)[:,n],'--r');

        elif lorentz==1:
                
                fig_lor_fit,ax_lor_fit= plt.subplots(1,3, figsize=(18,6), constrained_layout=True)
                fig_lor_fit.canvas.manager.set_window_title('Lorenztian Fitting results')
            
                for n,clr2 in enumerate(clr):
                    if n==int(len(clr)/2):
                        ax_lor_fit[0].plot(Shift_rng1,Intensity_1rng[:,n],color=clr2, label='Data')
                        ax_lor_fit[0].plot(Shift_rng1,np.transpose(I_G_fit)[:,n],'--r', label='Lorentzian Fit');
            
                        ax_lor_fit[1].plot(Shift_rng2,Intensity_2rng[:,n],color=clr2, label='Data')
                        ax_lor_fit[1].plot(Shift_rng2,np.transpose(I_D_fit)[:,n],'--r', label='Lorentzian Fit');
    
                        ax_lor_fit[2].plot(Shift_rng3,Intensity_3rng[:,n],color=clr2, label='Data')
                        ax_lor_fit[2].plot(Shift_rng3,np.transpose(I_2D_fit)[:,n],'--r', label='Lorentzian Fit');
                    else:
                        ax_lor_fit[0].plot(Shift_rng1,Intensity_1rng[:,n],color=clr2)
                        ax_lor_fit[0].plot(Shift_rng1,np.transpose(I_G_fit)[:,n],'--r');
            
                        ax_lor_fit[1].plot(Shift_rng2,Intensity_2rng[:,n],color=clr2)
                        ax_lor_fit[1].plot(Shift_rng2,np.transpose(I_D_fit)[:,n],'--r');
    
                        ax_lor_fit[2].plot(Shift_rng3,Intensity_3rng[:,n],color=clr2)
                        ax_lor_fit[2].plot(Shift_rng3,np.transpose(I_2D_fit)[:,n],'--r');
                    
                        
                    
                ax_lor_fit[1].set_title('Fitting results D',fontsize=fs)  
                ax_lor_fit[2].set_title('Fitting results 2D',fontsize=fs)  


    
        for axx in ax_lor_fit:
            axx.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
            axx.set_ylabel('Intensity / a.u',fontsize=fs)  
            axx.legend(fontsize=fs)
        ax_lor_fit[0].set_title('Fitting results G',fontsize=fs)
        fig_lor_fit.suptitle('Lorentzian fits of '+name[z],fontsize=fs+3)

        fig_lor_fit.savefig('Lorentz_fit_'+name[z]+imgtype)
#%% #Raman shift RBM
    if rbm==1:

        if z==0:
            fig_RBM,ax_RBM= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_RBM.canvas.manager.set_window_title('Raman shift RBM modes')
        
        npeaks=max([len(n) for n in PeaksLoc1])
        PeaksLoc  = np.array(flatten(PeaksLoc1))

        binsShift=np.arange(min(PeaksLoc)-width, max(PeaksLoc) + width, width)
  
        [histRBM, edgesRBM]=np.histogram(PeaksLoc,binsShift)
        peaksRBM, _ = find_peaks(np.concatenate(([min(histRBM)],histRBM,[min(histRBM)])), distance=10)
        peaksRBM=peaksRBM-1
        f_peaksRBM=edgesRBM[peaksRBM]
        
        groups=dict.fromkeys(f_peaksRBM)
        meanRBM[name[z]]=dict.fromkeys(f_peaksRBM)
        stdRBM[name[z]]=dict.fromkeys(f_peaksRBM)

        for n in groups:
            groups[n]=[]
        for p in PeaksLoc:
            try:
                diff_p=abs(p-f_peaksRBM)
                ni=np.where(diff_p==min(diff_p))[0][0]
                groups[f_peaksRBM[ni]].append(p)
            except ValueError:
                continue
        
        for n in groups:
            meanRBM[name[z]][n]=np.mean(groups[n])
            stdRBM[name[z]][n]=np.std(groups[n])
        
        lb=[str(round(meanRBM[name[z]][n],2))+ ' $\pm$ ' +str(round(stdRBM[name[z]][n],2)) for n in stdRBM[name[z]]]
            
        ax_RBM.hist(PeaksLoc,binsShift,color=clr[round(len(clr)/2)],label=name[z]+': '+', '.join(map(str, lb))+' $cm^{-1}$',
                      alpha=0.5,ec='k',align='left',density=den)
        
        ax_RBM.legend(fontsize=fs-2)  
        ax_RBM.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_RBM.set_ylabel('counts',fontsize=fs)
        ax_RBM.tick_params(axis="both", labelsize=fs)   
        ax_RBM.set_title('Raman shift $RBM$ modes',fontsize=fs+2);
        ax_RBM.xaxis.set_major_locator(MaxNLocator(integer=True))
        fig_peaks.savefig('Peak Data.svg')
        fig_RBM.savefig('RBM.svg')




#%% #Correlations

    if z==0:
        fig_Corr,ax_Corr= plt.subplots(1,3, figsize=(15,5), constrained_layout=True)
        fig_Corr.canvas.manager.set_window_title('Correlations')

    if nt==1:# In the case where there are separate G+/G- peaks, use G+
        #Linear fit Intensity Ratio vs G+
        fit_I21_G=np.polyfit(I21,center_Gplus,1)
        #Plot Intensity Ratio vs G+
        ax_Corr[0].plot(I21,center_Gplus,marker='o',ls='',color=clr[round(len(clr)/2)],label=name[z]+', slope='+str(round(fit_I21_G[0],3)))
      
        
        #Linear fit 2D vs Gplus
        fit_2D_G=np.polyfit(center_Gplus,center_2D,1)
        #Plot 2D vs G
        ax_Corr[2].plot(center_Gplus,center_2D,marker='o',ls='',color=clr[round(len(clr)/2)],label=name[z]+', slope='+str(round(fit_2D_G[0],3)))
        xx2=np.linspace(min(center_Gplus),max(center_Gplus),100)

    else:
        fit_I21_G=np.polyfit(I21,center_G,1)
        ax_Corr[0].plot(I21,center_G,marker='o',ls='',color=clr[round(len(clr)/2)],label=name[z]+', slope='+str(round(fit_I21_G[0],3)))
        
        #Linear fit 2D vs G
        fit_2D_G=np.polyfit(center_G,center_2D,1)
        #Plot 2D vs G
        ax_Corr[2].plot(center_G,center_2D,marker='o',ls='',color=clr[round(len(clr)/2)],label=name[z]+', slope='+str(round(fit_2D_G[0],3)))
        xx2=np.linspace(min(center_G),max(center_G),100)


    xx=np.linspace(min(I21),max(I21),100)
    ax_Corr[0].plot(xx,fit_I21_G[0]*xx+fit_I21_G[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
    fit_I21_D=np.polyfit(I21,center_D,1)    
    ax_Corr[1].plot(I21,center_D,marker='o',ls='',color=clr[round(len(clr)/2)],label=name[z]+', slope='+str(round(fit_I21_D[0],3)))
    ax_Corr[1].plot(xx,fit_I21_D[0]*xx+fit_I21_D[1],lw=2,ls='-',color=clr[round(len(clr)/3)])




    ax_Corr[2].plot(xx2,fit_2D_G[0]*xx2+fit_2D_G[1],lw=2,ls='-',color=clr[round(len(clr)/3)])

        
    xlab=['$I_{d}/I_{g}$','$I_{d}/I_{g}$','Raman Shift $G$ band /$ cm^{-1}$']
    ylab=['Raman Shift $G$ band /$ cm^{-1}$','Raman Shift $D$ band /$ cm^{-1}$','Raman Shift $2D$ band /$ cm^{-1}$']
    titles=['$G$ band Shift vs. $I_{d}/I_{g}$','$D$ band Shift vs. $I_{d}/I_{g}$','$2D$ band Shift vs. $G$ band Shift']
    for i,ax in enumerate(ax_Corr):
        ax.legend(fontsize=fs-3)
        ax.tick_params(axis="both", labelsize=fs)
        ax.set_xlabel(xlab[i],fontsize=fs)
        ax.set_ylabel(ylab[i],fontsize=fs)
        ax.set_title(titles[i],fontsize=fs)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        

#%% #2D Map

    if maps==1:
        if use_leng!=1:
            ratio=col/rows
        else:
            ratio=xlen/ylen
        try:
            fig_Map.clf()
        except NameError:
            pass
        fig_Map=plt.figure('2D Map', figsize=(7*ratio,6), constrained_layout=True)
        ax_Map= fig_Map.add_subplot(111)
        if map_var=='G':
            if nt==1:
                var2D=np.array(center_Gplus)
            else:
                var2D=center_G

        elif map_var=='I':
            var2D=I21
        elif map_var=='D':
            var2D=center_D
        elif map_var=='2D':
            var2D=center_2D
        else:
            print('map_var must be I, G, D, or 2D')
            break
        
        spec_label={'I':'$I_{d}/I_{g}$',
                    'G':'$G$ band',
                    'D':'$D$ band',
                    '2D':'$2D$ band'}
        if map_var=='I':
            unit=''
        else:
            unit=' Raman Shift /$ cm^{-1}$'
        
        G2D=var2D.reshape([rows,col])
        
        if use_leng!=1:
            y1=np.linspace(0,rows-1,rows)
            x1=np.linspace(0,col-1,col)
            ylab='y pixel'
            xlab='x pixel'
        else:
            y1=np.linspace(0,ylen,rows)
            x1=np.linspace(0,xlen,col)
            ylab='y length, ($\mu m$)'
            xlab='x length, ($\mu m$)'        
        
        h=ax_Map.pcolormesh(x1,y1,G2D)
        clbr=plt.colorbar(h,ax=ax_Map)
        clbr.set_label(spec_label[map_var] +unit, rotation=270,labelpad=20,fontsize=fs)
        clbr.ax.tick_params(labelsize=fs) 
        ax_Map.set_xlabel(xlab,fontsize=fs)
        ax_Map.set_ylabel(ylab,fontsize=fs)
        ax_Map.set_title('2D Map of '+spec_label[map_var]+' for '+name[z],fontsize=fs)

        ax_Map.tick_params(axis="both", labelsize=fs)
        sz=rows*col
        pixelwidthx=xlen/(2*int(col))
        pixelwidthy=ylen/(2*int(rows))
        x1a=np.linspace(0-pixelwidthx,xlen+pixelwidthx,int(col)+1)
        y1a=np.linspace(0-pixelwidthy,ylen+pixelwidthy,int(rows)+1)
        xx,y1=np.meshgrid(x1,y1)
        xx=xx.flatten()
        y1=y1.flatten()
    
        fig_Map.savefig(name[z]+'_2DMap_'+map_var+imgtype)

    
    
        pkn=5
        line, = ax_Map.plot(xx, y1, 'b',picker=pkn)
        line.set_visible(False)

     


#%% Save shifts in .csv 


    if lorentz==0 and nt==1:
        data_all={'Intensity_G-':Int_Gmin,'Shift_G-':center_Gmin,'FWHM_G-':FWHM_Gmin,'Intensity_G+':Int_Gplus,
                  'Shift_G+':center_Gplus,'FWHM_G+':FWHM_Gplus,'Intensity_D':Int_D,'Shift_D':center_D,'Intensity_2D':Int_2D,'Shift_2D':center_2D}
    elif lorentz==0 and nt == 0:
        data_all={'Intensity_G':Int_G,'Shift_G':center_G,'Intensity_D':Int_D,'Shift_D':center_D,'Intensity_2D':Int_2D,'Shift_2D':center_2D}
    elif lorentz==1 and nt==1:        
        data_all={'Intensity_G-':Int_Gmin,'Shift_G-':center_Gmin,'FWHM_G-':FWHM_Gmin,'Intensity_G+':Int_Gplus,
                  'Shift_G+':center_Gplus,'FWHM_G+':FWHM_Gplus,'Intensity_D':Int_D,'Shift_D':center_D,
                  'FWHM_D':FWHM_D,'Intensity_2D':Int_2D,'Shift_2D':center_2D,'FWHM_2D':FWHM_2D}
    elif lorentz==1 and nt==0:        
        data_all={'Intensity_G':Int_G,'Shift_G':center_G,'FWHM_G':FWHM_G,'Intensity_D':Int_D,
                  'Shift_D':center_D,'FWHM_D':FWHM_D,'Intensity_2D':Int_2D,'Shift_2D':center_2D,'FWHM_2D':FWHM_2D}       
        
    data_all['Intensity Ratio']=I21  
    
   # if rbm==1:
    #    data_all['Intensity_RBM']=PeaksInt
     #   data_all['Shift_RBM']=PeaksLoc
        
    data_all_df=pd.DataFrame(data_all)
    data_all_df.index.name='Spectra #'
    data_all_df.to_csv(name[z]+'_results.csv')


#%%% Save Figures in Scaleable Vector Format
if raw==1:
    fig_raw.savefig('Raw Data'+imgtype)
if norm==1:
    fig_norm.savefig('Normalized Data'+imgtype)
    
fig_avg.savefig('Avg Spectra'+imgtype)

if rng==1:
    fig_rng.savefig('SpecRangeIR'+imgtype)

fig_IR.savefig('IntesityRatio'+imgtype)
fig_2D.savefig('2DBand'+imgtype)
fig_G.savefig('Gband'+imgtype)
fig_D.savefig('Dband'+imgtype)

if correlations==1:
    fig_Corr.savefig('Correlations'+imgtype)


if nt==1:
    fig_G_I.savefig('IntesityRatioG+G-'+imgtype)
    

#%% Pick Event for 2D Map

def onpick(event):
    global pts
    global ptsM
    global ind
    global figFC
    global axFC
    global data_spec
    
    if 'ind' in globals():
        del ind
    
    
    try:
        pts
    except NameError:
        pts = None
    
    if pts!=None:
        if plt.fignum_exists('Spectra'):    
            pass
        else:
            for p in pts:
                for pp in p:
                    pp.remove()

    if plt.fignum_exists('Spectra'):    
        pass
    else:
        figFC = plt.figure(num='Spectra',figsize=(10*1.5,5*1.5))
        axFC = figFC.add_subplot(111)
        data_spec=pd.DataFrame()
        data_spec['Raman Shift(cm-1)']=Shift
        pts=[]
        
           
    if event.artist!=line:  #check that you clicked on the object you wanted
        return True     
    if not len(event.ind):  #check the index is valid
        return True
    
    ind = event.ind[0]
    
    if map_var=='I':
        unit=''
    else:
        unit=' $cm^{-1}$'
        
    
    pts.append(ax_Map.plot(xx[ind],y1[ind],'ro'))

    if use_leng!=1:
        labp='x:'+str(int(xx[ind]))+', y:'+str(int(y1[ind]))
    else:
        labp='x:'+str(round_sig(xx[ind],3))+'$\mu m$, y:'+str(round_sig(y1[ind],3))+'$\mu m$'
        

    data_spec[labp]=Intensity_norm[:,ind]
     
    #Plot spectra for selected pixel
    if use_leng!=1:
        axFC.plot(Shift,Intensity_norm[:,ind],
                  label=labp+', '+spec_label[map_var]+': '+str(round_sig(var2D[ind],5))+unit,
                  markersize=1)
        print('\nPixel Location x:',int(xx[ind]),' y:',int(y1[ind]),'\n')

    else:
        axFC.plot(Shift,Intensity_norm[:,ind], label=labp+', '+spec_label[map_var]+': '+str(round_sig(var2D[ind],5))+unit, markersize=1)
        print('\nPixel Location :'+str(round_sig(xx[ind],3))+' $\mu m$, y:'+str(round_sig(y1[ind],3))+' $\mu m$ \n')

    axFC.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
    axFC.set_ylabel('Intensity',fontsize=fs)
    axFC.legend(fontsize=fs-1)
    axFC.set_title('Selected Normalized Spectra on '+name[z],fontsize=fs+1)

    axFC.tick_params(axis="both", labelsize=fs)   
    print('Intensity Ratio is '+str(round_sig(I21[ind],5))+'\n')
    if nt==1:
        print('G+ band Shift is '+str(round_sig(center_Gplus[ind],5))+' cm-1 \n')
    else:
        print('G band Shift is '+str(center_G[ind])+' cm-1 \n')
    print('D band Shift is '+str(center_D[ind])+' cm-1 \n')
    print('2D band Shift is '+str(center_2D[ind])+' cm-1  \n')
    
    data_spec.to_csv('Selected_Spectra.csv',index=False)
    fig_Map.canvas.draw() 
    figFC.canvas.draw() 
    figFC.savefig('SelectedSpectra'+imgtype)
    return True


    
if maps==1:    
    fig_Map.canvas.mpl_connect('pick_event', onpick)    

