# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:48:31 2023

@author: matte
"""

# load packages
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

from tqdm import tqdm

from scipy.signal import find_peaks

from analysisInputs import analysisInputs 
from RBMinputs import RBMinputs
from LorentzInputsDouble1 import LorentzInputsDouble1
from LorentzInputs import LorentzInputs
from PlotsInputs import PlotsInputs

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
#%%%%
#Process Data defined in GUI

def process_data(folder_selected, selected_files, file_name, labels, delim):
    # Perform your data processing operations using the passed variables
    total=len(selected_files)
    
    for z in range(0,total):

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#Import data #%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%

        data = pd.read_csv(selected_files[z],delimiter=delim,decimal='.',header=None);
        data=data.to_numpy()
        Shift = np.array(data[:,0]); #Raman shift is the first column of the file
        Intensity=data[:,1:]; #Raman intensities are all data from column 2 til end
    
        Intensity_av_raw=np.mean(Intensity,axis=1); #Calculate average spectra
        
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
            
        
     #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#   
     #%%#%%#%%#Average Raw Spectra #%%#%%#%%#%%
     #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%       
     
        if z==0:
           fig_avg_raw,ax_avg_raw= plt.subplots(div[0],div[1], figsize=(div[1]*5,div[0]*5), constrained_layout=True)
           fig_avg_raw.canvas.manager.set_window_title('Average Raw Spectra')
           if total==1:
               ax_avg_raw=[ax_avg_raw]
           else:
               ax_avg_raw=ax_avg_raw.flatten()  
               
        ax_avg_raw[z].plot(Shift,Intensity_av_raw,color=clr[round(len(clr)/2)]);
        
        fs=16
        ax_avg_raw[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_avg_raw[z].set_ylabel('Intensity / a.u',fontsize=fs)
        ax_avg_raw[z].tick_params(axis="both", labelsize=fs)   
        ax_avg_raw[z].set_title(labels[z]+': Average raw spectra');    
        ax_avg_raw[z].set_ylim((0,np.max(Intensity_av_raw)))
        ax_avg_raw[z].set_xlim((0,np.max(Shift)))
        
    lorentz, nt, rbm, map_, cols, rows, xlen, ylen, normLow, normHigh, band1Low, band1High, band1Name, band2Low, band2High, band2Name, band3Low, band3High, band3Name, use_leng = analysisInputs()
    
    #Close current figure
    plt.close('all')
    
    for z in range(0,total):

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#Import data #%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%

        data = pd.read_csv(selected_files[z],delimiter=delim,decimal='.',header=None);
        data=data.to_numpy()
        Shift = np.array(data[:,0]); #Raman shift is the first column of the file
        Intensity=data[:,1:]; #Raman intensities are all data from column 2 til end
        color=cmap[z];
        clr=color(np.linspace(0,1,len(Intensity.transpose())))
        
        
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
        

        if z==0:
           fig_avg,ax_avg= plt.subplots(div[0],div[1], figsize=(div[1]*5,div[0]*5), constrained_layout=True)
           fig_avg.canvas.manager.set_window_title('Average Normalized Spectra')
           if total==1:
               ax_avg=[ax_avg]
           else:
               ax_avg=ax_avg.flatten()  
               
        ax_avg[z].plot(Shift,Intensity_av,color=clr[round(len(clr)/2)]);
        
        fs=16
        ax_avg[z].set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_avg[z].set_ylabel('Intensity / a.u',fontsize=fs)
        ax_avg[z].tick_params(axis="both", labelsize=fs)   
        ax_avg[z].set_title(labels[z]+': Average normalized spectra');    
        ax_avg[z].set_ylim((0,np.max(Intensity_av)))
        ax_avg[z].set_xlim((0,np.max(Shift)))
         
#%%Additonal User Interfaces for more inputs (RBMs, Lorentzian fitting parameters, plotting options)


    # RBM inputs (if RBM=1)
    if rbm == 1:
        RBMregion_Low,RBMregion_High,Prom=RBMinputs(total)
    
    #Fitting parameters (only if lorentz=1)
    if lorentz == 1:
        if nt == 1:
            Init_peak1plus,Init_peak1_min,Init_peak2,Init_peak3= LorentzInputsDouble1(band1Name,band2Name,band3Name,total)
        else:
            Init_peak1plus,Init_peak2,Init_peak3= LorentzInputs(band1Name,band2Name,band3Name,total)
    else:
        nt=0 # You can only have the split peaks option if a lorentzian fit is selected

    # Ploting options/Output plots:
    width, width_fw, width_int, raw, norm, rng, peaks, correlation1, correlation2, correlation3, correlation4, correlation5, map1, map2, map3, map4, saveFigs, nameFigs, imgtype, dens, fs, maxI21 = PlotsInputs(band1Name,band2Name,band3Name)
    
    #Create folder, if it doesnt exist for output data
    if not os.path.exists(folder_selected + '/'+nameFigs+'_Results'):
        os.makedirs(folder_selected + '/'+nameFigs+'_Results')
        
    #Prevent user from inputing 0 for histogram widths
    if width==0:
        width=3; 
    if width_fw==0:
        width_fw=5; 
    if width_int==0:
        width_int=0.1
    
    if lorentz==0:
        map1 = 0
        map2 = 0
        map3 = 0
        
        
    # Close all figures
    plt.close('all')
    
         
#%%ANALYSIS CODE STARTS HERE
    for z in range(0,total):

#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#Import data #%%#%%#%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%%

        data = pd.read_csv(selected_files[z],delimiter=delim,decimal='.',header=None);
        data=data.to_numpy()
        Shift = np.array(data[:,0]); #Raman shift is the first column of the file
        Intensity=data[:,1:]; #Raman intensities are all data from column 2 til end
        color=cmap[z];
        clr=color(np.linspace(0,1,len(Intensity.transpose())))
        
        
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
    
    ###  Define spectral range for G, D and 2D
        ind1Low=np.where(abs(Shift-band1Low)==min(abs(Shift-band1Low)))[0][0]
        ind1High=np.where(abs(Shift-band1High)==min(abs(Shift-band1High)))[0][0];
        
        ind2Low=np.where(abs(Shift-band2Low)==min(abs(Shift-band2Low)))[0][0]
        ind2High=np.where(abs(Shift-band2High)==min(abs(Shift-band2High)))[0][0];
        
        ind3Low=np.where(abs(Shift-band3Low)==min(abs(Shift-band3Low)))[0][0]
        ind3High=np.where(abs(Shift-band3High)==min(abs(Shift-band3High)))[0][0];
        
    ### Define Intensity and shift vectors
        Intensity_1rng=Intensity_norm[ind1Low:ind1High,:];
        Intensity_2rng=Intensity_norm[ind2Low:ind2High,:];
        Intensity_3rng=Intensity_norm[ind3Low:ind3High,:];
        
        Shift_rng1=Shift[ind1Low:ind1High];
        Shift_rng2=Shift[ind2Low:ind2High];
        Shift_rng3=Shift[ind3Low:ind3High];
        
        Int_1plus=np.array([]);
        Int_2=np.array([]);
        Int_3=np.array([]);
        
        center_1plus=np.array([]);
        center_2=np.array([]);
        center_3=np.array([]);
        
        FWHM_1plus=np.array([]);
        FWHM_2=np.array([]);
        FWHM_3=np.array([]);
        
        
    ### Initialize Lorentztian fitting varables, if selected
        if lorentz==1:
            I_2_fit=[]
            I_3_fit=[]
            if nt==1:
                I_1_fit=[]
                center_1min=[]
                Int_1min=[]
                FWHM_1min=[]
            else:
                I_1_fit=[]
                
                
        for n in tqdm(range(0,len(Intensity_norm[0,:]))):
            
            
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
    #%%#%%#%%#%% Method 1 (Find bands using Max Intensity) #%%#%%#%%#%%#%%
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%            
           
           if lorentz==0:
               #peak 1/G: 
               Int_1plus=np.append(Int_1plus,max(Intensity_1rng[:,n])-min(Intensity_1rng[:,n]));
               indMax1=np.where(Intensity_1rng[:,n]==max(Intensity_1rng[:,n]))[0][0];
               center_1plus=np.append(center_1plus,Shift_rng1[indMax1]);
               
               #peak 2/D: 
               Int_2=np.append(Int_2,max(Intensity_2rng[:,n])-min(Intensity_2rng[:,n]));   
               indMax2=np.where(Intensity_2rng[:,n]==max(Intensity_2rng[:,n]))[0][0];
               center_2=np.append(center_2,Shift_rng2[indMax2]);
               
               #peak 3/2D: 
               Int_3=np.append(Int_3,max(Intensity_3rng[:,n])-min(Intensity_3rng[:,n]));
               indMax3=np.where(Intensity_3rng[:,n]==max(Intensity_3rng[:,n]))[0][0];
               center_3=np.append(center_3,Shift_rng3[indMax3]);
               
               
               
               
   
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
    #%%#%%#%%#%% Method 2 (Lorentzian Fitting) #%%#%%#%%#%%#%%
    #%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%% 
                     
           if lorentz==1:
               
               #D band / band #2
               #Define Initial guesses
               InitGuess_2=[Init_peak2[z][2]*((Init_peak2[z][1]/2)**2),
                            Init_peak2[z][0],
                            (Init_peak2[z][1]/2)**2,
                            0];
               
               def fit_func2(x, g1,g2,g3,g4):
                   return (g1/((x-g2)**2+g3)+g4)# lorentz
               
               gamma_2 ,pcov2= curve_fit(fit_func2, Shift_rng2, Intensity_2rng[...,n], 
                                          InitGuess_2,bounds=([0,Shift[ind2Low],0,0],
                                        [np.inf,Shift[ind2High],0.25*(Shift[ind2High]-Shift[ind2Low])**2,np.inf]),
                                        maxfev=10000)
               
               
               I_2_fit.append(fit_func2(Shift_rng2, *gamma_2))
               
               if r2_score(I_2_fit[n],Intensity_2rng[...,n])>=0.5:
                   center_2=np.append(center_2,gamma_2[1]);
                   FWHM_2=np.append(FWHM_2,2*np.sqrt(gamma_2[2]));
                   Int_2=np.append(Int_2,gamma_2[0]/gamma_2[2]); 
               else:
                   center_2=np.append(center_2,np.nan);
                   FWHM_2=np.append(FWHM_2,np.nan);
                   Int_2=np.append(Int_2,np.nan); 
        
               #2D band / band #3
               #Define Initial guesses
               InitGuess_3=[Init_peak3[z][2]*((Init_peak3[z][1]/2)**2),
                            Init_peak3[z][0],
                            (Init_peak3[z][1]/2)**2,
                            0];    
          
               
               gamma_3 ,pcov2= curve_fit(fit_func2, Shift_rng3, Intensity_3rng[...,n], 
                                          InitGuess_3, bounds=([0,Shift[ind3Low],0,0],
                                        [np.inf,Shift[ind3High],0.25*(Shift[ind3High]-Shift[ind3Low])**2,np.inf]),
                                        maxfev=10000)
               
               I_3_fit.append(fit_func2(Shift_rng3, *gamma_3))
               
               if r2_score(I_3_fit[n],Intensity_3rng[...,n])>=0.5:
    
                   center_3=np.append(center_3,gamma_3[1]);
                   FWHM_3=np.append(FWHM_3,2*np.sqrt(gamma_3[2]));
                   Int_3=np.append(Int_3,gamma_3[0]/gamma_3[2]);  
               else:
                   center_3=np.append(center_3,np.nan);
                   FWHM_3=np.append(FWHM_3,np.nan);
                   Int_3=np.append(Int_3,np.nan);      
                   
       #%%G+ and G- fitting   
               if nt==1:
                   
                   #Define initial guesses
                   InitGuess_1=[Init_peak1_min[z][2]*((Init_peak1_min[z][1]/2)**2),
                                Init_peak1_min[z][0],
                                (Init_peak1_min[z][1]/2)**2,
                                Init_peak1plus[z][2]*((Init_peak1plus[z][1]/2)**2),
                                Init_peak1plus[z][0],
                                (Init_peak1plus[z][1]/2)**2,
                                0.1]
                   def fit_func(x, g1,g2,g3,g4,g5,g6,g7):
                       return (g1/((x-g2)**2+g3)+g4/((x-g5)**2+g6)+g7)# lorentz
        
                   gamma_1 ,pcov2= curve_fit(fit_func, Shift_rng1, Intensity_1rng[...,n], InitGuess_1, 
                                             bounds=([0,Shift[ind1Low],0,0,Shift[ind1Low],0,0],
                                           [np.inf,Shift[ind1High],0.25*(Shift[ind1High]-Shift[ind1Low])**2,np.inf,
                                            Shift[ind1High],0.25*(Shift[ind1High]-Shift[ind1Low])**2,np.inf]),
                                           maxfev=10000)
                                       
                   I_1_fit.append(fit_func(Shift_rng1, *gamma_1))
                   
                   if r2_score(I_1_fit[n],Intensity_1rng[...,n])>=0.5:
                       if gamma_1[1]<gamma_1[4]:
                              center_1min=np.append(center_1min,gamma_1[1])
                              FWHM_1min=np.append(FWHM_1min,2*np.sqrt(gamma_1[2]))
                              Int_1min=np.append(Int_1min,gamma_1[0]/gamma_1[2])
                              center_1plus=np.append(center_1plus,gamma_1[4])
                              FWHM_1plus=np.append(FWHM_1plus,2*np.sqrt(gamma_1[5]))
                              Int_1plus=np.append(Int_1plus,gamma_1[3]/gamma_1[5])
                       else:
                             center_1plus=np.append(center_1plus,gamma_1[1])
                             FWHM_1plus=np.append(FWHM_1plus,2*np.sqrt(gamma_1[2]))
                             Int_1plus=np.append(Int_1plus,gamma_1[0]/gamma_1[2])
                             center_1min=np.append(center_1min,gamma_1[4])
                             FWHM_1min=np.append(FWHM_1min,2*np.sqrt(gamma_1[5]))
                             Int_1min=np.append(Int_1min,gamma_1[3]/gamma_1[5])
                   else:
                             center_1plus=np.append(center_1plus,np.nan)
                             FWHM_1plus=np.append(FWHM_1plus,np.nan)
                             Int_1plus=np.append(Int_1plus,np.nan)
                             center_1min=np.append(center_1min,np.nan)
                             FWHM_1min=np.append(FWHM_1min,np.nan)
                             Int_1min=np.append(Int_1min,np.nan)
                             
                    
                             
               else:
                    #Lorentzian G fit if no G+/G- fit
                    InitGuess_1=[Init_peak1plus[z][2]*((Init_peak1plus[z][1]/2)**2),
                                 Init_peak1plus[z][0],
                                 (Init_peak1plus[z][1]/2)**2,
                                 0];
                   
                    
                    gamma_1 ,pcov2= curve_fit(fit_func2, Shift_rng1, Intensity_1rng[...,n], 
                                     InitGuess_1, bounds=([0,Shift[ind1Low],0,0],
                                   [np.inf,Shift[ind1High],0.25*(Shift[ind1High]-Shift[ind1Low])**2,np.inf]),
                                   maxfev=10000)
                    
                    I_1_fit.append(fit_func2(Shift_rng1, *gamma_1))
                    
                    if r2_score(I_1_fit[n],Intensity_1rng[...,n])>=0.5:
                        center_1plus=np.append(center_1plus,gamma_1[1]);
                        FWHM_1plus=np.append(FWHM_1plus,2*np.sqrt(gamma_1[2]));
                        Int_1plus=np.append(Int_1plus,gamma_1[0]/gamma_1[2]); 
                    else:
                        center_1plus=np.append(center_1plus,np.nan);
                        FWHM_1plus=np.append(FWHM_1plus,np.nan);
                        Int_1plus=np.append(Int_1plus,np.nan);   
                
                
        #Calculate mean values and standard deviations       
        if nt==1:
            
            #If peaks are split, use peak with max intensity
            maxInt_1=np.maximum(Int_1plus,Int_1min)
            I21=Int_2/maxInt_1
            I31=Int_3/maxInt_1

            center_1min[I21>maxI21]=np.nan
            Int_1min[I21>maxI21]=np.nan


  
            peak1min_av=np.nanmean(center_1min)
            peak1min_std=np.nanstd(center_1min)  
        else:
            I21=Int_2/Int_1plus;
            I31=Int_3/Int_1plus
            
        Int_1plus[I21>maxI21]=np.nan
        Int_2[I21>maxI21]=np.nan
        Int_3[I21>maxI21]=np.nan
       
        center_1plus[I21>maxI21]=np.nan
        center_2[I21>maxI21]=np.nan
        center_3[I31>maxI21]=np.nan
        I21[I21>maxI21]=np.nan    
        I31[I31>maxI21]=np.nan    

        peak1plus_av=np.nanmean(center_1plus)
        peak1plus_std=np.nanstd(center_1plus)              
        I21_av=np.nanmean(I21);
        I21_error=np.nanstd(I21);
        I31_av=np.nanmean(I31);
        I31_error=np.nanstd(I31);
        peak2_av=np.nanmean(center_2)
        peak2_std=np.nanstd(center_2)
        peak3_av=np.nanmean(center_3)
        peak3_std=np.nanstd(center_3)
        
        if lorentz == 1:
            FWHM_1plus_av=np.nanmean(FWHM_1plus);
            FWHM_1plus_std=np.nanstd(FWHM_1plus);
            FWHM_3_av=np.nanmean(FWHM_3);
            FWHM_3_std=np.nanstd(FWHM_3);
            FWHM_2_av=np.nanmean(FWHM_2);
            FWHM_2_std=np.nanstd(FWHM_2);  
            if nt == 1:
                FWHM_1min_av=np.nanmean(FWHM_1min);
                FWHM_1min_std=np.nanstd(FWHM_1min);
                
                
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
                    ax_peaks[z].set_title(labels[z]+' :Region Peaks');
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
           ax_raw[z].set_title(labels[z]+': Raw data');
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
            ax_norm[z].set_title(labels[z]+': Normalized spectra');    
            ax_norm[z].set_ylim((0,np.max(Intensity_norm)))
            ax_norm[z].set_xlim((0,np.max(Shift)))                               
            

    #%%Average spectra
        if z==0:
           fig_avg,ax_avg= plt.subplots(1,1, figsize=(6,6), constrained_layout=True)
           fig_avg.canvas.manager.set_window_title('Average Spectra')
               
        
        ax_avg.plot(Shift,Intensity_av,color=clr[round(len(clr)/2)],label=labels[z]);
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
            ax_rng[z].set_title(labels[z]+': Range for $I_{'+band2Name+'}/I_{'+band1Name+'}$ Calculation');    
            ax_rng[z].set_ylim((0,1))

    #%%%# %% Lorentz fit results: plots the spectra together with the best fit
        if lorentz==1:
                    
            fig_lor_fit,ax_lor_fit= plt.subplots(1,3, figsize=(18,6), constrained_layout=True)
            fig_lor_fit.canvas.manager.set_window_title(labels[z]+' Lorenztian Fitting results')
        
            for n,clr2 in enumerate(clr):
                if n==int(len(clr)/2):
                    ax_lor_fit[0].plot(Shift_rng1,Intensity_1rng[:,n],color=clr2, label='Data')
                    ax_lor_fit[0].plot(Shift_rng1,np.transpose(I_1_fit)[:,n],'--r', label='Lorentzian Fit');
        
                    ax_lor_fit[1].plot(Shift_rng2,Intensity_2rng[:,n],color=clr2, label='Data')
                    ax_lor_fit[1].plot(Shift_rng2,np.transpose(I_2_fit)[:,n],'--r', label='Lorentzian Fit');

                    ax_lor_fit[2].plot(Shift_rng3,Intensity_3rng[:,n],color=clr2, label='Data')
                    ax_lor_fit[2].plot(Shift_rng3,np.transpose(I_3_fit)[:,n],'--r', label='Lorentzian Fit');
                else:
                    ax_lor_fit[0].plot(Shift_rng1,Intensity_1rng[:,n],color=clr2)
                    ax_lor_fit[0].plot(Shift_rng1,np.transpose(I_1_fit)[:,n],'--r');
        
                    ax_lor_fit[1].plot(Shift_rng2,Intensity_2rng[:,n],color=clr2)
                    ax_lor_fit[1].plot(Shift_rng2,np.transpose(I_2_fit)[:,n],'--r');

                    ax_lor_fit[2].plot(Shift_rng3,Intensity_3rng[:,n],color=clr2)
                    ax_lor_fit[2].plot(Shift_rng3,np.transpose(I_3_fit)[:,n],'--r');
                    
            fig_lor_fit.suptitle('Lorentzian fits of '+labels[z],fontsize=fs+3)            
            ax_lor_fit[0].set_title('Fitting results '+band1Name,fontsize=fs)               
            ax_lor_fit[1].set_title('Fitting results '+band2Name,fontsize=fs)  
            ax_lor_fit[2].set_title('Fitting results '+band3Name,fontsize=fs)  

            for axx in ax_lor_fit:
                axx.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
                axx.set_ylabel('Intensity / a.u',fontsize=fs)  
                axx.legend(fontsize=fs)
            
            #Save Lorentzian Fit for each file
            os.chdir(folder_selected + '/'+nameFigs+'_Results/')
            fig_lor_fit.savefig('Lorentz_fit_'+labels[z]+imgtype)
            
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%
#%%#%%#%%Histograms #%%#%%#%%#%%%
#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%#%%

#%% Histogram intensity ratio
    
        #Create Figure 
        if z==0:
            fig_IR,ax_IR= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_IR.canvas.manager.set_window_title('Intensity Ratio Peak2/Peak1')
            fig_IR2,ax_IR2= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_IR2.canvas.manager.set_window_title('Intensity Ratio Peak3/Peak1')
    
        
    
        #Define Bin Size
        binsInt=np.arange(np.nanmin(I21)- width_int/2, np.nanmax(I21) + width_int/2, width_int)
             
        ax_IR.hist(I21,binsInt,color=clr[round(len(clr)/2)],
                   label=labels[z]+': $I_{'+band2Name+'}/I_{'+band1Name+'}$='+str(round(I21_av,4))+'$\pm$'
                   +str(round(I21_error,4)),alpha=0.5,ec='k',align='left',density=dens)
    
        ax_IR.legend(fontsize=fs-2)  
        ax_IR.set_xlabel('$I_{'+band2Name+'}/I_{'+band1Name+'}$',fontsize=fs)
        ax_IR.set_ylabel('counts',fontsize=fs)
        ax_IR.tick_params(axis="both", labelsize=fs)   
        ax_IR.set_title('Intensity ratio: $I_{'+band2Name+'}/I_{'+band1Name+'}$',fontsize=fs+2);
        
        
        #Now, do the same for the Intensity ration between band 3 and band 1
        binsInt=np.arange(np.nanmin(I31)- width_int/2, np.nanmax(I31) + width_int/2, width_int)
             
        ax_IR2.hist(I31,binsInt,color=clr[round(len(clr)/2)],
                   label=labels[z]+': $I_{'+band3Name+'}/I_{'+band1Name+'}$='+str(round(I31_av,4))+'$\pm$'
                   +str(round(I31_error,4)),alpha=0.5,ec='k',align='left',density=dens)
    
        ax_IR2.legend(fontsize=fs-2)  
        ax_IR2.set_xlabel('$I_{'+band3Name+'}/I_{'+band1Name+'}$',fontsize=fs)
        ax_IR2.set_ylabel('counts',fontsize=fs)
        ax_IR2.tick_params(axis="both", labelsize=fs)   
        ax_IR2.set_title('Intensity ratio: $I_{'+band3Name+'}/I_{'+band1Name+'}$',fontsize=fs+2);
    
    
    #%% Raman shift G mode/ peak 1
        
        if z==0:
            fig_peak1,ax_peak1= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_peak1.canvas.manager.set_window_title('Raman shift '+band1Name+' band')
            if lorentz == 1:
                fig_peak1_FWHMplus,ax_peak1_FWHMplus= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                
                if nt==1:
                    fig_peak1_I,ax_peak1_I= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                    fig_peak1_I.canvas.manager.set_window_title('Intensity ratio '+band1Name+' modes')
                    fig_peak1_FWHMmin,ax_peak1_FWHMmin= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                    fig_peak1_FWHMplus.canvas.manager.set_window_title('FWHM '+band1Name+'+ mode')
                    fig_peak1_FWHMmin.canvas.manager.set_window_title('FWHM '+band1Name+'- mode')

                else:
                    fig_peak1_FWHMplus.canvas.manager.set_window_title('FWHM '+band1Name+' band')

        if nt == 0:
            
    
            binsShift=np.arange(np.nanmin(center_1plus)-width, np.nanmax(center_1plus) + width, width)
            
            ax_peak1.hist(center_1plus,binsShift,color=clr[round(len(clr)/2)],
                          label=labels[z]+': $'+band1Name+'$ band='+str(round(peak1plus_av,2))+'$\pm$'
                          +str(round(peak1plus_std,2))+' $cm^{-1}$',
                          alpha=0.5,ec='k',align='left',density=dens)
            
            
            ax_peak1.set_title('Raman shift $'+band1Name+'$ band',fontsize=fs+2);
        else:
            binsShift=np.arange(np.nanmin(center_1min) - width/2, np.nanmax(center_1plus) + width/2, width)
            ax_peak1.hist(center_1plus,binsShift, color=clr[round(len(clr)/3)],
                      label=labels[z]+': $'+band1Name+'^{+}$ mode='+str(round(peak1plus_av,2))+'$\pm$'
                      +str(round(peak1plus_std,2))+' $cm^{-1}$',
                      alpha=0.5,ec='k',align='left',density=dens)
            ax_peak1.hist(center_1min,binsShift,color=clr[round(len(clr)/2)],
                      label=labels[z]+': $'+band1Name+'^{-}$ mode='+str(round(peak1min_av,2))+'$\pm$'
                      +str(round(peak1min_std,2))+' $cm^{-1}$',
                      alpha=0.5,ec='k',align='left',density=dens)
            ax_peak1.set_title('Raman shift $'+band1Name+'$ modes',fontsize=fs+2);
            
            
    
            Iplus_minus=np.array(Int_1plus)/np.array(Int_1min);
            Iplus_minus[Iplus_minus>maxI21]=np.nan
            binsInt=np.arange(np.nanmin(Iplus_minus)-width_int/2, np.nanmax(Iplus_minus) + width_int/2, width_int)
            ax_peak1_I.hist(Iplus_minus,binsInt,color=clr[round(len(clr)/2)],
                      label=labels[z]+': $I_{'+band1Name+'^{+}}/I_{'+band1Name+'^{-}}$ ='+str(round(np.nanmean(Iplus_minus),2))+'$\pm$'
                      +str(round(np.nanstd(Iplus_minus),2)),
                      alpha=0.5,ec='k',align='left',density=dens)
            ax_peak1_I.legend(fontsize=fs-2)  
            ax_peak1_I.set_xlabel('$I_{G^{+}}/I_{G^{-}}$',fontsize=fs)
            ax_peak1_I.set_ylabel('counts',fontsize=fs)
            ax_peak1_I.tick_params(axis="both", labelsize=fs)
            ax_peak1_I.set_title('Intensity ratio '+band1Name+' modes',fontsize=fs+2);
    
    
    
        ax_peak1.legend(fontsize=fs-2)  
        ax_peak1.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_peak1.set_ylabel('counts',fontsize=fs)
        ax_peak1.tick_params(axis="both", labelsize=fs)  
        
       # Histogram of FWHM
        if lorentz==1:
            if nt == 0:
                binsFWHM=np.arange(np.nanmin(FWHM_1plus) - width_fw/2, np.nanmax(FWHM_1plus) + width_fw/2, width_fw)
                ax_peak1_FWHMplus.hist(FWHM_1plus, binsFWHM, color=clr[round(len(clr)/2)],
                                label=labels[z]+': FWHM $'+band1Name+'$ ='+
                                str(round(FWHM_1plus_av,2))+'$\pm$'+str(round(FWHM_1plus_std,2)), 
                                alpha=0.5,ec='k',align='left',density=dens)
                ax_peak1_FWHMplus.set_title('FWHM $'+band1Name+'$ band',fontsize=fs+2);

            else:
                binsFWHM=np.arange(np.nanmin(FWHM_1plus) - width_fw/2, np.nanmax(FWHM_1plus) + width_fw/2, width_fw)
                ax_peak1_FWHMplus.hist(FWHM_1plus, binsFWHM, color=clr[round(len(clr)/2)],
                                label=labels[z]+': FWHM $'+band1Name+'^{+}$ ='+
                                str(round(FWHM_1plus_av,2))+'$\pm$'+str(round(FWHM_1plus_std,2)), 
                                alpha=0.5,ec='k',align='left',density=dens)
                
                binsFWHM=np.arange(np.nanmin(FWHM_1min) - width/2, np.nanmax(FWHM_1min) + width/2, width)
                ax_peak1_FWHMmin.hist(FWHM_1min, binsFWHM, color=clr[round(len(clr)/2)],
                                label=labels[z]+': FWHM $'+band1Name+'^{-}$ ='+
                                str(round(FWHM_1min_av,2))+'$\pm$'+str(round(FWHM_1min_std,2)), 
                                alpha=0.5,ec='k',align='left',density=dens)
                
                
                ax_peak1_FWHMmin.legend(fontsize=fs-2)  
                ax_peak1_FWHMmin.set_xlabel('FWHM / $cm^{-1}$',fontsize=fs)
                ax_peak1_FWHMmin.set_ylabel('counts',fontsize=fs)
                ax_peak1_FWHMmin.tick_params(axis="both", labelsize=fs) 
                ax_peak1_FWHMmin.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax_peak1_FWHMmin.set_title('FWHM $'+band1Name+'^{-}$ band',fontsize=fs+2);
                ax_peak1_FWHMplus.set_title('FWHM $'+band1Name+'^{+}$ band',fontsize=fs+2);

                
            ax_peak1_FWHMplus.legend(fontsize=fs-2)  
            ax_peak1_FWHMplus.set_xlabel('FWHM / $cm^{-1}$',fontsize=fs)
            ax_peak1_FWHMplus.set_ylabel('counts',fontsize=fs)
            ax_peak1_FWHMplus.tick_params(axis="both", labelsize=fs)  
            ax_peak1_FWHMplus.xaxis.set_major_locator(MaxNLocator(integer=True))

    #%% #Raman shift D mode/peak 2
    
        if z==0:
            fig_peak2,ax_peak2= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_peak2.canvas.manager.set_window_title('Raman shift '+band2Name+' band')
            if lorentz==1:
                fig_peak2_FWHM,ax_peak2_FWHM= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_peak2_FWHM.canvas.manager.set_window_title('FWHM '+band2Name+' band')
        
    
    
        binsShift=np.arange(np.nanmin(center_2)-width/2, np.nanmax(center_2) + width/2, width)
    
        ax_peak2.hist(center_2,binsShift,color=clr[round(len(clr)/2)],
                  label=labels[z]+': $'+band2Name+'$ band='+str(round(peak2_av,2))+'$\pm$'
                  +str(round(peak2_std,2))+' $cm^{-1}$',
                  alpha=0.5,ec='k',align='left',density=dens)
        
        
        ax_peak2.legend(fontsize=fs-2)  
        ax_peak2.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_peak2.set_ylabel('counts',fontsize=fs)
        ax_peak2.tick_params(axis="both", labelsize=fs)   
        ax_peak2.set_title('Raman shift $'+band2Name+'$ band',fontsize=fs+2);
        ax_peak2.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        
        
        # Histogram of FWHM
        if lorentz==1:
            binsFWHM=np.arange(np.nanmin(FWHM_2) - width/2, np.nanmax(FWHM_2) + width/2, width)
            ax_peak2_FWHM.hist(FWHM_2, binsFWHM, color=clr[round(len(clr)/2)],
                            label=labels[z]+': FWHM $'+band2Name+'$ ='+
                            str(round(FWHM_2_av,2))+'$\pm$'+str(round(FWHM_2_std,2)), 
                            alpha=0.5,ec='k',align='left',density=dens)
            
            ax_peak2_FWHM.legend(fontsize=fs-2)  
            ax_peak2_FWHM.set_xlabel('FWHM / $cm^{-1}$',fontsize=fs)
            ax_peak2_FWHM.set_ylabel('counts',fontsize=fs)
            ax_peak2_FWHM.tick_params(axis="both", labelsize=fs)
            ax_peak2_FWHM.set_title('FWHM $'+band2Name+'$ band',fontsize=fs+2);
            ax_peak2_FWHM.xaxis.set_major_locator(MaxNLocator(integer=True))

                 
                 
    #%% #Raman shift 2D mode / peak 3
        if z==0:
            fig_peak3,ax_peak3= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
            fig_peak3.canvas.manager.set_window_title('Raman shift '+band3Name+' band')
            if lorentz==1:
                fig_peak3_FWHM,ax_peak3_FWHM= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_peak3_FWHM.canvas.manager.set_window_title('FWHM '+band3Name+' band')
    
    
    
        binsShift=np.arange(np.nanmin(center_3)-width/2, np.nanmax(center_3) + width/2, width)
        
        ax_peak3.hist(center_3,binsShift,color=clr[round(len(clr)/2)],
                   label=labels[z]+': $2D$ band='+str(round(peak3_av,2))+'$\pm$'
                   +str(round(peak3_std,2))+' $cm^{-1}$',
                   alpha=0.5,ec='k',align='left',density=dens)
        
        
        ax_peak3.legend(fontsize=fs-2)  
        ax_peak3.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        ax_peak3.set_ylabel('counts',fontsize=fs)
        ax_peak3.tick_params(axis="both", labelsize=fs)   
        ax_peak3.set_title('Raman shift $'+band3Name+'$ band',fontsize=fs+2);
        ax_peak3.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        
        # Histogram of FWHM
        if lorentz==1:
            binsFWHM=np.arange(np.nanmin(FWHM_3) - width/2, np.nanmax(FWHM_3) + width/2, width)
            ax_peak3_FWHM.hist(FWHM_3, binsFWHM, color=clr[round(len(clr)/2)],
                            label=labels[z]+': FWHM $'+band3Name+'$ ='+
                            str(round(FWHM_3_av,2))+'$\pm$'+str(round(FWHM_3_std,2)), 
                            alpha=0.5,ec='k',align='left',density=dens)
            
            ax_peak3_FWHM.legend(fontsize=fs-2)  
            ax_peak3_FWHM.set_xlabel('FWHM / $cm^{-1}$',fontsize=fs)
            ax_peak3_FWHM.set_ylabel('counts',fontsize=fs)
            ax_peak3_FWHM.tick_params(axis="both", labelsize=fs)
            ax_peak3_FWHM.set_title('FWHM $'+band3Name+'$ band',fontsize=fs+2);
            ax_peak3_FWHM.xaxis.set_major_locator(MaxNLocator(integer=True))
            
                 
    #%% #Raman shift RBM Histograms
        if rbm==1:
    
            if z==0:
                fig_RBM,ax_RBM= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_RBM.canvas.manager.set_window_title('Raman shift RBM modes')
                meanRBM=dict.fromkeys(labels)
                stdRBM=dict.fromkeys(labels)    
                meanRBM_Int=dict.fromkeys(labels)
                stdRBM_Int=dict.fromkeys(labels)  
            
         #   npeaks=max([len(n) for n in PeaksLoc1])
            PeaksLoc  = np.array(flatten(PeaksLoc1))
            PeaksInt2= np.array(flatten(PeaksInt))
            
            binsShift=np.arange(np.nanmin(PeaksLoc)-width/2, np.nanmax(PeaksLoc) + width/2, width)
      
            [histRBM, edgesRBM]=np.histogram(PeaksLoc,binsShift)
            peaksRBM, _ = find_peaks(np.concatenate(([min(histRBM)],histRBM,[min(histRBM)])), distance=10)
            peaksRBM=peaksRBM-1
            f_peaksRBM=edgesRBM[peaksRBM]
            
            groups=dict.fromkeys(f_peaksRBM)
            groups_Int=dict.fromkeys(f_peaksRBM)
            meanRBM[labels[z]]=dict.fromkeys(f_peaksRBM)
            stdRBM[labels[z]]=dict.fromkeys(f_peaksRBM)
            meanRBM_Int[labels[z]]=dict.fromkeys(f_peaksRBM)
            stdRBM_Int[labels[z]]=dict.fromkeys(f_peaksRBM)
    
            for n in groups:
                groups[n]=[]
                groups_Int[n]=[]
            for jj, p in enumerate(PeaksLoc):
                try:
                    diff_p=abs(p-f_peaksRBM)
                    ni=np.where(diff_p==min(diff_p))[0][0]
                    groups[f_peaksRBM[ni]].append(p)
                    groups_Int[f_peaksRBM[ni]].append(PeaksInt2[jj])
                except ValueError:
                    continue
            
            for n in groups:
                meanRBM[labels[z]][n]=np.mean(groups[n])
                stdRBM[labels[z]][n]=np.std(groups[n])
                meanRBM_Int[labels[z]][n]=np.mean(groups_Int[n])
                stdRBM_Int[labels[z]][n]=np.std(groups_Int[n])
            
            lb=[str(round(meanRBM[labels[z]][n],2))+ ' $\pm$ ' +str(round(stdRBM[labels[z]][n],2)) for n in stdRBM[labels[z]]]
                
            ax_RBM.hist(PeaksLoc,binsShift,color=clr[round(len(clr)/2)],label=labels[z]+': '+', '.join(map(str, lb))+' $cm^{-1}$',
                          alpha=0.5,ec='k',align='left',density=dens)
            
            ax_RBM.legend(fontsize=fs-2)  
            ax_RBM.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
            ax_RBM.set_ylabel('counts',fontsize=fs)
            ax_RBM.tick_params(axis="both", labelsize=fs)   
            ax_RBM.set_title('Raman shift $RBM$ modes',fontsize=fs+2);
            ax_RBM.xaxis.set_major_locator(MaxNLocator(integer=True))



#%% #Correlations

        # Position Band 1 (G) vs. Intensity Ratio (Band 2/Band 1)
        if correlation1 == 1:
            if z==0:
                fig_Corr1,ax_Corr1= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_Corr1.canvas.manager.set_window_title(band1Name +' band vs I('+ band2Name +')/I_' +band1Name+ ')')
    
            if nt==1:
                # In the case where there are separate G+/G- peaks, use the most intense one
                mask=Int_1plus==maxInt_1
                maxCenter1=np.where(mask,center_1plus,center_1min)
                
                idx = np.isfinite(I21) & np.isfinite(maxCenter1)
                #Linear fit Intensity Ratio vs G+
                fit_I21_1=np.polyfit(np.array(I21)[idx],np.array(maxCenter1)[idx],1)
                #Plot Intensity Ratio vs G+
                ax_Corr1.plot(I21,maxCenter1,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_I21_1[0],3)))         
                ax_Corr1.set_title('$'+band1Name +'$ band Shift vs. $I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr1.set_xlabel('$I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr1.set_ylabel('$'+band1Name +'$ band Shift / $cm^{-1}$',fontsize=fs) 
            else:
                idx = np.isfinite(I21) & np.isfinite(center_1plus)
                #Linear fit Intensity Ratio vs G+
                fit_I21_1=np.polyfit(np.array(I21)[idx],np.array(center_1plus)[idx],1)
                #Plot Intensity Ratio vs G+
                ax_Corr1.plot(I21,center_1plus,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_I21_1[0],3)))         
                ax_Corr1.set_title('$'+band1Name +'$ band Shift vs. $I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr1.set_xlabel('$I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr1.set_ylabel('$'+band1Name +'$ band Shift / $cm^{-1}$',fontsize=fs) 
                    
            ax_Corr1.legend(fontsize=fs-3)
            ax_Corr1.tick_params(axis="both", labelsize=fs)
            ax_Corr1.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            xx=np.linspace(min(I21),max(I21),100)
            ax_Corr1.plot(xx,fit_I21_1[0]*xx+fit_I21_1[1],lw=2,ls='-',color=clr[round(len(clr)/3)])

# Position Band 2 (D) vs. Intensity Ratio (Band 2/Band 1)
        if correlation2 == 1:
                if z==0:
                    fig_Corr2,ax_Corr2= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                    fig_Corr2.canvas.manager.set_window_title(band2Name +' band vs I('+ band2Name +')/I_' +band1Name+ ')')
                      
                #Linear fit 2D vs Gplus
                idx = np.isfinite(I21) & np.isfinite(center_2)
                fit_I21_2=np.polyfit(np.array(I21)[idx],np.array(center_2)[idx],1)
                #Plot 2D vs G
                ax_Corr2.plot(I21,center_2,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_I21_2[0],3)))         
                ax_Corr2.set_title('$'+band2Name +'$ band Shift vs. $I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr2.set_xlabel('$I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr2.set_ylabel('$'+band2Name +'$ band Shift / $cm^{-1}$',fontsize=fs)         
                ax_Corr2.legend(fontsize=fs-3)
                ax_Corr2.tick_params(axis="both", labelsize=fs)
                ax_Corr2.yaxis.set_major_locator(MaxNLocator(integer=True))
                
                #Plot Trendline
                xx=np.linspace(min(I21),max(I21),100)
                ax_Corr2.plot(xx,fit_I21_2[0]*xx+fit_I21_2[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
               
        # Position Band 3 (2D) vs. Intensity Ratio (Band 2/Band 1)
        if correlation3 == 1:
                if z==0:
                    fig_Corr3,ax_Corr3= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                    fig_Corr3.canvas.manager.set_window_title(band3Name +' band vs I('+ band2Name +')/I_' +band1Name+ ')')
                      
                #Linear fit 2D vs Gplus
                idx = np.isfinite(I21) & np.isfinite(center_3)
                fit_I21_3=np.polyfit(np.array(I21)[idx],np.array(center_3)[idx],1)
                #Plot 2D vs G
                ax_Corr3.plot(I21,center_3,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_I21_3[0],3)))         
                ax_Corr3.set_title('$'+band3Name +'$ band Shift vs. $I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr3.set_xlabel('$I_{'+band2Name +'}/I_{'+band1Name+'}$',fontsize=fs)
                ax_Corr3.set_ylabel('$'+band3Name +'$ band Shift / $cm^{-1}$',fontsize=fs)         
                ax_Corr3.legend(fontsize=fs-3)
                ax_Corr3.tick_params(axis="both", labelsize=fs)
                ax_Corr3.yaxis.set_major_locator(MaxNLocator(integer=True))
                
                
                #Plot Trendline
                xx=np.linspace(min(I21),max(I21),100)
                ax_Corr3.plot(xx,fit_I21_3[0]*xx+fit_I21_3[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
                
                
        # Position Band 3 (2D) vs. Position Band 1 (G)
        if correlation4 == 1:
            if z==0:
                fig_Corr4,ax_Corr4= plt.subplots(1,1, figsize=(8,6), constrained_layout=True)
                fig_Corr4.canvas.manager.set_window_title(band3Name +' band vs '+ band1Name +' band')
    
            if nt==1:
                # In the case where there are separate G+/G- peaks, use the most intense one
                mask=Int_1plus==maxInt_1
                maxCenter1=np.where(mask,center_1plus,center_1min)
                idx = np.isfinite(center_3) & np.isfinite(maxCenter1)
                #Linear fit Intensity Ratio vs G+
                fit_1_3=np.polyfit(np.array(maxCenter1)[idx],np.array(center_3)[idx],1)
                #Plot peak 3 + vs. peak 1 +
                ax_Corr4.plot(maxCenter1, center_3,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_1_3[0],3)))         
                ax_Corr4.set_title(band3Name+' band Shift vs. $'+band1Name +'$ band Shift',fontsize=fs)
                ax_Corr4.set_ylabel('$'+band3Name +'$ band Shift / $cm^{-1}$',fontsize=fs)
                ax_Corr4.set_xlabel('$'+band1Name +'$ band Shift / $cm^{-1}$',fontsize=fs)
                xx=np.linspace(min(center_1plus),max(center_1plus),100)

            else:
                idx = np.isfinite(center_3) & np.isfinite(center_1plus)
                #Linear fit Intensity Ratio vs G
                fit_1_3=np.polyfit(np.array(center_1plus)[idx],np.array(center_3)[idx],1)
                #Plot peak 3 + vs. peak 1 
                ax_Corr4.plot(center_1plus, center_3,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_1_3[0],3)))         
                ax_Corr4.set_title(band3Name+' band Shift vs. $'+band1Name +'$ band Shift',fontsize=fs)
                ax_Corr4.set_ylabel('$'+band3Name +'$ band Shift / $cm^{-1}$',fontsize=fs)
                ax_Corr4.set_xlabel('$'+band1Name +'$ band Shift / $cm^{-1}$',fontsize=fs)
                xx=np.linspace(min(center_1plus),max(center_1plus),100)
                    
            ax_Corr4.legend(fontsize=fs-3)
            ax_Corr4.tick_params(axis="both", labelsize=fs)
            ax_Corr4.yaxis.set_major_locator(MaxNLocator(integer=True))
            
            
            ax_Corr4.plot(xx,fit_1_3[0]*xx+fit_1_3[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
            
        #  FWHM vs. Position Bands
        if correlation5 == 1 and lorentz == 1:
            if z==0:
                fig_Corr5,ax_Corr5= plt.subplots(1,3, figsize=(15,5), constrained_layout=True)
                fig_Corr5.canvas.manager.set_window_title('FWHM vs. Position of Each Band')
            
            # Peak 1 FWHM vs. Position
            if nt==1:
                # In the case where there are separate G+/G- peaks, use the most intense one
                mask=Int_1plus==maxInt_1
                maxCenter1=np.where(mask,center_1plus,center_1min)
                maxFWHM1=np.where(mask,FWHM_1plus,FWHM_1min)

                idx = np.isfinite(maxCenter1) & np.isfinite(maxFWHM1)
                #Linear fit G+ FWHM vs. Position
                fit_1_FWHM=np.polyfit(maxCenter1[idx],maxFWHM1[idx],1)
                #Plot G+ FWHM vs. Position
                ax_Corr5[0].plot(maxCenter1, maxFWHM1,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_1_FWHM[0],3)))         
                ax_Corr5[0].set_title('$'+band1Name+'$ band FWHM vs. Shift',fontsize=fs)
                ax_Corr5[0].set_ylabel('$'+band1Name +'$ band Shift',fontsize=fs)
                ax_Corr5[0].set_xlabel('$'+band1Name +'$ band FWHM / $cm^{-1}$',fontsize=fs)
                xx=np.linspace(min(maxCenter1),max(maxCenter1),100)

            else:
                idx = np.isfinite(center_1plus) & np.isfinite(FWHM_1plus)
                #Linear fit G+ FWHM vs. Position
                fit_1_FWHM=np.polyfit(np.array(center_1plus)[idx],np.array(FWHM_1plus)[idx],1)
                #Plot G+ FWHM vs. Position
                ax_Corr5[0].plot(center_1plus, FWHM_1plus,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_1_FWHM[0],3)))         
                ax_Corr5[0].set_title('$'+band1Name+'$ band FWHM vs. Shift',fontsize=fs)
                ax_Corr5[0].set_ylabel('$'+band1Name +'$ band Shift',fontsize=fs)
                ax_Corr5[0].set_xlabel('$'+band1Name +'$ band FWHM / $cm^{-1}$',fontsize=fs)
                xx=np.linspace(min(center_1plus),max(center_1plus),100)
                    
                    
            #Plot Trendline
            ax_Corr5[0].plot(xx,fit_1_FWHM[0]*xx+fit_1_FWHM[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
            
            # Peak 2 FWHM vs. Position
            idx = np.isfinite(center_2) & np.isfinite(FWHM_2)
            #Linear fit G+ FWHM vs. Position
            fit_2_FWHM=np.polyfit(np.array(center_2)[idx],np.array(FWHM_2)[idx],1)
            #Plot G+ FWHM vs. Position
            ax_Corr5[1].plot(center_2, FWHM_2,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_2_FWHM[0],3)))         
            ax_Corr5[1].set_title('$'+band2Name+'$ band FWHM vs. Shift',fontsize=fs)
            ax_Corr5[1].set_ylabel('$'+band2Name +'$ band Shift',fontsize=fs)
            ax_Corr5[1].set_xlabel('$'+band2Name +'$ band FWHM / $cm^{-1}$',fontsize=fs)
            #Plot Trendline
            xx=np.linspace(min(center_2),max(center_2),100)       
            ax_Corr5[1].plot(xx,fit_2_FWHM[0]*xx+fit_2_FWHM[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
       
            # Peak 3 FWHM vs. Position
            idx = np.isfinite(center_3) & np.isfinite(FWHM_3)
            #Linear fit G+ FWHM vs. Position
            fit_3_FWHM=np.polyfit(np.array(center_3)[idx],np.array(FWHM_3)[idx],1)
            #Plot G+ FWHM vs. Position
            ax_Corr5[2].plot(center_3, FWHM_3,marker='o',ls='',color=clr[round(len(clr)/2)],label=labels[z]+', slope='+str(round(fit_3_FWHM[0],3)))         
            ax_Corr5[2].set_title('$'+band3Name+'$ band FWHM vs. Shift',fontsize=fs)
            ax_Corr5[2].set_ylabel('$'+band3Name +'$ band Shift',fontsize=fs)
            ax_Corr5[2].set_xlabel('$'+band3Name +'$ band FWHM / $cm^{-1}$',fontsize=fs)
            #Plot Trendline
            xx=np.linspace(min(center_3),max(center_3),100)       
            ax_Corr5[2].plot(xx,fit_3_FWHM[0]*xx+fit_3_FWHM[1],lw=2,ls='-',color=clr[round(len(clr)/3)])
              
            fig_Corr5.suptitle('FWHM vs. Position of Each Band',fontsize=fs+3)
            for axx in ax_Corr5:     
                axx.legend(fontsize=fs-3)
                axx.tick_params(axis="both", labelsize=fs)
                axx.yaxis.set_major_locator(MaxNLocator(integer=True))
                
#%% #2D Maps: Plot 2D Maps if selected
        if map_==1:
            #Change folder to save 2D maps
            os.chdir(folder_selected + '/'+nameFigs+'_Results/')

            if use_leng!=1:
                ratio=cols/rows
            else:
                ratio=xlen/ylen
                
            if use_leng!=1:
                y1=np.linspace(0,rows-1,rows)
                x1=np.linspace(0,cols-1,cols)
                ylab='y pixel'
                xlab='x pixel'
            else:
                y1=np.linspace(0,ylen,rows)
                x1=np.linspace(0,xlen,cols)
                ylab='y length, ($\mu m$)'
                xlab='x length, ($\mu m$)'    

            xx,yy=np.meshgrid(x1,y1)
            xx=xx.flatten()
            yy=yy.flatten()
            
            
            #%%Intensity Ratio peak 2 - peak 1 2D Map
            
            try:
                fig_MapI21.clf()
            except NameError:
                pass
                fig_MapI21 = plt.figure('Intensity Ratio '+band2Name+'-'+band1Name+' Map', figsize=(7*ratio,6), constrained_layout=True)
                ax_MapI21 = fig_MapI21.add_subplot(111)
                
            var2D_I21=I21
            spec_labelI21='$I_{'+band2Name+'}/I_{'+band1Name+'}$'
            unitI21=''
            var2D_I21_r=var2D_I21.reshape([rows, cols])
           
            #Make 2D Map of Intensity Data            
            h=ax_MapI21.pcolormesh(x1,y1,var2D_I21_r)
            clbr=plt.colorbar(h,ax=ax_MapI21)
            clbr.set_label(spec_labelI21, rotation=270,labelpad=20,fontsize=fs)
            clbr.ax.tick_params(labelsize=fs) 
            ax_MapI21.set_xlabel(xlab,fontsize=fs)
            ax_MapI21.set_ylabel(ylab,fontsize=fs)
            ax_MapI21.set_title('Map of '+spec_labelI21+' for '+labels[z],fontsize=fs)        
            ax_MapI21.tick_params(axis="both", labelsize=fs)            
            fig_MapI21.savefig(labels[z]+'_2DMap_I21'+imgtype)
               
        
            pkn=5
            lineI21, = ax_MapI21.plot(xx, yy, 'b',picker=pkn)
            lineI21.set_visible(False) 
            
            #%%peak 1 position 2D Map
            
            try:
                fig_Map_peak1.clf()
            except NameError:
                pass
                fig_Map_peak1 = plt.figure(band1Name+' shift Map', figsize=(7*ratio,6), constrained_layout=True)
                ax_Map_peak1 = fig_Map_peak1.add_subplot(111)
                
            if nt==1:
                mask=Int_1plus==maxInt_1
                maxCenter1=np.where(mask,center_1plus,center_1min)
                var2D_peak1=maxCenter1
            else:
                var2D_peak1=center_1plus

                
            spec_label_peak1='$'+band1Name+'$ shift'
                
            unit_peak1=' $cm^{-1}$'        
            var2D_peak1_r=var2D_peak1.reshape([rows, cols])
           
            #Make 2D Map of Intensity Data            
            h=ax_Map_peak1.pcolormesh(x1,y1,var2D_peak1_r)
            clbr=plt.colorbar(h,ax=ax_Map_peak1)
            clbr.set_label(spec_label_peak1, rotation=270,labelpad=20,fontsize=fs)
            clbr.ax.tick_params(labelsize=fs) 
            ax_Map_peak1.set_xlabel(xlab,fontsize=fs)
            ax_Map_peak1.set_ylabel(ylab,fontsize=fs)
            ax_Map_peak1.set_title('Map of '+spec_label_peak1+' for '+labels[z],fontsize=fs)        
            ax_Map_peak1.tick_params(axis="both", labelsize=fs)            
            fig_Map_peak1.savefig(labels[z]+'_2DMap_'+band1Name+imgtype)
               
        
            pkn=5
            line_peak1, = ax_Map_peak1.plot(xx, yy, 'b',picker=pkn)
            line_peak1.set_visible(False) 
            
            #%%peak 2 position 2D Map
            
            try:
                fig_Map_peak2.clf()
            except NameError:
                pass
                fig_Map_peak2 = plt.figure(band2Name+' shift Map', figsize=(7*ratio,6), constrained_layout=True)
                ax_Map_peak2 = fig_Map_peak2.add_subplot(111)
                
            var2D_peak2=center_2
                
            spec_label_peak2='$'+band2Name+'$ shift'
                
            unit_peak2=' $cm^{-1}$'        
            var2D_peak2_r=var2D_peak2.reshape([rows, cols])
           
            #Make 2D Map of Intensity Data            
            h=ax_Map_peak2.pcolormesh(x1,y1,var2D_peak2_r)
            clbr=plt.colorbar(h,ax=ax_Map_peak2)
            clbr.set_label(spec_label_peak2, rotation=270,labelpad=20,fontsize=fs)
            clbr.ax.tick_params(labelsize=fs) 
            ax_Map_peak2.set_xlabel(xlab,fontsize=fs)
            ax_Map_peak2.set_ylabel(ylab,fontsize=fs)
            ax_Map_peak2.set_title('Map of '+spec_label_peak2+' for '+labels[z],fontsize=fs)        
            ax_Map_peak2.tick_params(axis="both", labelsize=fs)            
            fig_peak2.savefig(labels[z]+'_2DMap_'+band2Name+imgtype)
               
        
            pkn=5
            line_peak2, = ax_Map_peak2.plot(xx, yy, 'b',picker=pkn)
            line_peak2.set_visible(False) 
            
            #%%peak 3 position 2D Map
            
            try:
                fig_Map_peak3.clf()
            except NameError:
                pass
                fig_Map_peak3 = plt.figure(band3Name+' shift Map', figsize=(7*ratio,6), constrained_layout=True)
                ax_Map_peak3 = fig_Map_peak3.add_subplot(111)
                
            var2D_peak3=center_3
                
            spec_label_peak3='$'+band3Name+'$ shift'
                
            unit_peak3=' $cm^{-1}$'        
            var2D_peak3_r=var2D_peak3.reshape([rows, cols])
           
            #Make 2D Map of Intensity Data            
            h=ax_Map_peak3.pcolormesh(x1,y1,var2D_peak3_r)
            clbr=plt.colorbar(h,ax=ax_Map_peak3)
            clbr.set_label(spec_label_peak3, rotation=270,labelpad=20,fontsize=fs)
            clbr.ax.tick_params(labelsize=fs) 
            ax_Map_peak3.set_xlabel(xlab,fontsize=fs)
            ax_Map_peak3.set_ylabel(ylab,fontsize=fs)
            ax_Map_peak3.set_title('Map of '+spec_label_peak3+' for '+labels[z],fontsize=fs)        
            ax_Map_peak3.tick_params(axis="both", labelsize=fs)            
            fig_peak3.savefig(labels[z]+'_2DMap_'+band3Name+imgtype)
               
        
            pkn=5
            line_peak3, = ax_Map_peak3.plot(xx, yy, 'b',picker=pkn)
            line_peak3.set_visible(False) 
            #%%peak 1 FWHM 2D Map            
            if map1==1:
                try:
                    fig_Map1.clf()
                except NameError:
                    pass
                    fig_Map1 = plt.figure('FWHM '+band1Name+' band Map', figsize=(7*ratio,6), constrained_layout=True)
                    ax_Map1 = fig_Map1.add_subplot(111)
                    
                
                if nt==1:
                    mask=Int_1plus==maxInt_1
                    maxFWHM1=np.where(mask,FWHM_1plus,FWHM_1min)
                    var2D_1=maxFWHM1

                else:
                    var2D_1=FWHM_1plus

                spec_label1='$'+band1Name+'$ FWHM'
                unit1=' $cm^{-1}$'
                
                var2D_1_r=var2D_1.reshape([rows, cols])
                
               
                #Make 2D Map of Intensity Data            
                h=ax_Map1.pcolormesh(x1,y1,var2D_1_r)
                clbr=plt.colorbar(h,ax=ax_Map1)
                clbr.set_label(spec_label1, rotation=270,labelpad=20,fontsize=fs)
                clbr.ax.tick_params(labelsize=fs) 
                ax_Map1.set_xlabel(xlab,fontsize=fs)
                ax_Map1.set_ylabel(ylab,fontsize=fs)
                ax_Map1.set_title('Map of '+spec_label1+' for '+labels[z],fontsize=fs)        
                ax_Map1.tick_params(axis="both", labelsize=fs)            
                fig_Map1.savefig(labels[z]+'_2DMap_'+band1Name+'_FWHM'+imgtype)
                   
            
                pkn=5
                line1, = ax_Map1.plot(xx, yy, 'b',picker=pkn)
                line1.set_visible(False) 

            #%%peak 2 FWHM 2D Map            
            if map2==1:
                try:
                    fig_Map2.clf()
                except NameError:
                    pass
                    fig_Map2 = plt.figure('FWHM '+band2Name+' band Map', figsize=(7*ratio,6), constrained_layout=True)
                    ax_Map2 = fig_Map2.add_subplot(111)
                    
                var2D_2=FWHM_2
                

                spec_label2='$'+band2Name+'$ FWHM'
                unit2=' $cm^{-1}$'
                
                var2D_2_r=var2D_2.reshape([rows, cols])
                
               
                #Make 2D Map of Intensity Data            
                h=ax_Map2.pcolormesh(x1,y1,var2D_2_r)
                clbr=plt.colorbar(h,ax=ax_Map2)
                clbr.set_label(spec_label1, rotation=270,labelpad=20,fontsize=fs)
                clbr.ax.tick_params(labelsize=fs) 
                ax_Map2.set_xlabel(xlab,fontsize=fs)
                ax_Map2.set_ylabel(ylab,fontsize=fs)
                ax_Map2.set_title('Map of '+spec_label2+' for '+labels[z],fontsize=fs)        
                ax_Map2.tick_params(axis="both", labelsize=fs)            
                fig_Map2.savefig(labels[z]+'_2DMap_'+band2Name+'_FWHM'+imgtype)
                   
            
                pkn=5
                line2, = ax_Map2.plot(xx, yy, 'b',picker=pkn)
                line2.set_visible(False)     
                
            #%%peak 3 FWHM 2D Map            
            if map3==1:
                try:
                    fig_Map3.clf()
                except NameError:
                    pass
                    fig_Map3 = plt.figure('FWHM '+band3Name+' band Map', figsize=(7*ratio,6), constrained_layout=True)
                    ax_Map3 = fig_Map3.add_subplot(111)
                    
                var2D_3=FWHM_3
                

                spec_label3='$'+band3Name+'$ FWHM'
                unit3=' $cm^{-1}$'
                
                var2D_3_r=var2D_3.reshape([rows, cols])
                
               
                #Make 2D Map of Intensity Data            
                h=ax_Map3.pcolormesh(x1,y1,var2D_3_r)
                clbr=plt.colorbar(h,ax=ax_Map3)
                clbr.set_label(spec_label3, rotation=270,labelpad=20,fontsize=fs)
                clbr.ax.tick_params(labelsize=fs) 
                ax_Map3.set_xlabel(xlab,fontsize=fs)
                ax_Map3.set_ylabel(ylab,fontsize=fs)
                ax_Map3.set_title('Map of '+spec_label3+' for '+labels[z],fontsize=fs)        
                ax_Map3.tick_params(axis="both", labelsize=fs)            
                fig_Map3.savefig(labels[z]+'_2DMap_'+band2Name+'_FWHM'+imgtype)
                   
            
                pkn=5
                line3, = ax_Map3.plot(xx, yy, 'b',picker=pkn)
                line3.set_visible(False)       
                
            if map4==1:
                try:
                    fig_Map4.clf()
                except NameError:
                    pass
                    fig_Map4 = plt.figure('Intensity Ratio '+band3Name+'-'+band1Name+' Map', figsize=(7*ratio,6), constrained_layout=True)
                    ax_Map4 = fig_Map4.add_subplot(111)
                    
                var2D_4=I31
                spec_label4='$I_{'+band3Name+'}/I_{'+band1Name+'}$'
                unit4=''
                var2D_4_r=var2D_4.reshape([rows, cols])
               
                #Make 2D Map of Intensity Data            
                h=ax_Map4.pcolormesh(x1,y1,var2D_4_r)
                clbr=plt.colorbar(h,ax=ax_Map4)
                clbr.set_label(spec_label4, rotation=270,labelpad=20,fontsize=fs)
                clbr.ax.tick_params(labelsize=fs) 
                ax_Map4.set_xlabel(xlab,fontsize=fs)
                ax_Map4.set_ylabel(ylab,fontsize=fs)
                ax_Map4.set_title('Map of '+spec_label4+' for '+labels[z],fontsize=fs)        
                ax_Map4.tick_params(axis="both", labelsize=fs)            
                fig_Map4.savefig(labels[z]+'_2DMap_I31'+imgtype)
                   
            
                pkn=5
                line4, = ax_Map4.plot(xx, yy, 'b',picker=pkn)
                line4.set_visible(False) 
                
            
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%%%%%% Output data %%%%%%%%%%%%%%%
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Save data in txt file                             
        if rbm == 1:
            if lorentz == 0:
                T = pd.DataFrame({
                    'Intensity_' + band1Name: Int_1plus,
                    'Shift_' + band1Name: center_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3,
                    'Intensity_RBM': PeaksInt,
                    'Shift_RBM': PeaksLoc1
                })
                
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'_avg', 'Intensity_' + band1Name+'_std',
                                                 'Shift_' + band1Name+'_avg', 'Shift_' + band1Name+'_std',                                       
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std',
                                                 'Intensity_RBM_avg','Intensity_RBM_std',
                                                 'Shift_RBM_avg','Shift_RBM_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'_avg': peak1plus_av,
                    'Shift_' + band1Name+'_std': peak1plus_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error,                    
                    'Intensity_RBM_avg': [meanRBM_Int[labels[z]][n] for n in meanRBM_Int[labels[z]]],
                    'Intensity_RBM_std': [stdRBM_Int[labels[z]][n] for n in stdRBM_Int[labels[z]]],
                    'Shift_RBM_avg': [meanRBM[labels[z]][n] for n in meanRBM[labels[z]]],
                    'Shift_RBM_std': [stdRBM[labels[z]][n] for n in stdRBM[labels[z]]]
                    }
            elif lorentz == 1 and nt == 1:
                T = pd.DataFrame({
                    'Intensity_' + band1Name + '-': Int_1min,
                    'Shift_' + band1Name + '-': center_1min,
                    'FWHM_' + band1Name + '-': FWHM_1min,
                    'Intensity_' + band1Name + '+': Int_1plus,
                    'Shift_' + band1Name + '+': center_1plus,
                    'FWHM_' + band1Name + '+': FWHM_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'FWHM_' + band2Name: FWHM_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3,
                    'FWHM_' + band3Name: FWHM_3,
                    'Intensity_RBM': PeaksInt,
                    'Shift_RBM': PeaksLoc1
                })
                
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'+_avg', 'Intensity_' + band1Name+'+_std',
                                                 'Shift_' + band1Name+'+_avg', 'Shift_' + band1Name+'+_std',
                                                 'FWHM_' + band1Name+'+_avg', 'FWHM_' + band1Name+'+_std',
                                                 'Intensity_' + band1Name+'-_avg', 'Intensity_' + band1Name+'-_std',
                                                 'Shift_' + band1Name+'-_avg', 'Shift_' + band1Name+'-_std',
                                                 'FWHM_' + band1Name+'-_avg', 'FWHM_' + band1Name+'-_std'
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'FWHM_' + band2Name+'_avg', 'FWHM_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'FWHM_' + band3Name+'_avg', 'FWHM_' + band3Name+'_std',
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std',
                                                 'Intensity_RBM_avg','Intensity_RBM_std',
                                                 'Shift_RBM_avg','Shift_RBM_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'+_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'+_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'+_avg': peak1plus_av,
                    'Shift_' + band1Name+'+_std': peak1plus_std,
                    'FWHM_' + band1Name+'+_avg': FWHM_1plus_av,
                    'FWHM_' + band1Name+'+_std': FWHM_1plus_std,
                    'Intensity_' + band1Name+'-_avg': np.nanmean(Int_1min),
                    'Intensity_' + band1Name+'-_std':np.nanstd(Int_1min),
                    'Shift_' + band1Name+'-_avg': peak1min_av,
                    'Shift_' + band1Name+'-_std': peak1min_std,
                    'FWHM_' + band1Name+'-_avg': FWHM_1min_av,
                    'FWHM_' + band1Name+'-_std': FWHM_1min_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'FWHM_' + band2Name+'_avg': FWHM_2_av,
                    'FWHM_' + band2Name+'_std': FWHM_2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'FWHM_' + band3Name+'_avg': FWHM_3_av,
                    'FWHM_' + band3Name+'_std': FWHM_3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error, 
                    'Intensity_RBM_avg': [meanRBM_Int[labels[z]][n] for n in meanRBM_Int[labels[z]]],
                    'Intensity_RBM_std': [stdRBM_Int[labels[z]][n] for n in stdRBM_Int[labels[z]]],
                    'Shift_RBM_avg': [meanRBM[labels[z]][n] for n in meanRBM[labels[z]]],
                    'Shift_RBM_std': [stdRBM[labels[z]][n] for n in stdRBM[labels[z]]]
                    }
            elif lorentz == 1 and nt == 0:
                T = pd.DataFrame({
                    'Intensity_' + band1Name: Int_1plus,
                    'Shift_' + band1Name: center_1plus,
                    'FWHM_' + band1Name: FWHM_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'FWHM_' + band2Name: FWHM_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3,
                    'FWHM_' + band3Name: FWHM_3,
                    'Intensity_RBM': PeaksInt,
                    'Shift_RBM': PeaksLoc1
                })
                
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'_avg', 'Intensity_' + band1Name+'_std',
                                                 'Shift_' + band1Name+'_avg', 'Shift_' + band1Name+'_std',
                                                 'FWHM_' + band1Name+'_avg', 'FWHM_' + band1Name+'_std',
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'FWHM_' + band2Name+'_avg', 'FWHM_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'FWHM_' + band3Name+'_avg', 'FWHM_' + band3Name+'_std',                                                 
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std',
                                                 'Intensity_RBM_avg','Intensity_RBM_std',
                                                 'Shift_RBM_avg','Shift_RBM_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'_avg': peak1plus_av,
                    'Shift_' + band1Name+'_std': peak1plus_std,
                    'FWHM_' + band1Name+'_avg': FWHM_1plus_av,
                    'FWHM_' + band1Name+'_std': FWHM_1plus_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'FWHM_' + band2Name+'_avg': FWHM_2_av,
                    'FWHM_' + band2Name+'_std': FWHM_2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'FWHM_' + band3Name+'_avg': FWHM_3_av,
                    'FWHM_' + band3Name+'_std': FWHM_3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error, 
                    'Intensity_RBM_avg': [meanRBM_Int[labels[z]][n] for n in meanRBM_Int[labels[z]]],
                    'Intensity_RBM_std': [stdRBM_Int[labels[z]][n] for n in stdRBM_Int[labels[z]]],
                    'Shift_RBM_avg': [meanRBM[labels[z]][n] for n in meanRBM[labels[z]]],
                    'Shift_RBM_std': [stdRBM[labels[z]][n] for n in stdRBM[labels[z]]]
                    }
                
                
        else:
            if lorentz == 0:
                T = pd.DataFrame({
                    'Intensity_' + band1Name: Int_1plus,
                    'Shift_' + band1Name: center_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3
                })
                
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'_avg', 'Intensity_' + band1Name+'_std',
                                                 'Shift_' + band1Name+'_avg', 'Shift_' + band1Name+'_std',                                       
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'_avg': peak1plus_av,
                    'Shift_' + band1Name+'_std': peak1plus_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error 
                    }
            elif lorentz == 1 and nt == 1:
                T = pd.DataFrame({
                    'Intensity_' + band1Name + '-': Int_1min,
                    'Shift_' + band1Name + '-': center_1min,
                    'FWHM_' + band1Name + '-': FWHM_1min,
                    'Intensity_' + band1Name + '+': Int_1plus,
                    'Shift_' + band1Name + '+': center_1plus,
                    'FWHM_' + band1Name + '+': FWHM_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'FWHM_' + band2Name: FWHM_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3,
                    'FWHM_' + band3Name: FWHM_3
                })
                
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'+_avg', 'Intensity_' + band1Name+'+_std',
                                                 'Shift_' + band1Name+'+_avg', 'Shift_' + band1Name+'+_std',
                                                 'FWHM_' + band1Name+'+_avg', 'FWHM_' + band1Name+'+_std',
                                                 'Intensity_' + band1Name+'-_avg', 'Intensity_' + band1Name+'-_std',
                                                 'Shift_' + band1Name+'-_avg', 'Shift_' + band1Name+'-_std',
                                                 'FWHM_' + band1Name+'-_avg', 'FWHM_' + band1Name+'-_std'
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'FWHM_' + band2Name+'_avg', 'FWHM_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'FWHM_' + band3Name+'_avg', 'FWHM_' + band3Name+'_std',
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'+_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'+_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'+_avg': peak1plus_av,
                    'Shift_' + band1Name+'+_std': peak1plus_std,
                    'FWHM_' + band1Name+'+_avg': FWHM_1plus_av,
                    'FWHM_' + band1Name+'+_std': FWHM_1plus_std,
                    'Intensity_' + band1Name+'-_avg': np.nanmean(Int_1min),
                    'Intensity_' + band1Name+'-_std':np.nanstd(Int_1min),
                    'Shift_' + band1Name+'-_avg': peak1min_av,
                    'Shift_' + band1Name+'-_std': peak1min_std,
                    'FWHM_' + band1Name+'-_avg': FWHM_1min_av,
                    'FWHM_' + band1Name+'-_std': FWHM_1min_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'FWHM_' + band2Name+'_avg': FWHM_2_av,
                    'FWHM_' + band2Name+'_std': FWHM_2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'FWHM_' + band3Name+'_avg': FWHM_3_av,
                    'FWHM_' + band3Name+'_std': FWHM_3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error 
                    }
                
            elif lorentz == 1 and nt == 0:
                T = pd.DataFrame({
                    'Intensity_' + band1Name: Int_1plus,
                    'Shift_' + band1Name: center_1plus,
                    'FWHM_' + band1Name: FWHM_1plus,
                    'Intensity_' + band2Name: Int_2,
                    'Shift_' + band2Name: center_2,
                    'FWHM_' + band2Name: FWHM_2,
                    'Intensity_' + band3Name: Int_3,
                    'Shift_' + band3Name: center_3,
                    'FWHM_' + band3Name: FWHM_3
                })
                if z==0:
                    avg_df=pd.DataFrame(columns=['file', 'Intensity_' + band1Name+'_avg', 'Intensity_' + band1Name+'_std',
                                                 'Shift_' + band1Name+'_avg', 'Shift_' + band1Name+'_std',
                                                 'FWHM_' + band1Name+'_avg', 'FWHM_' + band1Name+'_std',
                                                 'Intensity_' + band2Name+'_avg', 'Intensity_' + band2Name+'_std',
                                                 'Shift_' + band2Name+'_avg', 'Shift_' + band2Name+'_std',
                                                 'FWHM_' + band2Name+'_avg', 'FWHM_' + band2Name+'_std',
                                                 'Intensity_' + band3Name+'_avg', 'Intensity_' + band3Name+'_std',
                                                 'Shift_' + band3Name+'_avg', 'Shift_' + band3Name+'_std',
                                                 'FWHM_' + band3Name+'_avg', 'FWHM_' + band3Name+'_std',
                                                 'IntensityRatio_' + band2Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band2Name+'/'+band1Name+'_std',
                                                 'IntensityRatio_' + band3Name+'/'+band1Name+'_avg', 'IntensityRatio_' + band3Name+'/'+band1Name+'_std'])
                avg_df.loc[z]={
                    'file':file_name[z],
                    'Intensity_' + band1Name+'_avg': np.nanmean(Int_1plus),
                    'Intensity_' + band1Name+'_std':np.nanstd(Int_1plus),
                    'Shift_' + band1Name+'_avg': peak1plus_av,
                    'Shift_' + band1Name+'_std': peak1plus_std,
                    'FWHM_' + band1Name+'_avg': FWHM_1plus_av,
                    'FWHM_' + band1Name+'_std': FWHM_1plus_std,
                    'Intensity_' + band2Name+'_avg': np.nanmean(Int_2),
                    'Intensity_' + band2Name+'_std':np.nanstd(Int_2),
                    'Shift_' + band2Name+'_avg': peak2_av,
                    'Shift_' + band2Name+'_std': peak2_std,
                    'FWHM_' + band2Name+'_avg': FWHM_2_av,
                    'FWHM_' + band2Name+'_std': FWHM_2_std,
                    'Intensity_' + band3Name+'_avg': np.nanmean(Int_3),
                    'Intensity_' + band3Name+'_std':np.nanstd(Int_3),
                    'Shift_' + band3Name+'_avg': peak3_av,
                    'Shift_' + band3Name+'_std': peak3_std,
                    'FWHM_' + band3Name+'_avg': FWHM_3_av,
                    'FWHM_' + band3Name+'_std': FWHM_3_std,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_avg':I21_av,
                    'IntensityRatio_' + band2Name+'/'+band1Name+'_std':I21_error,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_avg':I31_av,
                    'IntensityRatio_' + band3Name+'/'+band1Name+'_std':I31_error 
                    }
        
        T['Intensity Ratio_2_1']=I21  
        T['Intensity Ratio_3_1']=I31
        T.index.name='Spectra #'
        
        nameT = folder_selected+'/'+nameFigs+'_Results/'+labels[z] + '_results.csv'
        T.to_csv(nameT)
        if z==total-1:
            name_avg = folder_selected+'/'+nameFigs+'_Results/'+nameFigs + '_Average_results.csv'
            avg_df.to_csv(name_avg,index=False)

                    
            
 #%%% Save Figures in Scaleable Vector Format
    if saveFigs==1:
        os.chdir(folder_selected + '/'+nameFigs+'_Results/')
        fig_avg.savefig(nameFigs+'_AverageSpectra'+imgtype)
        fig_IR.savefig(nameFigs+'_I' + band2Name + band1Name + imgtype)
        fig_IR2.savefig(nameFigs+'_I' + band3Name + band1Name + imgtype)
        fig_peak3.savefig(nameFigs+'_Shift_'+band3Name + imgtype)
        fig_peak1.savefig(nameFigs+'_Shift_'+band1Name + imgtype)
        fig_peak2.savefig(nameFigs+'_Shift_'+band2Name + imgtype)
        
        if rbm==1:
            fig_RBM.savefig('RBM.svg')
            if peaks == 1:
               fig_peaks.savefig('Peak Data.svg')
               
        if raw==1:
            fig_raw.savefig(nameFigs+'_RawSpectra'+imgtype)
        if norm==1:
            fig_norm.savefig(nameFigs+'_NormalizedSpectra'+imgtype)
        if rng==1:
            fig_rng.savefig(nameFigs+'_SpectralFeatureRanges'+imgtype)
    
       
        if correlation1==1:
            fig_Corr1.savefig(nameFigs+'_'+band1Name+'_vs_I'+ band2Name + band1Name + imgtype)
        if correlation2==1:            
            fig_Corr2.savefig(nameFigs+'_'+band2Name+'_vs_I'+ band2Name + band1Name + imgtype)
        if correlation3==1:            
            fig_Corr3.savefig(nameFigs+'_'+band3Name+'_vs_I'+ band2Name + band1Name + imgtype)
        if correlation4==1:            
            fig_Corr4.savefig(nameFigs+'_shift_'+ band3Name + 'vs_shift_'+ band1Name+imgtype)   
        if correlation5==1 and lorentz==1:            
            fig_Corr5.savefig(nameFigs+'_FWHM_vs_shift'+imgtype)    
   
        if lorentz == 1:
            fig_peak2_FWHM.savefig(nameFigs+'_FWHM_'+band2Name + imgtype)
            fig_peak3_FWHM.savefig(nameFigs+'_FWHM_'+band3Name + imgtype)

            if nt==1:
                fig_peak1_I.savefig(nameFigs+'_IntensityRatio_'+band1Name+'+'+band1Name+'-'+imgtype)
                fig_peak1_FWHMplus.savefig(nameFigs+'_FWHM_'+band1Name +'+'+ imgtype)
                fig_peak1_FWHMplus.savefig(nameFigs+'_FWHM_'+band1Name +'-'+ imgtype)
            else:
                fig_peak1_FWHMplus.savefig(nameFigs+'_FWHM_'+band1Name + imgtype)
                
                
    #%% Pick Event for 2D Map
    

    def onpick(event, fig, ax, var2D, line, spec_label, unit,red_dots_dict):
        global data_spec
        global figFC
        global axFC
        global prev_Fig
        
        fig = event.canvas.figure
        ax = event.artist.axes
    
        if 'ind' in globals():
            del ind
        try: 
            prev_Fig
        except NameError:
            prev_Fig=None
            
        if plt.fignum_exists('Spectra') and prev_Fig != fig:
            # Close Selected Spectra plot if clicked on a different 2D Map
            plt.close(figFC)
            figFC = None
            axFC = None
    
        if not plt.fignum_exists('Spectra'):
            figFC = plt.figure(num='Spectra', figsize=(10 * 1.5, 5 * 1.5))
            axFC = figFC.add_subplot(111)
            data_spec = pd.DataFrame()
            data_spec['Raman Shift(cm-1)'] = Shift
            if  prev_Fig == fig:
                prev_red_dots = red_dots_dict.get(fig, [])
                for dot in prev_red_dots:
                    dot.remove()
                red_dots_dict[fig] = []
                
        
        if prev_Fig is not None and prev_Fig != fig:
            prev_red_dots_dict = getattr(prev_Fig.canvas, 'red_dots_dict', {})
            if prev_Fig in prev_red_dots_dict:
                for p in prev_red_dots_dict[prev_Fig]:
                    p.remove()
                del prev_red_dots_dict[prev_Fig]
                prev_Fig.canvas.draw()

    
        if event.artist != line or not len(event.ind):
            return True
    
        ind = event.ind[0]

        if fig in red_dots_dict:
            red_dots = red_dots_dict[fig]
        else:
            red_dots = []
        red_dot = ax.plot(xx[ind], yy[ind], 'ro')
        red_dots.append(red_dot[0])
        red_dots_dict[fig] = red_dots
    
        if data_spec is not None:
            labp = 'x:' + str(int(xx[ind])) + ', y:' + str(int(yy[ind]))
            data_spec[labp] = Intensity_norm[:, ind]
                
       # pts.append(ax.plot(xx[ind],yy[ind],'ro'))

        if use_leng!=1:
            labp='x:'+str(int(xx[ind]))+', y:'+str(int(yy[ind]))
        else:
            labp='x:'+str(round_sig(xx[ind],3))+'$\mu m$, y:'+str(round_sig(yy[ind],3))+'$\mu m$'
            

        data_spec[labp]=Intensity_norm[:,ind]
         
        
        print('\nSpectrum Index: ',ind,'\n')

        #Plot spectra for selected pixel
        if use_leng!=1:
            axFC.plot(Shift,Intensity_norm[:,ind],
                      label=labp+', '+spec_label +': '+str(round_sig(var2D[ind],5))+unit,
                      markersize=1)
            print('Pixel Location x:',int(xx[ind]),' y:',int(yy[ind]),'\n')

        else:
            axFC.plot(Shift,Intensity_norm[:,ind], label=labp+', '+spec_label+': '+str(round_sig(var2D[ind],5))+unit, markersize=1)
            print('Pixel Location :'+str(round_sig(xx[ind],3))+' $\mu m$, y:'+str(round_sig(yy[ind],3))+' $\mu m$ \n')

        axFC.set_xlabel('Raman shift / $cm^{-1}$',fontsize=fs)
        axFC.set_ylabel('Intensity',fontsize=fs)
        axFC.legend(fontsize=fs-1)
        axFC.set_title('Selected Normalized Spectra on '+labels[z],fontsize=fs+1)

        axFC.tick_params(axis="both", labelsize=fs)   
        print('Intensity Ratio is '+str(round_sig(I21[ind],5))+'\n')
        if nt==1:
            print('G+ band Shift is '+str(round_sig(center_1plus[ind],5))+' cm-1 \n')
        else:
            print('G band Shift is '+str(center_1plus[ind])+' cm-1 \n')
        print('D band Shift is '+str(center_2[ind])+' cm-1 \n')
        print('2D band Shift is '+str(center_3[ind])+' cm-1  \n')
        
        data_spec.to_csv('Selected_Spectra.csv',index=False)
        fig.canvas.draw()

        if figFC is not None:
            figFC.canvas.draw()
            figFC.savefig('SelectedSpectra' + imgtype)
        prev_Fig = fig
        fig.canvas.red_dots_dict = red_dots_dict

        
        return True
    
    


        
    if map_==1: 
        red_dots_dict={}
        fig_MapI21.canvas.red_dots_dict = {}
        fig_Map_peak1.canvas.red_dots_dict = {}
        fig_Map_peak2.canvas.red_dots_dict = {}
        fig_Map_peak3.canvas.red_dots_dict = {}

        fig_MapI21.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_MapI21, ax_MapI21, var2D_I21, lineI21, spec_labelI21, unitI21,red_dots_dict))
        fig_Map_peak1.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map_peak1, ax_Map_peak1, var2D_peak1, line_peak1, spec_label_peak1, unit_peak1,red_dots_dict))
        fig_Map_peak2.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map_peak2, ax_Map_peak2, var2D_peak2, line_peak2, spec_label_peak2, unit_peak2,red_dots_dict))
        fig_Map_peak3.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map_peak3, ax_Map_peak3, var2D_peak3, line_peak3, spec_label_peak3, unit_peak3,red_dots_dict))

        if map1==1:
            fig_Map1.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map1, ax_Map1, var2D_1, line1, spec_label1, unit1, red_dots_dict))
        if map2==1:
            fig_Map2.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map2, ax_Map2, var2D_2, line2, spec_label2, unit2, red_dots_dict))
        if map3==1:
            fig_Map3.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map3, ax_Map3, var2D_3, line3, spec_label3, unit3, red_dots_dict))
        if map4==1:
            fig_Map4.canvas.mpl_connect('pick_event', lambda event: onpick(event, fig_Map4, ax_Map4, var2D_4, line4, spec_label4, unit4, red_dots_dict)) 

