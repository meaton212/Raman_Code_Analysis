# Raman Analysis Python Code


## Summary

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relavent spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample, namely:

* G band
* D band
* 2D band
* Intensity ratio between the D and G bands

It has the option of finding these bands through the maximum value in a defined range, or through the fitting a lorentizan distribution to the spectral data

Additionally, the G band can be split into G+ and G- bands

The user also has the option to locate and plot the statistical distribution of the RBM modes (optional)

The code will also compare several spectral features in the following correlation plots:
* G band vs. Intensity ratio between the D and G bands
* D band vs. Intensity ratio between the D and G bands
* 2D band vs. G band

Additionally, the Raman code will plot the spectra for each file in several plots, if the user wishes
* Raw Spectra
* Normalized Spectra
* Average Spectra for each file
* Spectral Range for Intensity Ratio

Finally, the code can also plot a heat map of any of the spectral features in 2D for any spatially distributed Raman spectra. Only the final file will be interactive for the user so they can select specifc pixels in the map and overlay different pixels' spectra in a single plot.

The user can plot Raman maps of specific spectral features for samples where the spatial arrangement is relevant, such as in the case of controlled patterning of nanomaterials.  


## Data format

Data can be saved in any text based file, such as .txt, .csv, .dpt, etc.... The file type is defined by the user in the type_ variable in the input section

Currently, the data files must have no header, and only contain the raw data

The first column in the data is the raman shift (in cm<sup>-1</sup>) and must be the same for all the spectra.

Each additional column is the intensity data corresponding the the raman shift in the first column. Each column is a unique spectra.

The data can be separated by any delimeter. The delimeter is defined in the delim variable by the user in the input section


## Inputs
User Input section starts on Line 47
The following variables must be changed by the user before running the code

* path: directory where the files to be analyzed are stored.
* file_name: list of filenames, without file extension, that use desires to be processed, in the order that the user wants them to be processed
* total: Number of files in filenames that you want analyzed. Number must be equal to or less than the length of file_name
* delim: delimiter of the data in the data file
* type_: file type/extension for the data file (i.e. '.txt', '.csv', etc...)
* name: list of names used in the legend for each of the files in file_names
* imgtype: the type/extension for the output images (i.e. '.svg', '.jpg', '.png')
* dens: (True or False) Select True for histograms to be normalized, select False for non-normalized histograms

### Spectral Bounds
* normLow & normHigh: Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the maximum intensity value will be used to normalize the spectra
* band1Low & band1High: Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the maximum intensity value will be used to find the G band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.
* band2Low & band2High: Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the maximum intensity value will be used to find the D band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.
* band3Low & band3High: Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the maximum intensity value will be used to find the 2D band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.

### RBM
* rbm: set to 1 if Radial breathing mode (RBM) analysis is desired. Set to 0 to ignore. Note, that if the number of spectra is very high, the computing is going to slow down significantly
* RBMregion_Low & RBMregion_Low: Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the RBM analysis will occur
* Prom: This value sets the max limit at which peaks will be considered. This can be a single number, or a list of values that different for each file. Local maxima with a maximum height lower than Prom will be discarded.

### Lorentzian Fitting and G band Splitting Options
* lorentz: Set to 1 to fit a lorentzian model each of your spectral bands (G band, D band, 2D band). Set to 0, to find the band based on solely the maximum intensity value in defined bounds.
* nt: Set to 1 to split G band into G<sup>+</sup> and G<sup>-</sup> bands. Set to 0 to only find single G band.
* Gmin_init: Initial Guess for the lorentzian fit parameters (3) of the G<sup>-</sup> band, only if nt=1. There sould be a list of 3 parameters for each analyzed file.
* Gplus_init: Initial Guess for the lorentzian fit parameters (3) of the G<sup>+</sup> band, if nt=1, or G band, if nt=0.There sould be a list of 3 parameters for each analyzed file.
* D_init: Initial Guess for the lorentzian fit parameters (3) of the D band. There sould be a list of 3 parameters for each analyzed file.
* init_2D: Initial Guess for the lorentzian fit parameters (3) of the 2D band. There sould be a list of 3 parameters for each analyzed file.

### Plotting Options

User can decided whether or not to produce and save specific plots by assigning 1 or 0 to these variables
* raw: Set to 1 to plot all raw spectra for each of the files and save the plot. Set to 0 to ignore.
* norm: Set to 1 to plot all normalized spectra for each of the files and save the plot. Set to 0 to ignore.
* rng: Set to 1 to plot the spectral regions chosen for  G, D and 2D bands in intensity calculation for each of the files and save the plot. Set to 0 to ignore.
* peaks: Set to 1 to plot found peaks in the RBM region and save the plot. Set to 0 to ignore.
* correlations: Set to 1 to Plot correlations between peaks shifts and I<sub>d</sub>/I<sub>g</sub> and save the plot.

### Plot Format
* fs: Fontsize for the texts used in all the plots
* width: Width of histograms

### 2D Maps
* maps: Set to 1 to plot a 2D heatmap of the selected variable defined by map_var. Set to 0 to ignore.
* map_var: Spectral variable chosen for 2D maps. Must be 'I' for intensity ratio, 'G' for G band (or G<sup>+</sup> if nt=1), 'D' for D band, or '2D' for 2D band. 
* rows: Number of rows in 2D map. Note rows x col must equal total number of spectra in file.
* col: number of columns in 2D map. Note rows x col must equal total number of spectra in file.
* use_length: Set to 1 if you want to use the user defined dimensions for the x and y axes (in μm). If set to 0, then the x and y axes will count the pixel numbers.
* xlen: length of x axis in μm
* ylen: length of y axis in μm



## Outputs

### Output Text Files


### Images
Image type is determined by the extension defined in the imgtype variable

##### Histograms
Histograms of the followin







