# Raman Analysis Python Code


## Summary

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relavent spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample, namely:

* G band
* D band
* 2D band
* Intensity ratio between the D and G bands (I<sub>d</sub>/I<sub>g</sub>)

It has the option of finding these bands through the maximum value in a defined range, or through the fitting a lorentizan distribution to the spectral data

Additionally, the G band can be split into G+ and G- bands

The user also has the option to locate and plot the statistical distribution of the RBM modes (optional)

The code will also compare several spectral features in the following correlation plots:
* G band vs. Intensity ratio between the D and G bands (I<sub>d</sub>/I<sub>g</sub>)
* D band vs. I<sub>d</sub>/I<sub>g</sub>
* 2D band vs. G band

Additionally, the Raman code will plot the spectra for each file in several plots, if the user wishes
* Raw Spectra
* Normalized Spectra
* Average Spectra for each file
* Spectral Range for Intensity Ratio

Finally, the code can also plot a heat map of any of the spectral features in 2D for any spatially distributed Raman spectra. Only the final file will be interactive for the user so they can select specifc pixels in the map and overlay different pixels' spectra in a single plot.

The user can plot Raman maps of specific spectral features for samples where the spatial arrangement is relevant, such as in the case of controlled patterning of nanomaterials.  

## Python Packages

To run this code, make sure the following Python packages are installed:
* pandas
* numpy
* os
* matplotlib
* scipy
* tqdm

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

##### Spectral Feature Data
After analyzing all the data, the code will output a .csv file for each file analyzed. The resulting file will contain the spectral feature data for each spectrum analyzed in the file. The columns are: 
* Spectra #: The index for each spectrum
* Intensity_G: The intensity at the G band shift (or G<sup>+</sup> if nt = 1)
* Shift_G: The Raman shift of the G band in cm<sup>-1</sup> (or G<sup>+</sup> if nt = 1)
* Intensity_D: The intensity at the D band shift
* Shift_D: The Raman shift of the D band in cm<sup>-1</sup>
* Intensity_2D: The intensity at the 2D band shift
* Shift_2D: The Raman shift of the 2D band in cm<sup>-1</sup>
* Intensity Ratio: The intensity ratio between the D and G bands I<sub>d</sub>/I<sub>g</sub>

The file is saved as "name_results.csv"

##### Selected Spectra Data from 2D Map
If maps = 1, then the 2D heat map of the spectral feature defined in the map_var variable will be mapped in a 2D heat map.

For the last file analyzed, the user can select specfic pixels from the map, and plot the spectra of the selected pixels in a new figure. The user can select multiple pixels and overlay the spectra. 

The selected spectra of the most recently overlayed selected pixel spectra will be output in a .csv file containing the full spectra selected. The first column will be the raman shift in cm<sup>-1</sup> and each subsequent column are the intensity values for each selected spectrum, with the header giving the x and y positions of the pixel.

The data is saved as "Selected_Spectra.csv"


### Figure Images
Image type is determined by the extension defined in the imgtype variable

##### Histograms
Histograms of the distributions for all curves for each file for the following spectral data will be saved:
* Intensity ratio, I<sub>d</sub>/I<sub>g</sub> as "IntensityRatio.imgtype"
* G band (or G<sup>+</sup> and G<sup>-</sup>, if nt = 1) as "Gband.imgtype"
* D band as "Dband.imgtype"
* 2D band as "2DBand.imgtype"
* RBM (if rbm = 1) as "RBM.imgtype"
* G band Intensity ratio, I<sub>G<sup>+</sup></sub>/I<sub>G<sup>-</sup></sub> as "IntensityRatioG+G-.imgtype", if nt = 1.

Each histogram will contain the data for all files analysized in a single histogram

##### Correlations
If correlations = 1, a single graphic with 3 subplots of the following relationships will be saved:
* G band vs. I<sub>d</sub>/I<sub>g</sub>
* D band vs. I<sub>d</sub>/I<sub>g</sub>
* 2D band vs. G band

The data for each spectra is plotted as a point. The data for each file is fit to a linear regression, with the linear regression line overlayed with the data. The data for each file is plot on the same graph, with a different color for each file. 

##### Spectra
The plots of spectra (intensity vs. raman shift (cm<sup>-1</sup>)) will be saved as the following:
* Average spectra of each file. A single curve for each file on a single plot. Saved as "Avg Spectra.imgtype"
* raw data (if raw = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. Saved as "Raw Data.imgtype"
* normalized spectra (if norm = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. Saved as "Normalized Data.imgtype"
* The spectral regions chosen for G, D and 2D bands in intensity calculation (if rng=1). Each spectra in a file is plot on a single subplot, with a subplot for each file. Saved as "SpecRangeIR.imgtype"
* The RBM regions with the selected peaks chosen for each spectra marked (if rbm = 1 and peaks = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. Saved as "Peak Data.imgtype"
* The Lorentian fit curves overlayed over normalized spectra for the G band, D band, and 2D band (if lorentz = 1). Each band is plot in a subplot, with each spectra for each region in each subplot. There is a separate graphic for each file. Saved as "Lorentz_fit_"ame.imgtype"

##### 2D Map
If maps = 1, then the 2D heat map of the spectral feature defined in the map_var variable will be saved for each file as "name_2DMap_map_var.imgtype So, for a file with filename "p-G 2.txt" and name "patterned graphene" with imgtype ".svg" and map_var "I" will save as "patterned graphene_2DMap_I.svg"

For the last file analyzed, the user can select specfic pixels from the map, and plot the spectra of the selected pixels in a new figure. The user can select multiple pixels and overlay the spectra. The Most recent figure will be saved as "SelectedSpectra.svg"



## Example Data

Example Data is provided to test out various aspects of the code

### Example Data Set 1
Tests out multiple (2) Raman files to comare the spectral data between them, as well as lorentzian fitting,  G band splitting, and RBMs

Before running code, make sure to change the path in the Input Data Section,  (Starting on line 47) if running in a differnt folder

Example:

path= '../Example Data/Example_Data_1/' (line 54)

leave default file_names

Check to see that the following variables are correctly selected
* Set maps = 0 (Line 157)
* Set lorentz = 1 (Line 124)
* Set nt = 1 (Line 125)
* Set rbm = 1 (Line 109)

Leave the rest as the default values.

### Example Data Set 2
Tests out 2D Mapping Feature for a patterned graphene sample with 32 x 32 datapoints

Before running code, make sure to change the path in the Input Data Section (Starting on line 47)

Example:

Uncomment the variable for  path (line 55)
path='../Example Data/Example_Data_2/' (line 54)

Uncomment the variable for filename below (line 62)
* file_name=['mapExample.txt']


Uncomment the variable for name below (line 71)
* name= ['patterned graphene']

Change the following Variables
* Set maps = 1 (Line 157)
* Set lorentz = 0 (Line 124)
* Set nt = 0 (Line 125)
* Set rbm = 0 (Line 109)

Leave the rest as the default data.




