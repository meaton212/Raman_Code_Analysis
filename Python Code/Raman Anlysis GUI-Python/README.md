# Raman Analysis Python Code w/ GUI
Date: July 19, 2023
Python Code: RamanAnalysis_nanomaterials_GUI.py, analysisInputs.py, LorentzInputs.py, LorentzInputsDouble1.py, PlotsInputs.py, raman_data_processing.py, RBMinputs.py

Cite as: [DOI: 10.26434/chemrxiv-2023-5nhdj](https://chemrxiv.org/engage/chemrxiv/article-details/645cc378fb40f6b3ee660d89)

All codes and sample data can be found at IMDEA nanoscience institutional repository: https://repositorio.imdeananociencia.org/handle/20.500.12614/3316
Contact: natalia.martin@imdea.org

## Summary

This code is developed to help material scientist in the analysis of Raman spectra of nanomaterials (graphene, carbon nanotubes, TMDCs etc). In the paper Automated Statistical analysis of Raman spectra of nanomaterials, the user can find a detailed summary of the code structure and its capabilities. This document provides instructions to use the codes and the graphical user interface. 

This code processes Raman spectral data files consisting of multiple spectra and analyzes each curve individually for their relevant spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample. The total number of bands analyzed are 3 (referred as peak 1, 2 and 3 in the code)

It has the option of finding these bands through the maximum value in a defined range (method 1), or through the fitting of a lorentizan distribution to the spectral data (method 2). Additionally, the peak 1 can be split into two bands if the lorentzian fitting method is chosen (for example to capture G<sup>+</sup> and G<sup>-</sup> splitting in carbon-based nanomaterials). The user also has the option to locate and plot the statistical distribution of the RBM modes (radial breathing modes).

The code will also compare several spectral features in the following correlation plots (selected by the user):

*	Position of the 3 bands vs the intensity ratio I<sub>2</sub>/I<sub>1</sub>
*	Peak 3 position vs. peak 1 position

If lorentzian peak fitting is performed (method 2), additional correlations are plotted:
*	Position vs FWHM of peaks 1, 2 and 3

Additionally, the Raman code will plot the spectra for each file in several plots, if the user wishes
*	Raw Spectra
*	Normalized Spectra
*	Average Spectra for each file
*	Spectral Range for Intensity Ratio calculation
*	RBM peak identification results

Finally, the code can also plot a contour map of the spectral features in 2D for any spatially distributed Raman spectra. Specifically, color maps of the position of the 3 peaks, widths and the intensity ratios I<sub>2</sub>/I<sub>1</sub>, I<sub>3</sub>/I<sub>1</sub> will be plotted. The user can plot Raman maps of specific spectral features for samples where the spatial arrangement is relevant, such as in the case of controlled patterning of nanomaterials.

## Data Format

* Data can be saved in any text based file format, including, but not limited to .txt, .dpt, .csv, and .dat. The file type is defined by the user in the input section.
* Currently, the data files must have no header, and only contain the raw data.
* The first column in the data is the raman shift (in cm<sup>-1</sup>) and must be the same for all the spectra.
* Each additional column is the intensity data corresponding the raman shift in the first column. Each column is a unique spectra.
* The data must be separated by tab, space, comma (,), semicolon (;), or forward-slash (/) delimiter. The corresponding delimiter is selected by the user in the input section


## Inputs
A graphical interface has been created where the user should input their preferred analysis parameters. Be aware that all <b>the following files have to be saved together in the same folder</b> in order the program to work:

* RamanAnalysis_nanomaterials_GUI.py —> Main code where all the analysis is performed. This algorithm calls in the following functions that control the graphical interface:
* AnalysisInputs.py
* LorentzInputs.py
* LorentzInputsDouble1.py
* PlotInputs.py
* RBMInputs.py

### 1. Initializing:
Starting from the main code (RamanAnalysis_nanomaterials_GUI.py), upon pressing the RUN command the program will be initiated.
* A pop-up window will appear telling the user to select the file containing directory. Press Ok to continue
* Then, select the folder where the data files to be analyzed are stored.

### 2. Data Selection
* Choose the file type extension and the delimiter type from the dropdown menus, and press "Submit"
* In the next window, press "open files" and select the desired files that you want to analyze. No more than 10 files may be selected
* You can remove selected files by clicking the chosen file and clicking "Remove Selected"
* Once all files have been selected, press "Continue"
* In the next window, input the desired legend label for each of your selected files, and press "Save Labels"

### 3. Analysis Method

* 2 windows will pop up. First, an input panel where the user must add their preferred analysis parameters, and a second window containing a plot  of the average spectra of each of the samples. The plots are aimed to help the user to select the appropriate spectral range for their data

#### Panel 1: Analysis method
Select the analysis method. Method 1 is a simple approximation where peak location is defined as the point of maximum intensity within a range, while method 2 performs lorentzian peak fitting (width of the peaks are also retrieved with this method). Peak 1 can be splitted into 2 bands (for example G<sup>+</sup> and G<sup>-</sup> in carbon-based materials) by clicking the corresponding box. Radial breathing modes analysis can be selected here by clicking the corresponding option.

#### Panel 2: Mapping options
Select if Raman maps of the different spectral features are desired by clicking the "Mapping" box. In the left column, introduce the number of pixels in X and Y. If the dimensions of each dimension are known, click the "Use Dimension" box. In the right column, indicate the corresponding distances on X and Y for the map (in micrometers).

#### Panel 3: Peak information. 
A total of 3 peaks can be analyzed, and peak 1 can be split into 2 peaks (if the corresponding box is selected in panel 1, and lorentzian peak fitting is chosen). Input the lower and upper bounds of the spectral range (in cm<sup>-1</sup>) where the maximum intensity value will be used to find the peaks location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function. Moreover, the name of the peaks should be introduced here (ie G, D, 2D etc)
After pressing ‘OK’, a figure will be generated plotting the averaged and normalized spectra for each of the samples.

### 4. RBM inputs:
If the RBM analysis is chosen in previous step, a window will pop up asking for the spectral range where the peaks are located and the prominence value. This value sets the max limit at which peaks will be considered for the RBM analysis. Local maxima with a maximum height lower than Prominence will be discarded. Input a prominence value for each of the files analyzed separated by commas: Prom_sample1, Prom_sample2, Prom_sample3, … *Only relevant if RBM=1

### 5. Initial estimates for lorentzian fitting: 
If method 2 (lorentzian fit) is selected in the analysis method step, a window will open requesting initial estimates for the position, FWHM and maximum intensity of each of the peaks. Input a list of 3 parameters contained within brackets: [center, FWHM, and Max_intensity], in this order and separated by commas for each analyzed file. Initial guesses for different files should be separated by commas (,).

### 6. Plotting options: 
Last, the user can decide whether or not to produce and save specific plots in addition to the default ones. Click the desired boxes. 
Moreover, select the width of the bins for the histograms showing the calculated statistics of the samples. Different widths can be chosen for the histograms of Raman shift (in cm<sup>-1</sup>), FWHM (in cm<sup>-1</sup>) and intensity ratio.
Finally, click the corresponding box for saving the generated figures (in selected image format) in the current datafolder. Type the basename for the saved figures.










