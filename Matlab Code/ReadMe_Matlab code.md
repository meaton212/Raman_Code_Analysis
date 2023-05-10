**Raman Analysis Matlab code**

Date: 09-05-2023

Matlab code: RamanAnalysis\_nanomaterials.m

Cite as: XXX

**SUMMARY**

This code is developped to help material scientist in the analysis of Raman spectra of nanomaterials (graphene, carbon nanotubes, TMDCs etc). Within the code, anotations have been made to guide the user to its usage. Moreover, in the paper XXX, the user can find a detailed summary of the code structure and its capabilities. The code is included in a .m Matlab file that is ready to use.** Note that the axis labels and titles of the graphs are programmed for analysis of carbon-based materials. For analysis of other types of nanomaterial read last section in this file.

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relevant spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample. The total number of bands analysed are 3, which in the case of carbon based nanomaterials would be:

- G band (peak 1)
- D band (peak 2)
- 2D band (peak 3)

It has the option of finding these bands through the maximum value in a defined range (method 1), or through the fitting a lorentizan distribution to the spectral data (method 2). Additionally, the G band can be split into G+ and G- bands if the lorentzian fitting method is chosen. The user also has the option to locate and plot the statistical distribution of the RBM (radial breathing modes) modes (optional).

The code will also compare several spectral features in the following correlation plots:

- G band vs. Intensity ratio between the D and G bands (Id/Ig)
- D band vs. Id/Ig
- 2D band vs. Id/Ig
- 2D band position vs. G band position

`	`If lorentzian peak fitting is performed (method 2), additional correlations are plotted:

- Shift vs FWHM of G band (G+ in the case of split analysis)
- Shift vs FWHM of D band
- Shift vs FWHM of 2D band

Additionally, the Raman code will plot the spectra for each file in several plots, if the user wishes

- Raw Spectra
- Normalized Spectra
- Average Spectra for each file
- Spectral Range for Intensity Ratio

Finally, the code can also plot a heat map of the spectral features in 2D for any spatially distributed Raman spectra. Specifically, color maps of the position of the 3 peaks and the intenty ratio Id/Ig will be plotted. The user can plot Raman maps of specific spectral features for samples where the spatial arrangement is relevant, such as in the case of controlled patterning of nanomaterials.

**DATA FORMAT**

- Data can be saved in .txt or .dat format. The file type is defined by the user in the type variable in the input section
- Currently, the data files must have no header, and only contain the raw data. 
- The first column in the data is the raman shift (in cm-1) and must be the same for all the spectra.
- Each additional column is the intensity data corresponding the the raman shift in the first column. Each column is a unique spectra.
- The data must be separated by tab or space delimiter

**INPUTS**

User Input section goes from line 4 to line 86 of the code. Annotations have been made in green to guide the user in the completion of the inputs. The following variables must be changed by the user before running the code:

**Data Information (lines 4-20)**

- path: directory where the files to be analyzed are stored. Line 13
- file\_name: list of filenames, without file extension, that user desires to be analysed, in the order that the user wants them to be processed. Names should be separated by ; and written between ‘’. Line 14
- total: Number of files in filenames that you want analyzed. Number must be equal to or less than the length of file\_name. Line 16
- type: file type/extension for the data file (i.e. '.txt' or ‘.dat’). Line 17
- name: list of names used in the legend for each of the files in file\_names. Line 18

`	   `**Model used for analysis (lines 22-28)**

- rbm: set to 1 if Radial breathing mode (RBM) analysis is desired. Set to 0 to ignore. Note, that if the number of spectra is very high, the computing is going to slow down significantly
- lorentz: Set to 1 to fit a lorentzian model each of your spectral bands (G band, D band, 2D band). Set to 0, to find the band based on solely the maximum intensity value in defined bounds.
- nt: Set to 1 to split G band into G+ and G- bands. Set to 0 to only find single G band.
- map: Set to 1 to plot 2D maps of spectral features (positions and intesity ratio). Set to 0 to ignore.



`	  `**Spectral bounds: Normalization and peak identification (lines 29-53)**

- normLow & normHigh: Input the lower and upper bounds of the spectral range (in cm-1) where the maximum intensity value will be used to normalize the spectra
- band1Low & band1High: Input the lower and upper bounds of the spectral range (in cm-1) where the maximum intensity value will be used to find the G band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.
- band2Low & band2High: Input the lower and upper bounds of the spectral range (in cm-1) where the maximum intensity value will be used to find the D band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.
- band3Low & band3High: Input the lower and upper bounds of the spectral range (in cm-1) where the maximum intensity value will be used to find the 2D band location. If lorentzian fits are chosen, then the code will fit this shift range to a lorentzian distrubution function.
- RBMregion\_Low & RBMregion\_Low: Input the lower and upper bounds of the spectral range (in cm-1) where the RBM analysis will occur \*Only relevant if RBM=1
- Prom: This value sets the max limit at which peaks will be considered for the RBM analysis. Local maxima with a maximum height lower than Prom will be discarded. \*Only relevant if RBM=1

`	`**Fitting parameters (lines 54-66)** \*Only relevant if lorentz=1

- Gmin\_init: Initial Guess for the lorentzian fit parameters (3) of the G- band, only if nt=1. There sould be a list of 3 parameters (separated by a space) for each analyzed file. Initial guesses for different files should be separated by ;
- Gplus\_init: Initial Guess for the lorentzian fit parameters (3) of the G+ band, if nt=1, or G band, if nt=0. There sould be a list of 3 parameters (separated by a space) for each analyzed file. Initial guesses for different files should be separated by ;
- D\_init: Initial Guess for the lorentzian fit parameters (3) of the D band. There sould be a list of 3 parameters (separated by a space) for each analyzed file. Initial guesses for different files should be separated by ;
- init\_2D: Initial Guess for the lorentzian fit parameters (3) of the 2D band. There sould be a list of 3 parameters (separated by a space) for each analyzed file. Initial guesses for different files should be separated by ;

`	`**Mapping options (lines 67-75)** \*Only relevant if map=1

- col: number of columns in 2D map. Note rows x col must equal total number of spectra in file.
- raws: Number of rows in 2D map. Note rows x col must equal total number of spectra in file.
- colmum: length of x axis in μm
- rawsmum: length of y axis in μm

`	`**Plotting Options (lines 76-86)**

**	User can decided whether or not to produce and save specific plots by assigning 1 (yes) or 0 (no) to these variables

- raw: Set to 1 to plot all raw spectra for each of the files. Set to 0 to ignore.
- norm: Set to 1 to plot all normalized spectra for each of the files. Set to 0 to ignore.
- range: Set to 1 to plot the spectral regions chosen for G, D and 2D bands in intensity calculation for each of the files. Set to 0 to ignore.
- peaks: Set to 1 to plot found peaks in the RBM region and save the plot. Set to 0 to ignore.
- correlations: Set to 1 to Plot correlations between peaks shifts, FWHM and Id/Ig.
- width: width of Raman shift histograms (in cm-1)
- width\_fw: width of FWHM histograms (in cm-1). \*Only relevant if lorentz=1


**OUTPUTS**

**Output Text Files: spectral features results**

After analyzing all the data, the code will output a .csv file for each file analyzed . The resulting file will contain the spectral features results for each spectrum analyzed in the file. The file is saved as "name\_results.txt”

The columns are labelled and contain, depending on the analysis performed, intensity, shift (in cm-1), FWHM (in cm-1, only if lorentz=1) of G (G- and G+ if nt=1), D and 2D modes. The RBM modes result are included in the last columns, with intensity and shift for each RBM mode found. Note that different spectra within a file might have different number of RBM (depending on the value of prom defined). The code saves as many columns as the maximum number of RBMs found in the spèctra. In the case that one specfic spectrum has less RBMs, then the value will be blanck in the corresponding column.

**Figures**

Several figures will be produced, with appropriate titles, labels and legends. Succesive files will be plotted in different colors within the same graph. Figures will not be automatically saved, so the user has to select and saved the desired figures.

`	`**Spectra:** The plots of spectra (intensity vs. raman shift (cm-1)) will be saved as the following:

- Figure 1: Raw data (if raw = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. 
- Figure 2: Normalized spectra (if norm = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. 
- Figure 3: Average spectra of each file. A single curve for each file on a single plot.
- Figure 10: The spectral regions chosen for G, D and 2D bands in intensity calculation (if range=1). Each spectra in a file is plot on a single subplot, with a subplot for each file.
- Figures 101, 201 & 301: The Lorentian fit curves overlayed over normalized spectra for the G band, D band, and 2D band (if lorentz = 1). Each band is plot in a subplot, with each spectra for each region in each subplot. There is a separate graphic for each file.
- Figure 401: The RBM regions with the selected peaks chosen for each spectra marked (if rbm = 1 and peaks = 1). Each spectra in a file is plot on a single subplot, with a subplot for each file. 

`	`**Histograms:** Histograms of the distributions of shifts, intensity ratio and FWHM (if lorentz=1) for all spectra within every file. Different files will be plotted in the same graph to facilitate comparison and in different colors. Mean values and standard deviations are included in the legends of each plot. Each histogram will contain the data for all files analysized in a single histogram.

- Figure 11: Intensity ratio, Id/Ig 
- Figures 100, 200 & 300: Raman shift of  G (or G+ and G-, if nt = 1), D and 2D modes. 
- Figures 98, 199 and 299: FWHM of G (or G+ and G-, if nt = 1), D and 2D modes. 
- Figure 400: RBM peaks position (if rbm = 1) 
- Figure 99: G band Intensity ratio, IG+/IG- , if nt = 1.

`	`**Correlations:** If correlations = 1, the following figures will be ploted:

- Figures 1000-1002: Shift of G  (G+ if nt=1), D and 2D vs Intensity ratio Id/Ig
- Figure 1003: Shift of 2D band vs. G band shifts (G+ if nt=1)
- Figures 1004-1006: FWHM vs shift for G (G+ if nt=1), D and 2D bands (Only if lorentz=1)

The data for each spectra is plotted as a point. The data for each file is fit to a linear regression, with the linear regression line overlayed with the data and the fitting parameters included in the legend. The data for each file is plot on the same graph, with a different color for each file.

`	`**2D Map:** If map = 1, then 2D heat maps of the spectral features will be produced. 

- Figure 40: Id/Ig 
- Figure 41: Shift G band  (G+ if nt=1)
- Figure 42: Shift D band
- Figure 43: Shift 2D band
- Figure 44: Shift G- band  (G+ if nt=1)

**GENERALIZATION TO OTHER NANOMATERIALS:**

The code can be used to analyze any type of material, and it is not limited to carbon-based nanomaterials (however the axis labels and titles are currently set for the G, D, 2D and RBM modes of carbon-based materials). In order to generalize the code, the user can simply change the spectral range inputs (see INPUTS-Spectral bounds/Peak identification above) to the desired regions. Note that 3 peaks can be analysed, labelled peak 1, 2 and 3 in the code. The 3 of them can be analysed by the maximum value in a defined range (method 1) of by lorentzian peak fitting (method 2, lorentz=1). Peak 1 can be further fit into 2 lorentzian peaks (nt=1).

The user will need to manually change figure labels and titles.
