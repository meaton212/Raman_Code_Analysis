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



