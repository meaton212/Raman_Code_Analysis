# Raman Analysis Python Code w/ GUI
Date: July 19, 2023
Python Code: RamanAnalysis_GUI.py, analysisInputs.py, LorentzInputs.py, LorentzInputsDouble1.py, PlotsInputs.py, raman_data_processing.py, RBMinputs.py

Cite as: [DOI: 10.26434/chemrxiv-2023-5nhdj](https://chemrxiv.org/engage/chemrxiv/article-details/645cc378fb40f6b3ee660d89)

All codes and sample data can be found at IMDEA nanoscience institutional repository: https://repositorio.imdeananociencia.org/handle/20.500.12614/3316
Contact: natalia.martin@imdea.org

## Summary

This code is developped to help material scientist in the analysis of Raman spectra of nanomaterials (graphene, carbon nanotubes, TMDCs etc). In the paper Automated Statistical analysis of Raman spectra of nanomaterials, the user can find a detailed summary of the code structure and its capabilities. This document provides instructions to use the codes and the graphical user interface. 

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relevant spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample. The total number of bands analysed are 3 (referred as peak 1, 2 and 3 in the code)

It has the option of finding these bands through the maximum value in a defined range (method 1), or through the fitting of a lorentizan distribution to the spectral data (method 2). Additionally, the peak 1 can be split into two bands if the lorentzian fitting method is chosen (for example to capture G+ and G- splitting in carbon-based nanomaterials). The user also has the option to locate and plot the statistical distribution of the RBM modes (radial breathing modes).

The code will also compare several spectral features in the following correlation plots (selected by the user):

*	Position of the 3 bands vs the intensity ratio I2/I1
*	Peak 3 position vs. peak 1 position

If lorentzian peak fitting is performed (method 2), additional correlations are plotted:
*	Position vs FWHM of peaks 1, 2 and 3

Additionally, the Raman code will plot the spectra for each file in several plots, if the user wishes
*	Raw Spectra
*	Normalized Spectra
*	Average Spectra for each file
*	Spectral Range for Intensity Ratio calculation
*	RBM peak identification results

