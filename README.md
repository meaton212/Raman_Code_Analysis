# Raman Analysis Code

The code is available in Both Python and MATLAB. Please select the relevent folder for the specific code. Each code has its own README file outlining the important information specific to the code.

## General Summary

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relavent spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample, namely:

* G band
* D band
* 2D band
* Intensity ratio between the D and G bands

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

