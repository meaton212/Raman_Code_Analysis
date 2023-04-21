# Raman Analysis Python Code


## Summary

This code processes Raman spectral data files consisting of mulitiple spectra and analyzes each curve indivually for their relavent spectral features and provides the correspond statistical distributions of the position, widths and intensity of the main Raman bands within a sample, namely:
-G band
-D band
-2D band
-Intensity ratio between the D and G bands

It has the option of finding these bands through the maximum value in a defined range, or through the fitting a lorentizan distribution to the spectral data

Additionally, the G band can be split into G+ and G- bands

The user also has the option to locate and plot the statistical distribution of the RBM modes (optional)



The user can plot Raman maps of specific spectral features for samples where the spatial arrangement is relevant, such as in the case of controlled patterning of nanomaterials.  


## Data format

Currently, the data files must have no header, and only list the data in 
Delimeter can be defined in the 


## Inputs
User Input section starts on Line 47
The following variables must be changed by the user before running the code

path: directory where the files to be analyzed are stored.
file_name
