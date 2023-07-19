# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 12:33:55 2023

@author: matte
"""

import tkinter as tk
from tkinter import ttk


def analysisInputs():
    # Create the main window
    window = tk.Tk()
    window.geometry("800x600")  # Set the size of the window
    
    #Inital values for Selected ranges
    normMin=1500
    normMax=1675
    peak1Min=1450
    peak1Max=1700
    peak2Min=1250
    peak2Max=1400
    peak3Min=2400
    peak3Max=2800

    ###############################
    ### Select analysis routine ###
    ###############################

    panel1 = ttk.LabelFrame(window, text='1. Analysis method: Spectral features')
    panel1.place(x=10, y=5, width=340, height=150)

    # Method dropdown
    method_var = tk.StringVar()
    method_dropdown = ttk.Combobox(panel1, textvariable=method_var, state='readonly')
    method_dropdown['values'] = ['Method 1 (Max Peak)', 'Method 2 (Lorentz Fit)']
    method_dropdown.current(0)
    method_dropdown.place(x=10, y=10, width=250)

    # Split peak checkbox
    split_peak_var = tk.BooleanVar()
    split_peak_checkbox = ttk.Checkbutton(panel1, text='Split peak 1 (for method 2)', variable=split_peak_var)
    split_peak_checkbox.place(x=10, y=50, width=300)

    # RBMs checkbox
    rbm_var = tk.BooleanVar()
    rbm_checkbox = ttk.Checkbutton(panel1, text='Radial Breathing modes', variable=rbm_var)
    rbm_checkbox.place(x=10, y=80, width=300)

    ##########################
    ### Mapping options ######
    ##########################

    panel2 = ttk.LabelFrame(window, text='2. Mapping options')
    panel2.place(x=360, y=5, width=430, height=150)

    # Mapping checkbox
    mapping_var = tk.BooleanVar()
    mapping_checkbox = ttk.Checkbutton(panel2, text='Mapping', variable=mapping_var)
    mapping_checkbox.place(x=10, y=10, width=150)
    
    # Length checkbox
    length_var = tk.BooleanVar()
    length_checkbox = ttk.Checkbutton(panel2, text='Use Dimensions', variable=length_var)
    length_checkbox.place(x=210, y=10, width=175)

    # X in pixels
    x_pixels_label = ttk.Label(panel2, text='Pixels in X (Columns)')
    x_pixels_label.place(x=10, y=40)
    x_pixels_entry = ttk.Entry(panel2)
    x_pixels_entry.place(x=180, y=40, width=50)

    # Y in pixels
    y_pixels_label = ttk.Label(panel2, text='Pixels in Y (Rows)')
    y_pixels_label.place(x=10, y=80)
    y_pixels_entry = ttk.Entry(panel2)
    y_pixels_entry.place(x=180, y=80, width=50)

    # X in microns
    x_microns_label = ttk.Label(panel2, text='X (\u03BCm)')
    x_microns_label.place(x=250, y=40)
    x_microns_entry = ttk.Entry(panel2)
    x_microns_entry.place(x=310, y=40, width=50)

    # Y in microns
    y_microns_label = ttk.Label(panel2, text='Y (\u03BCm)')
    y_microns_label.place(x=250, y=80)
    y_microns_entry = ttk.Entry(panel2)
    y_microns_entry.place(x=310, y=80, width=50)

    #########################
    ### Peak information ###
    #########################

    panel3 = ttk.LabelFrame(window, text='3. Peak information')
    panel3.place(x=10, y=160, width=780, height=420)

    # Normalization
    normalization_label = ttk.Label(panel3, text='Normalization: Spectra will be normalized to the max height within this range',font="bold")
    normalization_label.place(x=10, y=10)

    norm_min_label = ttk.Label(panel3,    text='min (cm\u207B\u00B9)', justify=tk.RIGHT)
    norm_min_label.place(x=250, y=40)
    norm_min_entry = ttk.Entry(panel3)
    norm_min_entry.place(x=350, y=40, width=100)
    norm_min_entry.insert(tk.END, str(normMin))

    norm_max_label = ttk.Label(panel3, text='max (cm\u207B\u00B9)', justify=tk.RIGHT)
    norm_max_label.place(x=510, y=40)
    norm_max_entry = ttk.Entry(panel3)
    norm_max_entry.place(x=610, y=40, width=100)
    norm_max_entry.insert(tk.END, str(normMax))

    # Peak 1
    peak1_label = ttk.Label(panel3, text='Peak 1', justify=tk.LEFT,font="bold")
    peak1_label.place(x=10, y=70)

    peak1_range_label = ttk.Label(panel3, text='Indicate spectral range:')
    peak1_range_label.place(x=30, y=100)

    peak1_min_label = ttk.Label(panel3, text='min (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak1_min_label.place(x=250, y=100)
    peak1_min_entry = ttk.Entry(panel3)
    peak1_min_entry.place(x=350, y=100, width=100)
    peak1_min_entry.insert(tk.END, str(peak1Min))

    peak1_max_label = ttk.Label(panel3, text='max (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak1_max_label.place(x=510, y=100)
    peak1_max_entry = ttk.Entry(panel3)
    peak1_max_entry.place(x=610, y=100, width=100)
    peak1_max_entry.insert(tk.END, str(peak1Max))

    peak1_name_label = ttk.Label(panel3, text='Name of the peak:')
    peak1_name_label.place(x=30, y=140)
    peak1_name_entry = ttk.Entry(panel3)
    peak1_name_entry.place(x=350, y=140, width=360)
    peak1_name_entry.insert(tk.END, 'G')

    # Peak 2
    peak2_label = ttk.Label(panel3, text='Peak 2', justify=tk.LEFT, font="bold")
    peak2_label.place(x=10, y=170)

    peak2_range_label = ttk.Label(panel3, text='Indicate spectral range:')
    peak2_range_label.place(x=30, y=200)

    peak2_min_label = ttk.Label(panel3, text='min (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak2_min_label.place(x=250, y=200)
    peak2_min_entry = ttk.Entry(panel3)
    peak2_min_entry.place(x=350, y=200, width=100)
    peak2_min_entry.insert(tk.END, str(peak2Min))

    peak2_max_label = ttk.Label(panel3, text='max (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak2_max_label.place(x=510, y=200)
    peak2_max_entry = ttk.Entry(panel3)
    peak2_max_entry.place(x=610, y=200, width=100)
    peak2_max_entry.insert(tk.END, str(peak2Max))

    peak2_name_label = ttk.Label(panel3, text='Name of the peak:')
    peak2_name_label.place(x=30, y=240)
    peak2_name_entry = ttk.Entry(panel3)
    peak2_name_entry.place(x=350, y=240, width=360)
    peak2_name_entry.insert(tk.END, 'D')

    # Peak 3
    peak3_label = ttk.Label(panel3, text='Peak 3', justify=tk.LEFT, font="bold")
    peak3_label.place(x=10, y=270)

    peak3_range_label = ttk.Label(panel3, text='Indicate spectral range:')
    peak3_range_label.place(x=30, y=300)

    peak3_min_label = ttk.Label(panel3, text='min (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak3_min_label.place(x=250, y=300)
    peak3_min_entry = ttk.Entry(panel3)
    peak3_min_entry.place(x=350, y=300, width=100)
    peak3_min_entry.insert(tk.END, str(peak3Min))

    peak3_max_label = ttk.Label(panel3, text='max (cm\u207B\u00B9)', justify=tk.RIGHT)
    peak3_max_label.place(x=510, y=300)
    peak3_max_entry = ttk.Entry(panel3)
    peak3_max_entry.place(x=610, y=300, width=100)
    peak3_max_entry.insert(tk.END, str(peak3Max))

    peak3_name_label = ttk.Label(panel3, text='Name of the peak:')
    peak3_name_label.place(x=30, y=340)
    peak3_name_entry = ttk.Entry(panel3)
    peak3_name_entry.place(x=350, y=340, width=360)
    peak3_name_entry.insert(tk.END, '2D')

   ################
    ### OK button ###
    ################
    values = {}  # Dictionary to store the values

    def ok_button_callback():
        # Store the values in the dictionary
        values['lorentz'] = method_var.get() == 'Method 2 (Lorentz Fit)'
        values['nt'] = split_peak_var.get()
        values['rbm'] = rbm_var.get()

        values['map'] = mapping_var.get()
        if values['map'] == 1:
            values['cols'] = int(x_pixels_entry.get())
            values['rows'] = int(y_pixels_entry.get())
            values['use_leng'] = length_var.get()
            if values['use_leng']==1:
                values['colmum'] = float(x_microns_entry.get())
                values['rowsmum'] = float(y_microns_entry.get())
            else:
                values['colmum'] = 0
                values['rowsmum'] = 0
        else:
            values['cols'] = 0
            values['rows'] = 0
            values['colmum'] = 0
            values['rowsmum'] = 0
            values['use_leng'] = 0

        values['normLow'] = float(norm_min_entry.get())
        values['normHigh'] = float(norm_max_entry.get())

        values['band1Low'] = float(peak1_min_entry.get())
        values['band1High'] = float(peak1_max_entry.get())
        values['band1Name'] = peak1_name_entry.get()

        values['band2Low'] = float(peak2_min_entry.get())
        values['band2High'] =float(peak2_max_entry.get())
        values['band2Name'] = peak2_name_entry.get()

        values['band3Low'] = float(peak3_min_entry.get())
        values['band3High'] = float(peak3_max_entry.get())
        values['band3Name'] = peak3_name_entry.get()

        # Destroy the main window
        window.destroy()

    ok_button = ttk.Button(window, text='OK', command=ok_button_callback)
    ok_button.place(x=740, y=530, width=40, height=35)

    window.mainloop()
    
    # Return the values
    return (
        values['lorentz'], values['nt'], values['rbm'], values['map'], values['cols'], values['rows'],
        values['colmum'], values['rowsmum'], values['normLow'], values['normHigh'], values['band1Low'],
        values['band1High'], values['band1Name'], values['band2Low'], values['band2High'], values['band2Name'],
        values['band3Low'], values['band3High'], values['band3Name'], values['use_leng']
    )

    