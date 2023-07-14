# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 17:48:10 2023

@author: matte
"""





import tkinter as tk
from tkinter import ttk

def PlotsInputs(band1Name, band2Name, band3Name):
    width = 3
    width_fw = 5
    width_int = 0.1
    raw = False
    norm = False
    rng = False
    peaks = False
    correlation1 = False
    correlation2 = False
    correlation3 = False
    correlation4 = False
    correlation5 = False
    map1 = False
    map2 = False
    map3 = False
    map4 = False
    dens = False
    saveFigs = True
    nameFigs = "TestSample"
    imgtype='.svg'
    fs=16
    maxI21=100

    def finish_button():
        nonlocal width, width_fw, width_int, raw, norm, rng, peaks, correlation1, correlation2, correlation3, correlation4, correlation5, map1, map2, map3, map4, saveFigs, nameFigs, imgtype, dens, fs,maxI21
        width = float(w3.get())
        width_fw = float(w5.get())
        width_int = float(w7.get())
        raw = o1_var.get()
        norm = o2_var.get()
        rng = o3_var.get()
        peaks = o4_var.get()
        correlation1 = q1_var.get()
        correlation2 = q2_var.get()
        correlation3 = q3_var.get()
        correlation4 = q4_var.get()
        correlation5 = q5_var.get()
        map1 = r1_var.get()
        map2 = r2_var.get()
        map3 = r3_var.get()
        map4 = r4_var.get()
        saveFigs = f1_var.get()
        nameFigs = f3.get()
        imgtype=imgtype_var.get()
        dens=w8_var.get()
        fs=int(ff3.get())
        maxI21=float(ff5.get())
        fig.destroy()

    fig = tk.Tk()
    fig.geometry("955x630")
    fig.title("Plotting and saving options")

    p1 = tk.LabelFrame(fig, width=925, height=600, text="Plotting and saving options",font="bold")
    p1.place(x=15, y=10)

    p1a = tk.Label(fig, text="Figures plotted by default:", justify="left", anchor="w")
    p1a.place(x=30, y=40)

    p1b = tk.Label(fig, text="- Histograms of all spectral features", justify="left", anchor="w")
    p1b.place(x=200, y=80)

    p1c = tk.Label(fig, text=f"- Maps (if selected) of position and I({band2Name})/I({band1Name})", justify="left", anchor="w")
    p1c.place(x=200, y=120)

    p1d = tk.Label(fig, text="Select the additional graphs required:", justify="left", anchor="w")
    p1d.place(x=30, y=225)
    
    
    ff1 = tk.LabelFrame(fig, width=200, height=110, text="Plot settings")
    ff1.place(x=450, y=160)

    ff2 = tk.Label(ff1, text="font size:")
    ff2.place(x=10, y=5)

    ff3 = tk.Entry(ff1, width=5)
    ff3.place(x=125, y=5)
    ff3.insert(tk.END, str('16'))
    
    ff4 = tk.Label(ff1, text="Max I("+band2Name+")/I("+band1Name+"):")
    ff4.place(x=10, y=45)

    ff5 = tk.Entry(ff1, width=5)
    ff5.place(x=125, y=45)
    ff5.insert(tk.END, str(maxI21))

    w1 = tk.LabelFrame(fig, width=250, height=220, text="Width of histograms")
    w1.place(x=670, y=50)

    w2 = tk.Label(w1, text="Raman shift")
    w2.place(x=30, y=20)

    w3 = tk.Entry(w1, width=5)
    w3.place(x=150, y=20)
    w3.insert(tk.END, str(width))

    w4 = tk.Label(w1, text="FWHM")
    w4.place(x=65, y=60)

    w5 = tk.Entry(w1, width= 5)
    w5.place(x=150, y=60)
    w5.insert(tk.END, str(width_fw))

    w6 = tk.Label(w1, text="Intensity Ratio")
    w6.place(x=11, y=100)

    w7 = tk.Entry(w1, width=5)
    w7.place(x=150, y=100)
    w7.insert(tk.END, str(width_int))
    
    w8_var = tk.BooleanVar(value=False)
    w8 = tk.Checkbutton(w1, text="Normalize Histograms", variable=w8_var)
    w8.place(x=11, y=140)

    o = tk.LabelFrame(fig, width=300, height=200, text="Figures: general graphs")
    o.place(x=30, y=280)

    o1_var = tk.BooleanVar(value=False)
    o1 = tk.Checkbutton(o, text="Raw spectra", variable=o1_var)
    o1.place(x=18, y=8)

    o2_var = tk.BooleanVar(value=False)
    o2 = tk.Checkbutton(o, text="Normalised spectra", variable=o2_var)
    o2.place(x=18, y=37)

    o3_var = tk.BooleanVar(value=False)
    o3 = tk.Checkbutton(o, text="Spectral regions for peaks", variable=o3_var)
    o3.place(x=18, y=66)

    o4_var = tk.BooleanVar(value=False)
    o4 = tk.Checkbutton(o, text="RBM region peak identification", variable=o4_var)
    o4.place(x=18, y=95)

    q = tk.LabelFrame(fig, width=350, height=200, text="Figures: spectral features and correlations")
    q.place(x=350, y=280)

    q1_var = tk.BooleanVar(value=False)
    q1 = tk.Checkbutton(q, text=f"position {band1Name} vs I({band2Name})/I({band1Name})", variable=q1_var)
    q1.place(x=9, y=8)

    q2_var = tk.BooleanVar(value=False)
    q2 = tk.Checkbutton(q, text=f"position {band2Name} vs I({band2Name})/I({band1Name})", variable=q2_var)
    q2.place(x=9, y=38)

    q3_var = tk.BooleanVar(value=False)
    q3 = tk.Checkbutton(q, text=f"position {band3Name} vs I({band2Name})/I({band1Name})", variable=q3_var)
    q3.place(x=9, y=68)

    q4_var = tk.BooleanVar(value=False)
    q4 = tk.Checkbutton(q, text=f"position {band3Name} vs position {band1Name}", variable=q4_var)
    q4.place(x=9, y=98)

    q5_var = tk.BooleanVar(value=False)
    q5 = tk.Checkbutton(q, text="Position vs FWHM (all peaks)", variable=q5_var)
    q5.place(x=9, y=128)

    r = tk.LabelFrame(fig, width=210, height=200, text="Figures: Additional Maps")
    r.place(x=720, y=280)


    r1_var = tk.BooleanVar(value=False)
    r1 = tk.Checkbutton(r, text=f"FWHM {band1Name}", variable=r1_var)
    r1.place(x=9, y=8)

    r2_var = tk.BooleanVar(value=False)
    r2 = tk.Checkbutton(r, text=f"FWHM {band2Name}", variable=r2_var)
    r2.place(x=9, y=37)

    r3_var = tk.BooleanVar(value=False)
    r3 = tk.Checkbutton(r, text=f"FWHM {band3Name}", variable=r3_var)
    r3.place(x=9, y=66)

    r4_var = tk.BooleanVar(value=False)
    r4 = tk.Checkbutton(r, text=f"I({band3Name})/I({band1Name})", variable=r4_var)
    r4.place(x=9, y=95)

    f1_var = tk.BooleanVar(value=True)
    f1 = tk.Checkbutton(fig, text="Save Figures", variable=f1_var)
    f1.place(x=50, y=490)
    
    
    img_lab = tk.Label(fig, text="Image type:")
    img_lab.place(x=220, y=493)
    
    #Dropdown menu to select saved image type
    imgtype_var = tk.StringVar()
    imgtype_var.set('.svg')  # Default value
    imgtype_dropdown = ttk.Combobox(fig, textvariable=imgtype_var, values=['.svg', '.png', '.tiff', '.jpg'], width=4)
    imgtype_dropdown.place(x=325, y=493)

    f2 = tk.Label(fig, text="Basename:")
    f2.place(x=500, y=493)

    f3 = tk.Entry(fig, width=31)
    f3.place(x=600, y=493)
    f3.insert(tk.END, "TestSample")

    b = tk.Button(fig, text="FINISH", width=10, font=("TkDefaultFont", 10, "bold"), fg="white", bg="gray50", command=finish_button)
    b.place(x=400, y=550)

    fig.mainloop()

    return width, width_fw, width_int, raw, norm, rng, peaks, correlation1, correlation2, correlation3, correlation4, correlation5, map1, map2, map3, map4, saveFigs, nameFigs, imgtype, dens, fs, maxI21


