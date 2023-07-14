# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:46:15 2023

@author: matte
"""

import tkinter as tk
import numpy as np
import re
import ast

np.set_printoptions(suppress=True)

def str2array(s):
    # Remove space after [
    s=re.sub('\[ +', '[', s.strip())
    # Replace commas and spaces
    s=re.sub('[,\s]+', ', ', s)
    return np.array(ast.literal_eval(s))

def LorentzInputs(band1Name, band2Name, band3Name, total):
    Init_peak1plus = np.tile(np.array([1590, 30, 1]),(total,1))
    Init_peak2 = np.tile(np.array([1305, 40, 0.06]),(total,1))
    Init_peak3 = np.tile(np.array([2640, 50, 0.3]),(total,1))

    def ok_button():
        nonlocal Init_peak1plus, Init_peak2, Init_peak3
        Init_peak1plus = str2array(p5.get())
        Init_peak2 = str2array(p7.get())
        Init_peak3 = str2array(p9.get())
        fig.destroy()

    fig = tk.Tk()
    fig.geometry("1130x450")
    fig.title("Initial estimates for Lorentzian fitting")

    p = tk.LabelFrame(fig, width=1100, height=400,text="Initial estimates for Lorentzian fitting")
    p.place(x=17, y=16)

    p1 = tk.Label(p, text=f"- Initial estimates of spectral features of peaks 1, 2 and 3 should be provided for the {total} sample(s).",
                  justify="left", anchor="w",width=85)
    p1.place(x=9, y=30)

    p2 = tk.Label(p, text="- Format: 3 values per sample, separated by spaces. Values for different samples separated by ;",
                  justify="left", anchor="w",width=80)
    p2.place(x=9, y=60)

    p3 = tk.Label(p, text="[[center(1), FWHM(1), Max_intensity(1)], [center(2), FWHM(2), Max_intensity(2)], ...]",
                  justify="left", anchor="w",width=80)
    p3.place(x=102, y=90)

    p4 = tk.Label(p, text=f"Peak: {band1Name}", justify="left", anchor="w", font="bold", width = 78)
    p4.place(x=18, y=130)

    p5 = tk.Entry(p,width = 100)
    p5.place(x=40, y=170)
    p5.insert(tk.END, re.sub('[,\s]+', ', ', str(Init_peak1plus).replace('\n',',')))

    p6 = tk.Label(p, text=f"Peak: {band2Name}", justify="left", anchor="w", font="bold", width = 78)
    p6.place(x=18, y=210)

    p7 = tk.Entry(p,width = 100)
    p7.place(x=40, y=250)
    p7.insert(tk.END, re.sub('[,\s]+', ', ', str(Init_peak2).replace('\n',',')))

    p8 = tk.Label(p, text=f"Peak: {band3Name}", justify="left", anchor="w", font="bold", width = 78)
    p8.place(x=18, y=290)

    p9 = tk.Entry(p,width = 100)
    p9.place(x=40, y=330)
    p9.insert(tk.END, re.sub('[,\s]+', ', ', str(Init_peak3).replace('\n',',')))

    b = tk.Button(fig, text="OK", width=10, font=("TkDefaultFont", 10, "bold"),
                  fg="white", bg="gray50", command=ok_button)
    b.place(x=900, y=100)

    fig.mainloop()

    return Init_peak1plus, Init_peak2, Init_peak3

