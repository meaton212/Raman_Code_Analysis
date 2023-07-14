# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 12:28:35 2023

@author: matte
"""

import tkinter as tk
from tkinter import ttk
import numpy as np

def RBMinputs(n):
    RBMregion_Low = 200
    RBMregion_High = 400
    Prom = 0.01*np.ones(n)

    def ok_button():
        nonlocal RBMregion_Low, RBMregion_High, Prom
        RBMregion_Low = float(p4.get())
        RBMregion_High = float(p6.get())
        Prom = np.array(eval(p8.get()))
        fig.destroy()

    fig = tk.Tk()
    fig.geometry("850x450")
    fig.title("Radial breathing modes")

    p = tk.LabelFrame(fig, width=800, height=400, text="Radial breathing modes")
    p.place(x=17, y=23)

    p1 = tk.Label(p, text="Spectral Range: RBMs are the local intensity maxima within this range",
                  justify="left", anchor="w", width=60)
    p1.config(font=("TkDefaultFont", 10, "bold"))
    p1.place(x=15, y=30)

    p3 = tk.Label(p, text="min (cm\u207B\u00B9)", justify="left", anchor="w")
    p3.place(x=51, y=85)
    p4 = tk.Entry(p)
    p4.place(x=144, y=85)
    p4.insert(tk.END, str(RBMregion_Low))

    p5 = tk.Label(p, text="max (cm\u207B\u00B9)", justify="left", anchor="w")
    p5.place(x=415, y=85)
    p6 = tk.Entry(p)
    p6.place(x=515, y=85)
    p6.insert(tk.END, str(RBMregion_High))

    p7 = tk.Label(p, text="Prominence: minimum intensity (normalized) at which peaks will be considered.",
                  justify="left", anchor="w",  wraplength=800)
    p7.config(font=("TkDefaultFont", 10, "bold"))
    p7.place(x=15, y=150)

    p71 = tk.Label(p, text="(Local maxima with a maximum height lower than the prominence will be discarded)",
                    justify="left", anchor="w", wraplength=800, font=("TkDefaultFont", 9))
    p71.place(x=50, y=180)

    p8 = tk.Entry(p, width=60)
    p8.place(x=100, y=220)
    p8.insert(tk.END, ",".join(np.char.mod('%.2f', Prom)))
    
    p81 = tk.Label(p, text="Syntax: comma separated values: Prom(sample1),Prom(sample2),Prom(sample3), etc...",
                    justify="left", anchor="w", wraplength=800, font=("TkDefaultFont", 9))
    p81.place(x=50, y=260)

    b = tk.Button(fig, text="OK", width=10, font=("TkDefaultFont", 10, "bold"),
                  fg="white", bg="gray50", command=ok_button)
    b.place(x=350, y=360)

    fig.mainloop()

    return RBMregion_Low, RBMregion_High, Prom

