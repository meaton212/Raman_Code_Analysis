# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:28:15 2023

@author: matte
"""

# load packages
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
# Import the data_processing module
import raman_data_processing


#%% GUI

# Input data:
# The program expects multiple (up to 10) files in .txt or .dat format. In each file, the Raman shift
# #should be included in the first column, followed by the Intensity data in
# consecutive columns (as many spectra available, and not neccesarily the same number of spectra per file).

# The location and names of the files should be indicated in path and file_name
# variables, respectively. The total number of files should be included in
#  total.

def open_files():
    global selected_files
    files = filedialog.askopenfilenames(filetypes=my_filetypes)
    if len(files) + len(selected_files) > 10:
        messagebox.showerror(title="Error", message="Exceeded maximum file selection limit (10 files).")
        return
    selected_files.extend(files)
    for file in files:
        file_listbox.insert(tk.END, file)

def remove_selected():
    selected_indices = file_listbox.curselection()
    for index in reversed(selected_indices):
        file_listbox.delete(index)
        del selected_files[index]

def process_files():
    global selected_files, file_name, label_window
    if selected_files:
        # Create a new window for label entry
        label_window = tk.Toplevel(root)
        label_window.title("Label Entry")
        len1=175+len(selected_files)*60
        label_window.geometry("600x"+str(len1))

        label_title = tk.Label(label_window, text="Input the label for each file to be used in the figures", font=("Arial", 10, "bold"))
        label_title.pack(pady=10)

        label_frame = tk.Frame(label_window)
        label_frame.pack(pady=10)

        label_entries = []

        for file in selected_files:
            filename = file.split('/')[-1]
            label_label = tk.Label(label_frame, text=filename+': ')
            label_label.pack()

            label_entry = tk.Entry(label_frame, width=50)
            label_entry.pack()
            label_entries.append(label_entry)

        def save_labels():
            global labels
            labels = [entry.get() for entry in label_entries]
            for file, label in zip(selected_files, labels):
                file_name.append(file.split('/')[-1])
                print(f"File: {file}\tLabel: {label}")
            label_window.destroy()
            destroy_all_windows()
            
            # Call a function from data_processing.py and pass the variables
            raman_data_processing.process_data(folder_selected, selected_files, file_name, labels,delim)
            

        save_button = tk.Button(label_window, text="Save Labels", command=save_labels)
        save_button.pack(pady=10)

    else:
        print("No files selected.")

        destroy_all_windows()
        

def destroy_all_windows():
    for window in root.winfo_children():
        window.destroy()
    root.destroy()

root = tk.Tk()
root.withdraw()

messagebox.showinfo(title="Greetings", message='In the next window, select the path where the files to analyze are located')

folder_selected  = filedialog.askdirectory(initialdir='C:/', title="Please select folder with files to analyze:", parent=root)

os.chdir(folder_selected)

root.deiconify()
root.title('Data information')
root.geometry("800x200")

ftype = '.txt'


def submit():
    global ftype, my_filetypes, dtype, delim
    ftype = type_var.get()
    if ftype == '.txt':
        my_filetypes = [('Text files', '*.txt'),('all files', '.*')]
    elif ftype == '.csv':
        my_filetypes = [('CSV files', '*.csv'),('all files', '.*')]
    elif ftype == '.dpt':
        my_filetypes = [('DPT files', '*.dpt'),('all files', '.*')]
    elif ftype == '.dat':
        my_filetypes = [('DAT files', '*.dat'),('all files', '.*')]
    elif ftype == 'other':
        my_filetypes =  [('all files', '.*')]
    
    dtype=delimiter_var.get()
    if dtype=='comma (,)':
            delim=','
    elif dtype=='tab':
            delim='\t'
    elif dtype=='space':
            delim='\s+'
    elif dtype== 'semi-colon (;)':
            delim=';'
    elif dtype=='slash (/)':
            delim='/'
        
    root.destroy()

label2 = tk.Label(root, text='Please Select File Type (.txt, .csv, .dpt, or .dat):', font=("Arial", 10, "bold"))
label2.pack()
type_var = tk.StringVar(root)
type_var.set('.txt')  # Default value
type_dropdown = ttk.Combobox(root, textvariable=type_var, values=['.txt', '.csv', '.dpt', '.dat', 'other'], width=10)
type_dropdown.pack()

delimiter_label = tk.Label(root, text='Select Delimiter:', font=("Arial", 10, "bold"))
delimiter_label.pack()
delimiter_var = tk.StringVar(root)
delimiter_var.set('tab')  # Default value
delimiter_dropdown = ttk.Combobox(root, textvariable=delimiter_var, values=['tab','comma (,)', 'space', 'semi-colon (;)', 'slash (/)'], width=15)
delimiter_dropdown.pack()

submit_button = tk.Button(root, text='Submit', command=submit)
submit_button.pack()

root.mainloop()


selected_files = []
file_name = []

root = tk.Tk()
root.title("Multiple File Selection")
root.geometry("1250x650")

file_frame = tk.Frame(root)
file_frame.pack(pady=20)

scrollbar = tk.Scrollbar(file_frame)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

file_listbox = tk.Listbox(file_frame, selectmode=tk.MULTIPLE, height=20, width=125)
file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

scrollbar.config(command=file_listbox.yview)
file_listbox.config(yscrollcommand=scrollbar.set)

button_frame = tk.Frame(root)
button_frame.pack(pady=10)

open_button = tk.Button(button_frame, text="Open Files", command=open_files)
open_button.pack(side=tk.LEFT, padx=5)

remove_button = tk.Button(button_frame, text="Remove Selected", command=remove_selected)
remove_button.pack(side=tk.LEFT, padx=5)

process_button = tk.Button(button_frame, text="Continue", command=process_files)
process_button.pack(side=tk.LEFT, padx=5)

root.mainloop()