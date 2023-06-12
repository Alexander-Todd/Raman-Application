# from tqdm import tqdm
# import time

# total_imports = 15  # Total number of imports

# print("Loading, please wait...")
# for _ in tqdm(range(total_imports), desc="Imports"):
#     # Simulating import delay
#     time.sleep(1)
from sys import stdout
stdout.write("Importing GUI Modules")
stdout.flush()
from tkinter import filedialog as fd
import tkinter as tk
from tkinter import ttk
stdout.write("\rGUI modules imported ")
stdout.flush()
stdout.write("\nImporting Graphing Modules")
stdout.flush()
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import csv
import pickle
stdout.write("\rGraphing modules imported ")
stdout.flush()
stdout.write("\nImporting MCR Modules")
stdout.flush()
# for MCR-ALS
import pymcr
from pymcr.mcr import McrAR
from pymcr.regressors import OLS, NNLS
from pymcr.constraints import ConstraintNonneg, ConstraintNorm
# from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
stdout.write("\rMCR modules imported ")
stdout.flush()
stdout.write("\nImporting Raman Processing methods")
stdout.flush()
# ProcessingApplication.py is in the same directory as RamanProcessing.py
# This should work as long as this is the case
import newRamanProcessing as rp
stdout.write("\rRaman Processing methods imported ")
stdout.flush()
stdout.write("\n")
print("Program initiated")


SMALL_SIZE = 15
MEDIUM_SIZE = 20
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Setting up the window
root = tk.Tk()
root.title("Raman processing application")
root.geometry('550x800')
root.resizable(True, True)

# Setting up the tabs
tabControl = ttk.Notebook(root)
tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)
tabControl.add(tab1, text='Preprocessing')
tabControl.add(tab2, text='Spectral decomposition')
tabControl.pack(expand=1, fill="both")

# Create the frame for file selection
frame_files = ttk.Frame(master=tab1)
# Contains the user-changable variables, can be hidden
frame_processing = ttk.Frame(master=tab1, borderwidth=4, relief=tk.RAISED)
# Create the frame for graph settings
frame_display = ttk.Frame(master=tab1)

# The frames inside the processing frame, one for each process
# This is to allow finding, deleting, and reinstalling
# the labels and entry boxes
smooth_frame = ttk.Frame(master=frame_processing)
normalise_frame = ttk.Frame(master=frame_processing)
baseline_frame = ttk.Frame(master=frame_processing)
crs_frame = ttk.Frame(master=frame_processing)

# Spectral decomposition frame
filesdecomp_frame = ttk.Frame(master=tab2)
decomp_frame = ttk.Frame(master=tab2)
decomp_method_frame = ttk.Frame(master=decomp_frame)

#-----------------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------------

# The command for the "select files" button
def select_files():
    # Specify the file types the user will be allowed to select
    # For some reason it needs the extra comma or it breaks...?
    filetypes = (('text files', '*.txt'),)

    # Open file browse window where files can be selected from
    # filepath is stored
    filenames = fd.askopenfilenames(
        title='Open files',
        initialdir='C:/Users/1alex/OneDrive/Documents/Y5 Work/Project/Code/Raman-Processing-main/EthanolData/TIRF_100s_5um_below_quartz_ethanol_comparison_thirds',
        filetypes=filetypes)
    # If cancel is clicked
    if filenames is None:
       return

    # For each filepath selected, insert into bottom of the listbox
    for filepath in filenames:
        file_list_display.insert("end", filepath)

def select_component(widget, clear=True):
   filetypes = (('text files', '*.txt'),)
   filenames = fd.askopenfilenames(
      title='Open files',
      initialdir='C:/Users/1alex/OneDrive/Documents/Y5 Work/Project/Code/Raman-Processing-main/EthanolData/TIRF_100s_5um_below_quartz_ethanol_comparison_thirds',
      filetypes=filetypes)
   # If cancel is clicked
   if filenames is None:
      return   
   
   for filepath in filenames:
      if clear==True:
        widget.delete(0, tk.END)
      widget.insert("end", filepath)

# The command for the "remove files" button
def remove_files():
    file_list_display.delete(file_list_display.curselection())

def removeAll_files():
   file_list_display.delete(0,'end')

# The function in the "Show/hide variables" button
def show_hide():
    # Forget if present; if not present then show on the grid
    frame_processing.grid_forget() if frame_processing.winfo_manager() else frame_processing.grid(row=1, column=0)

def labelled_true():
   if labelledDataVar.get():
      frame_legend.grid(row=5, column=0, columnspan=3)
      frame_decomp_graph.grid(row=3, column=0, columnspan=2)
   else:
      frame_legend.grid_forget()
      frame_decomp_graph.grid_forget()


def get_legend():
    filepaths = list(file_list_display.get(0, 'end')) # get file paths from the listbox
    separator = legendChar.get()
    if len(separator) == 1:
       pass
    else:
       print("""Error with separator character.
             There must be a single character specified to denote where the file label should be collected from.
             Assigning consecutive numbers in the meantime.""")
       separator = False
    
    if separator:
      sampleNames = get_sample_names(filepaths, separator)
    else:
       sampleNames = []
       for count, file in enumerate(filepaths):
          sampleNames.append(str(count+1))
    return sampleNames

def updateLegendVals():
   sampleNames = get_legend()
   collectedEntry.delete(0, tk.END)
   collectedEntry.insert(tk.END, sampleNames)
   

# Lambda function for changing labels when the package is changed
# only runs when dropdowns are changed manually
def change_through_package(package, process, option):
    if option == False:
       checkbox_list[process].set(False)
       return checkbox_list
    else:
      for widget in dict_frames[process].winfo_children():
        widget.destroy()       
      if type(option) == tk.StringVar:  # can we keep the type the same so this isn't needed?
        option = option.get()
      for count, parameter in enumerate(processing_package[package][process][option]):
        dropdown_list[process].set(option)
        lab = tk.Label(master=dict_frames[process], text=parameter)
        lab.grid(column=0, row=count)
  
        # some parameters are None, don't want these on the grid
        if type(processing_package[package][process][option][parameter]): 
          ent = tk.Entry(master=dict_frames[process])
          ent.insert(tk.END, processing_package[package][process][option][parameter])
          ent.grid(column=1, row=count)

decomp_storage = []
def change_thu_decomp(method, dictionary):
   for widget in decomp_method_frame.winfo_children():
      widget.destroy()
   count = 0
   for key, value in dictionary[method].items():
      lab = tk.Label(master=decomp_method_frame, text=key)
      lab.grid(column=0, row=count)
      if type(value) == bool:
         start = tk.BooleanVar()
         widg = ttk.Checkbutton(decomp_method_frame, variable=start)
         widg.state(['!alternate'])
         if value == True:
            widg.state(['selected'])
      elif type(value) == int:
         start = tk.IntVar(value=value)
         widg = tk.Entry(decomp_method_frame)
         widg.insert(tk.END, start.get())
      elif type(value) == list:
         start = tk.StringVar(value=decomp_defaults[method][key])
         widg = tk.OptionMenu(decomp_method_frame, start, *decomp_option_dict[method][key])
      widg.grid(column=1, row=count)
      decomp_storage.append(start)
      count += 1

def get_params_from_file():
   filetypes = (('text files', '*.txt'),)
   params_file = fd.askopenfilenames(
      title='Select parameters file',
      initialdir='C:/Users/1alex/OneDrive/Documents/Y5 Work/Project/Code/Raman-Processing-main',
      filetypes=filetypes)
   # If cancel is clicked
   if params_file is None:
      return

   with open(params_file[0], 'r') as file:
      processing_dict={}
      decomp_dict={}      
      processing = True      
      reader=csv.reader(file)
      next(reader)  # Skip the header
      for row in reader:
         try:
            key, value = row
            if processing == True:
                if key != "Decomposition":
                  processing_dict[key] = value
                else:
                  processing = False
                  decomp_method = value
            else:
                try:
                   value = int(value)
                except ValueError:
                   pass
                decomp_dict[key] = value
         # ValueError when line is blank
         except ValueError:
            pass

      decomp_dict[decomp_method] = decomp_dict
      processing_package["Custom"] = processing_dict
      print(decomp_dict)
      print(processing_package)
      package_change("Custom")
      change_thu_decomp(decomp_method, decomp_dict)

# The function called by the lambda function of a dropdown; runs when changed manually
def change_manually(option, process):
    for widget in dict_frames[process].winfo_children():
      widget.destroy()
    for count, parameter in enumerate(process_menu[process][option]):
      lab = tk.Label(master=dict_frames[process], text=parameter)
      lab.grid(column=0, row=count)
      ent = tk.Entry(master=dict_frames[process])
      ent.insert(tk.END, process_menu[process][option][parameter])
      ent.grid(column=1, row=count)

# The function called by the lambda function of the package dropdown,
# is ran when the dropdown is changed
def package_change(package):
  package_option.set(package)
  for process in processing_package[package]:
    try:
      for option in processing_package[package][process]:
        change_through_package(package, process, option)
    except TypeError:
       checkbox_list[process].set(False)

# function to be called when the "plot graph" button is clicked
def plot_graph():
  file_paths = list(file_list_display.get(0, 'end')) # get file paths from the listbox
  graph_parameters['title'] = titleEntry.get()
  graph_parameters['x-axis'] = xEntry.get()
  graph_parameters['y-axis'] = yEntry.get()
  graph_parameters['show y=0'] = showaxis.get()
  graph_parameters['Show legend'] = showLeg.get()
  graph_parameters['figSizeX']= figureDimsEntryX.get()
  graph_parameters['figSizeY'] = figureDimsEntryY.get()
  graph_parameters['titleFontSize'] = titleFontEntry.get()
  graph_parameters['axesFontSize'] = axesTitleEntry.get()
  graph_parameters['axesValuesSize'] = axesValuesEntry.get()
  graph_parameters['legendFontSize'] = legendFontEntry.get()

  sample_names = get_legend()

  thepackage = get_processing_params()
  df, wn, intensities = rp.quickProcess(file_paths, sample_names, thepackage)
  # Set size of window
  try:
     xSize = int(graph_parameters["figSizeX"])
     ySize = int(graph_parameters["figSizeY"])
     plt.rcParams["figure.figsize"] = [xSize, ySize]
  except ValueError:
     plt.rcParams["figure.figsize"] = [15,9]
  # Set font size of the tick labels
  try:
    axesValuesSize = int(graph_parameters["axesValuesSize"])
    plt.rc('xtick', labelsize=axesValuesSize)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=axesValuesSize)    # fontsize of the tick labels
  except:
     plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
  # Set size of legend font  
  try:
    legendFont = int(graph_parameters["legendFontSize"])
    plt.rc('legend', fontsize=legendFont)    # legend fontsize
  except:
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

  ax = plt.gca()
  ax.xaxis.set_minor_locator(MultipleLocator(100))
  ax.tick_params(width=1.5)
  for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)  # change width
  if graph_parameters["show y=0"]:
    plt.axhline(y=0,color='0.2')
  plt.plot(np.transpose(wn), np.transpose(intensities), label=sample_names)
  plt.autoscale(enable=True, axis='x', tight=True)
  # Set title font
  try:
     titleFont = int(graph_parameters["titleFontSize"])
     plt.title(graph_parameters['title'], fontsize=titleFont)
  except ValueError:
     plt.title(graph_parameters['title'], fontsize=BIGGER_SIZE)
  # Set font size of axes titles
  try:
    axesTitles = int(graph_parameters["axesFontSize"])
    plt.xlabel('Wavenumbers (cm$^{-1}$)', fontsize=axesTitles)
    plt.ylabel('Intensity (a.u.)', fontsize=axesTitles)
  except ValueError:
    plt.xlabel('Wavenumbers (cm$^{-1}$)', fontsize=MEDIUM_SIZE)
    plt.ylabel('Intensity (a.u.)', fontsize=MEDIUM_SIZE)
  # Show legend if selected
  if graph_parameters['Show legend']:  # if true
    plt.legend()
  plt.show()

def get_custom():
   for frame in frame_processing.winfo_children():
      print(frame.winfo_class())
      if frame.winfo_class() == 'Frame':
         for widget in frame.winfo_children():
          try:
            print(widget, " = ", widget.get())
          except:
            print(widget, " says nope")
      else:
         pass

def get_decomp():
   fetched_pairs = {}
   for widget in decomp_method_frame.winfo_children():
      if widget.winfo_class() == 'Label':
         lab = widget.cget('text')
      else:
        if widget.winfo_class() == "Entry":
          val = widget.get()
        elif widget.winfo_class() == "TCheckbutton":
          val = widget.instate(['selected'])
        elif widget.winfo_class() == "Menubutton":
          val = widget.cget("text")
        pair = {lab: val}
        fetched_pairs.update(pair)
      
   method = decomp_option_menu.cget("text")
   fetched_sum = {method : fetched_pairs}
   print(fetched_sum)
   return fetched_sum
        
def multivariate_resolution(unlabelled, labelled, method_dict):
    constraints = []
    if method_dict["Sum to 1 constraint"] == True:
       constraints.append(ConstraintNorm())
    if method_dict["Non-negative constraint"] == True:
       constraints.append(ConstraintNonneg())
    mcrar = McrAR(max_iter=method_dict["Maximum iterations"],
                  c_regr=method_dict["Concentration regressor"],
                  tol_err_change=method_dict["Tolerance error change"],
                  tol_increase=method_dict["Tolerance increase"]
                  )
    if method_dict["Spectral regressor"] == "None":
       mcrar.fit(unlabelled, ST=labelled, st_fix=[0,1], verbose=True)
    else:
       mcrar = McrAR(st_regr=method_dict["Spectral regressor"])
       mcrar.fit(unlabelled, ST=labelled, verbose=True)
    return mcrar

def nonneg_elastic_net(unlabelled, labelled, method_dict):
   params = list(method_dict.values)
   # Not implemented

conc_results = None

def calculate_elements():
   ''' Decomposing spectra using two reference points.
   Runs when user clicks "run" button on decomposition tab.
   '''
   # get the method from dropdown of decompositions
   # This is given when button is pressed? No
   method = str(decomp_dropdown_default.get())

   mixture_filepaths = list(file_list_display.get(0, 'end'))
   labelled_filepaths = [str(a_file_entry.get()), str(b_file_entry.get())]
   mixture_samplenames = get_sample_names(mixture_filepaths)
   lab_samplenames = get_sample_names(labelled_filepaths)

   thepackage = str(package_option.get())
   the_dict = get_processing_params()
   lab_df, lab_wn, lab_int = rp.quickProcess(labelled_filepaths, lab_samplenames, the_dict)
   mix_df, mix_wn, mix_int = rp.quickProcess(mixture_filepaths, mixture_samplenames, the_dict)
   lab_int[1] = lab_int[1] - lab_int[0]  # Subtract 0% from 10%

   parameters = []
   for param in decomp_storage:
      parameters.append(param.get())

   # params = {"Maximum iterations": 50,
   #             "Spectral regressor": "OLS",
   #             "Concentration regressor": "OLS",
   #             "Tolerance error change": 1e-12,
   #             "Tolerance increase": 2,
   #             "Sum to 1 constraint": False,
   #             "Non-negative constraint": False}
   params = {"Maximum iterations": parameters[0],
            "Spectral regressor": parameters[1],
            "Concentration regressor": parameters[2],
            "Tolerance error change": 1e-12,
            "Tolerance increase": 2,
            "Sum to 1 constraint": parameters[3],
            "Non-negative constraint": False}

   global conc_results
   conc_results = decompositions[method](mix_int, lab_int, params)
      
   a_concs = [float(i) for i in conc_results.C_opt_[:,0]]
   b_concs = [float(i) for i in conc_results.C_opt_[:,1]]
   reference_quantity = float(reference_quantity_entry.get())
   # conc_estimate = []
   # for i in b_concs:
   #    conc_estimate.append(i*reference_quantity)
   conc_estimate = [np.array(i) * reference_quantity for i in b_concs]

   for count, result in enumerate(b_concs):
      print(mixture_samplenames[count], ":", result*reference_quantity)

   # If the data is labelled
   if labelledDataVar.get():
      coefficients = np.polyfit(mixture_samplenames, conc_estimate, 1)
      # If output stats = True
      if stats_var.get():
         mse = mean_squared_error(mixture_samplenames, conc_estimate)
         print("RMSE:", np.sqrt(mse))
         R_squared = r2_score(mixture_samplenames, conc_estimate)
         print("R2:", R_squared)
         print("PRE (%):", (1-R_squared)*100)
         print("Equation:", coefficients[0], "x + ", coefficients[1])
      

      x_of_fit = np.linspace(0, 10, 1000)
      y_of_fit = [coefficients[1] + coefficients[0]*x for x in x_of_fit]

      # plt.rcParams['figure.figsize'] = [18,10]
      # plt.axhline(y=0,color='0.2')
      # plt.scatter(mixture_samplenames, a_concs, label='Component A Estimates')
      if xy_var.get():
         endpoints = [0,reference_quantity]
         plt.plot(endpoints, endpoints, label="Perfect correlation", linestyle='dotted')
      plt.scatter(mixture_samplenames, conc_estimate, label='Concentration Estimates')
      if least_sq_var.get():
         plt.plot(x_of_fit, y_of_fit, label="Least squares fit")
      plt.autoscale(enable=True, axis='x', tight=True)
      the_units = str(units_var.get())
      plt.xlabel("Actual Concentration (%s)" % the_units)
      plt.ylabel("Estimated Concentration (%s)" % the_units)
      print(scale_var.get())
      if scale_var.get() == "Log":
         plt.xscale("log")
         plt.yscale("log")
      plt.legend()
      plt.tight_layout()
      plt.show()
   


def save_parameters():
  preprocessing = get_processing_params()
  # Get the decomposition parameters
  fetched_pairs = {}
  for widget in decomp_method_frame.winfo_children():
    if widget.winfo_class() == 'Label':
        lab = widget.cget('text')
    else:
      if widget.winfo_class() == "Entry":
        val = widget.get()
      elif widget.winfo_class() == "TCheckbutton":
        val = widget.instate(['selected'])
      elif widget.winfo_class() == "Menubutton":
        val = widget.cget("text")
      pair = {lab: val}
      fetched_pairs.update(pair)  
  decomp_method = decomp_option_menu.cget("text")

  f = fd.asksaveasfile(mode='w', defaultextension=".txt")
  if f is None:  # If cancel is clicked
     return

  # Get the fieldnames for the DictWriter from the nested dictionary
  fieldnames = ['Process', 'Parameters']

  # Create a DictWriter object and write the header row
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()

  # Write each person's data as a row to the CSV file
  for key, value in preprocessing.items():
      print(value)
      row = {'Process': key, 'Parameters': value}
      writer.writerow(row)
  row = {'Process': "Decomposition", 'Parameters': decomp_method}
  writer.writerow(row)
  for key, value in fetched_pairs.items():
      row = {'Process': key, 'Parameters': value}
      writer.writerow(row)
  f.close()

def save_preprocessed():
  file_paths = list(file_list_display.get(0, 'end')) # get file paths from the listbox
  sample_names = get_legend()
  the_package = get_processing_params()
  df, wn, intensities = rp.quickProcess(file_paths, sample_names, the_package)

  # If the spectra are normalised then values can be in order of 1e-4 which gets exported as 0.

  arraySave = []
  for count, array in enumerate(wn):
     arraySave.append(array)
     arraySave.append(intensities[count])

  arraySave = np.stack(np.array(arraySave), axis=-1)

  f = fd.asksaveasfile(defaultextension=".csv")
  if f is None:  # If cancel is clicked
     return
  
  np.savetxt(f, arraySave, delimiter=',', fmt='%.4e')


def pickle_save():
  preprocessing = get_processing_params()
  # Get the decomposition parameters
  fetched_pairs = {}
  for widget in decomp_method_frame.winfo_children():
    if widget.winfo_class() == 'Label':
        lab = widget.cget('text')
    else:
      if widget.winfo_class() == "Entry":
        val = widget.get()
      elif widget.winfo_class() == "TCheckbutton":
        val = widget.instate(['selected'])
      elif widget.winfo_class() == "Menubutton":
        val = widget.cget("text")
      pair = {lab: val}
      fetched_pairs.update(pair)  
  decomp_method = decomp_option_menu.cget("text")

  f = fd.asksaveasfile(mode='w', defaultextension=".txt")
  if f is None:  # If cancel is clicked
     return   
  
  with open(f, 'wb') as handle:
     pickle.dump()

def save_concentrations():
   filepaths = list(file_list_display.get(0, 'end'))
   f = fd.asksaveasfile(mode='w', defaultextension=".txt")
   if f is None:  # If cancel is clicked
      return
   headerlist = ["Water + Substrate Component", "Analyte Component", "Filepath"]
   writer = csv.writer(f)
   writer.writerow(headerlist)
   for count, file in enumerate(filepaths):
      writer.writerow([conc_results.C_opt_[count, 0], conc_results.C_opt_[count, 1], file])
   f.close()

def get_sample_names(filepaths, separator="%"):
   storage = []
   for filepath in filepaths:
       tail = filepath
       while True:
           head, sep, tail = tail.partition("\\")
           if len(sep) > 0:
               pass
           else:
               storage.append(head)
               break
        
   sample_names = []
   for sample_name in storage:
       head, tail, sep = sample_name.partition(separator)
       while len(head) > 0:
           try:
               head = float(head)
               sample_names.append(abs(head))  # Can't have -ve concs!
               break
           except ValueError:
               head = head[1:]
   return sample_names

def get_processing_params():
  fetched_sum = {}
  for process, frame in dict_frames.items():
    fetched_pairs = {}
    if checkbox_list[process].get() == True:
      for widget in frame.winfo_children():
          if widget.winfo_class() =='Label':
            lab = widget.cget("text")
          if widget.winfo_class() == 'Entry':
            ent = int(widget.get())
            pair = {lab: ent}
            fetched_pairs.update(pair)
      fetched_method = {dropdown_list[process].get(): fetched_pairs}
    else:
       fetched_method = False
    fetched_sum.update({process: fetched_method})

  return(fetched_sum)
       
#-----------------------------------------------------------------------------------
# Dictionary Creation
#-----------------------------------------------------------------------------------

# Dictionary of the smoothing methods, parameters, and default values
smoothing_options = {
  "Savitzky–Golay": {"Window": 3, "Polynomial": 0},
  "FFT": {"Fourier values":200},
  "Wavelet": {"Wavelet": "db29", "Wavelet level": 1},
  "Whittaker": {"Lambda": 50000, "d": 2},
  "None": {} # For the case where the user untics the box
}

# Dictionary of the normalisation methods, parameters, and default values
# Ip normalisation not here!
normalisation_options = {
    "Scale": None,
    "Maximum in range": {"Range start": 890, "Range end": 910},
    "Maximum whole array": None,
    "Whole array": None,
    "Single point": {"Wavenumber": 890},
    "Area": {"Area start": 890, "Area end": 910},
    "Interpolate area": {"Area start": 890, "Area end": 910},
    "Max in interpolation range": {"Normalisation indexes": [890, 910]},
    "Ip Normalisation": {"Normalisation indexes": 10},
    "None": {}
}

# Dictionary of the baseline methods, parameters, and default values
baseline_options = {
    "ALS": {"Lambda": 10**4, "p": 0.01, "Maximum iterations": 10},
    "airPLS": {"Lambda": 30, "p": 1, "Maximum iterations": 15},
    "FFT DFT": {"Fourier Values": 3},
    "FFT RDFT": {"Fourier Values": 3},
    "FFT Hermitian": {"Fourier Values": 3},
    "ModPoly": {"Polynomial": 3},
    "IModPoly": {"Polynomial": 3},
    "Zhang": {"Polynomial": 3},
    "SWiMA": {"Maximum iterations": 10},
    "RollingBall": {"Ball Height": 0.1, "Ball Width": 25},
    "Average": {"Window size": 100},
    "Manual": {"Mask": False, "Polynomial order": 20},
    "None": {}
}

# Dictionary of the CRS methods, parameters, and default values
crs_options = {
    "Modified Z-Score": {"Threshold score": 13, "Differential order": 4, "Width": 3, "Graphing": 0},
    "Wavenumber": {"Threshold wavenumber": 5},
    "Differential": {"Threshold differential": 5},
    "Consensus": {"Threshold wavenumber": 5, "Threshold differential": 5},
    "None": {}
}

# Dictionary containing previous dictionaries in operation order
process_menu = {
    "Remove Cosmic Ray Spikes": crs_options,
    "Smooth": smoothing_options,
    "Normalise": normalisation_options,
    "Baseline Correction": baseline_options    
}

# Create a dictionary of the frames that can be referenced when making labels etc.
dict_frames = {
  "Remove Cosmic Ray Spikes": crs_frame,
  "Smooth": smooth_frame,
  "Normalise": normalise_frame,
  "Baseline Correction": baseline_frame
}

# The preset processing packages including methods and parameters so that they can be displayed in the "variables" frame
# All need to have the same format even when it doesn't make immediate sense
# e.g. need to have "alignment indexes: alignment indexes: [num, num]" so the parameter is the required number of dictionaries deep
processing_package = {
    "Default": {
      "Remove Cosmic Ray Spikes": {"Modified Z-Score":{"Threshold score": 13, "Differential order": 4, "Width": 3, "Graphing": 0}},
      "Smooth": {"Whittaker":{"Lambda":28, "d":3}},
      "Normalise": {"Ip Normalisation":{"Normalisation indexes":1}},
      "Baseline Correction": {"airPLS": {"Lambda": 10, "porder": 1, "Maximum iterations": 50}},
    },
    "DNA Normalisation": {
      "Remove Cosmic Ray Spikes": {"wavenumber":{"threshold wavenumber":5}},
      "Smooth": {"FFT":{"Fourier Values":300}},
      "Normalise": {"area":{"Normalisation indexes":[770,790]}},
      "Baseline Correction": {"ALS":{"lam": 10**4.5, "p": 0.01, "niter": 10}}
    },
    "Heavy Normalisation": {
      "Remove Cosmic Ray Spikes": {"wavenumber":{"threshold wavenumber":5}},
      "Smooth": {"FFT":{"Fourier Values":300}},
      "Normalise": {"area":{"Normalisation indexes":[990,1015]}},
      "Baseline Correction": {"ALS":{"lam": 10**4.5, "p": 0.01, "niter": 10}}      
    }
}

decompositions = {
   "Multivariate Curve Resolution": multivariate_resolution,
   "Non-Negative Elastic Net": nonneg_elastic_net
}

regressor_options = ['OLS', 'NNLS', 'Linear Regressor', 'None']

# Store the values auto-selected when the decomp method dropdown is changed
decomp_defaults = {"Multivariate Curve Resolution":
                      {"Maximum iterations": 50,
                       "Spectral regressor":regressor_options[3],
                       "Concentration regressor":regressor_options[1],
                       "Sum to 1 constraint": False,
                       "Non-negative constraint": False}
}

# Stores all possible values, mainly to store options to pass to generated OptionMenus
# Reference defaults where possible rather than having random values
decomp_option_dict = {
   "Multivariate Curve Resolution":
        {"Maximum iterations": decomp_defaults["Multivariate Curve Resolution"]["Maximum iterations"],
         "Spectral regressor":regressor_options,
         "Concentration regressor":regressor_options,
         "Sum to 1 constraint":decomp_defaults["Multivariate Curve Resolution"]["Sum to 1 constraint"]},
   "Non-Negative Elastic Net":
        {"Regressors": False}
}

#-----------------------------------------------------------------------------------
# Making the widgets for the GUI and placing in the grid
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Frame Files
#-----------------------------------------------------------------------------

# Could just display filename in box but then we can't take the paths from the listbox
# Showing the whole path can make it easier to identify which files have been selected
file_x_scrollbar = ttk.Scrollbar(frame_files, orient='horizontal')
file_y_scrollbar = ttk.Scrollbar(frame_files, orient='vertical')
file_list_display = tk.Listbox(frame_files, width=60, xscrollcommand=file_x_scrollbar.set, yscrollcommand=file_y_scrollbar.set)

file_list_display.grid(row=1, column=0, rowspan=2, columnspan=2, sticky=tk.NSEW, padx=16)
file_x_scrollbar.grid(row=3, column=0, columnspan=2, sticky=tk.EW, padx=10, pady=5)
file_x_scrollbar.config(command = file_list_display.xview)
file_y_scrollbar.grid(row=1, column=1, rowspan=2, sticky="nse")
file_y_scrollbar.config(command= file_list_display.yview)

processing_label = tk.Label(frame_files, text="Processing Method")
processing_label.grid(row=6, column=0, sticky=tk.E, padx=5, pady=5)

# The dropdown for selecting the package, takes first key as default, keys as the options
default_package = "Default" # Make this more flexible, shouldn't rely on a package called "Default" existing
package_option = tk.StringVar()
package_option.set(default_package)
package_option_menu = tk.OptionMenu(frame_files, package_option, *processing_package.keys(), command=lambda chosen_package: package_change(chosen_package))
package_option_menu.grid(row=6, column=1)

# Label for "select files" that appears at the top of the page
file_label = tk.Label(frame_files, text="Select mixture files")
file_label.grid(column=0, row=0, sticky=tk.W, padx=5, pady=5)

# The "Open files" button
open_button = tk.Button(
    master=frame_files,
    text='Open Files',
    command=select_files
)
open_button.grid(row=1, column=2, sticky=tk.W)

# The "Remove file" button
remove_button = tk.Button(
    master=frame_files,
    text='Remove file',
    command=remove_files
)
remove_button.grid(row=2, column=2, sticky=tk.NW)

# The "Remove all files" button
removeAll_button = tk.Button(
    master=frame_files,
    text='Remove all files',
    command=removeAll_files
)
removeAll_button.grid(row=2, column=2, sticky=tk.SW)

# Button to show/hide the processing frame
hide_button = tk.Button(
    master=frame_files,
    text='Show/hide variables',
    command=show_hide
)
hide_button.grid(row=6, column=2, padx=5, pady=5)

frame_legend = ttk.Frame(master=frame_files)
units_var = tk.StringVar()
unitsLabel = tk.Label(master=frame_legend, text="Units (e.g. %)")
unitsEntry = tk.Entry(master=frame_legend, textvariable=units_var, width=2)
charRequestLabel = tk.Label(master=frame_legend, text="Collect values from before")
legendChar = tk.StringVar()
charRequestEntry = tk.Entry(master=frame_legend, textvariable=legendChar, width=2)
collectedLabel = tk.Label(master=frame_legend, text="Collected values")
collectedEntry = tk.Entry(master=frame_legend)
legendUpdateButton = tk.Button(
   master=frame_legend,
   text="Update", 
   command= lambda: updateLegendVals())
# Pack in the above widgets
unitsLabel.grid(row=0, column=0)
unitsEntry.grid(row=0, column=1, sticky='w')
charRequestLabel.grid(row=1, column=0)
charRequestEntry.grid(row=1, column=1, sticky='w')
collectedLabel.grid(row=2, column=0)
collectedEntry.grid(row=2, column=1, sticky='w')
legendUpdateButton.grid(row=2, column=3)

def character_limit(entry_text):
   if len(entry_text.get()) > 0:
      entry_text.set(entry_text.get()[-1])
legendChar.trace("w", lambda *args: character_limit(legendChar))

labelledDataVar = tk.IntVar()
labelled_checkbox = tk.Checkbutton(
   master = frame_files,
   text="Labelled data",
   variable=labelledDataVar,
   command= lambda: labelled_true()
)
labelled_checkbox.grid(row=4, column=0)


#-----------------------------------------------------------------------------
# Frame Processing
#-----------------------------------------------------------------------------

# Create header label for the processing frame
variables_label = tk.Label(frame_processing, text="Variables")
variables_label.grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)

# Instancing the menus and boxes as dictionaries and dropdown as a list
all_option_menus = {}
all_option_boxes = {}
dropdown_list = {}  # this is where I want to get my stringvars from
checkbox_list = {}

# For each process and method create a dropdown and checkbox
for process, methods in process_menu.items():
  dropdown_option = tk.StringVar()
  dropdown_option.set(process)
  dropdown_list[process] = dropdown_option
  checkbox_option = tk.BooleanVar(value=True)
  checkbox_list[process] = checkbox_option
  # Add to dictionary under the process name, lambda function is tied to each menu individually
  all_option_menus[process] = tk.OptionMenu(frame_processing, dropdown_option, *methods.keys(), command=lambda chosen_method, process=process: change_manually(chosen_method,process))
  all_option_boxes[process] = tk.Checkbutton(frame_processing, text=process, variable=checkbox_list[process])

# Add in all checkboxes, option menus, and the parameters associated with the auto-selected methods
for count, process in enumerate(process_menu.keys()):
  all_option_boxes[process].grid(row=1+(3*count), column=0, sticky=tk.W)
  all_option_menus[process].grid(row=2+(3*count), column=1)

  # installing the parameter frames into the grid
  dict_frames[process].grid(row=3+(3*count), column=1)
  
# Calling the change through package to start the program in the default state (1st package)
# Should have the dropdowns, labels, and entry boxes filled in with the appropriate values
package_change(default_package)
change_thu_decomp("Multivariate Curve Resolution", decomp_option_dict)

#-----------------------------------------------------------------------------
# Frame Display
#-----------------------------------------------------------------------------
# Option to import, otherwise this default is used
graph_parameters = {'title': False,
                    'x-axis': 'Wavenumbers (cm$^{-1}$)',
                    'y-axis': 'Intensity (AU)',
                    'show y=0': 0,
                    'Show legend': 0,
                    'figSizeX': False,
                    'figSizeY': False,
                    'titleFontSize': False,
                    'axesTitleSize': False,
                    'axesValuesSize': False,
                    'legendFontSize': False}

graph_button = tk.Button(
    master=tab1,
    text='Plot graph',
    command= lambda: plot_graph()
)

ppsave_button = tk.Button(
   master=tab1,
   text="Save preprocessed spectra",
   command= lambda: save_preprocessed()
)


# Splitting the graphing
frame_ppGraphLabels = ttk.Frame(master=frame_display, relief=tk.SOLID)
frame_ppFormat = ttk.Frame(master=frame_display, relief=tk.SOLID)
frame_ppFormat.grid_columnconfigure(0, weight=4)
frame_ppFormat.grid_columnconfigure(1, weight=1)
frame_ppFormat.grid_columnconfigure(2, weight=1)
frame_ppFormat.grid_columnconfigure(3, weight=1)

# Making all the labels and boxes for the graphing options
# This is long though isn't necessarily bad...
#--Graph labels---------------------------------------------------------------
# Overall heading label
ppGraphLabelsLabel = tk.Label(master=frame_ppGraphLabels, text="Graph labels", font=("Arial", 10, "bold"))
ppGraphLabelsLabel.grid(row=0, column=0, columnspan=3, pady=2, padx=60)     
# Label, entry box, and entry box insertion for figure title
titleLabel = tk.Label(master=frame_ppGraphLabels, text="Figure Title")
titleEntry = tk.Entry(master=frame_ppGraphLabels)
if graph_parameters["title"]:  # Value of False means the box is empty
    titleEntry.insert(tk.END, graph_parameters["title"])
# Label, entry box, and entry box insertion for x-axis
xLabel = tk.Label(master=frame_ppGraphLabels, text="x-axis")
xEntry = tk.Entry(master=frame_ppGraphLabels)
if graph_parameters["x-axis"]:
    xEntry.insert(tk.END, graph_parameters["x-axis"])
# Label, entry box, and entry box insertion for y-axis
yLabel = tk.Label(master=frame_ppGraphLabels, text="y-axis")
yEntry = tk.Entry(master=frame_ppGraphLabels)
if graph_parameters["y-axis"]:
    yEntry.insert(tk.END, graph_parameters["y-axis"])
# IntVar creation and linking with checkbox for x axis
showaxis = tk.IntVar(value=graph_parameters["show y=0"])
checkx = tk.Checkbutton(master=frame_ppGraphLabels, text="Show y=0", variable=showaxis)
# IntVar creation and linking with checkbox for legend
showLeg = tk.IntVar(value=graph_parameters["Show legend"])
checkLeg = tk.Checkbutton(master=frame_ppGraphLabels, text="Show legend", variable=showLeg)
# Fitting the labels and entries into the grid
# Could be done as for processing frame to allow for more labels and boxes
titleLabel.grid(row=1, column=0)
titleEntry.grid(row=1, column=1)
xLabel.grid(row=2, column=0)
xEntry.grid(row=2, column=1)
yLabel.grid(row=3, column=0)
yEntry.grid(row=3, column=1)
checkx.grid(row=4, column=1, sticky=tk.W)
checkLeg.grid(row=5, column=1, sticky=tk.W, pady=2)

#--Formatting options---------------------------------------------------------
ppFormatLabel = tk.Label(master=frame_ppFormat, text="Formatting options", font=("Arial", 10, "bold"))
ppFormatLabel.grid(row=0, column=0, columnspan=4, pady=2)
# Label, entry box, and entry box insertion for title font size
titleFontLabel = tk.Label(master=frame_ppFormat, text="Title font size")
titleFontEntry = tk.Entry(master=frame_ppFormat, width=8)
if graph_parameters["titleFontSize"]:  # Value of False means the box is empty
    titleFontEntry.insert(tk.END, graph_parameters["titleFontSize"])
# Label, entry box, and entry box insertion for figure dimensions
figureDimsLabel = tk.Label(master=frame_ppFormat, text="Figure dimensions")
figureDimsSeparator = tk.Label(master=frame_ppFormat, text="x")
figureDimsEntryX = tk.Entry(master=frame_ppFormat, width=3)
figureDimsEntryY = tk.Entry(master=frame_ppFormat, width=3)
if graph_parameters["figSizeX"]:  # Value of False means the box is empty
    figureDimsEntryX.insert(tk.END, graph_parameters["figSizeX"])
if graph_parameters["figSizeY"]:  # Value of False means the box is empty
    figureDimsEntryY.insert(tk.END, graph_parameters["figSizeY"])
# Label, entry box, and entry box insertion for axes title size
axesTitleLabel = tk.Label(master=frame_ppFormat, text="Axes title size")
axesTitleEntry = tk.Entry(master=frame_ppFormat, width=8)
if graph_parameters["axesTitleSize"]:  # Value of False means the box is empty
    axesTitleEntry.insert(tk.END, graph_parameters["axesTitleSize"])
# Label, entry box, and entry box insertion for axes values size
axesValuesLabel = tk.Label(master=frame_ppFormat, text="Axes values size")
axesValuesEntry = tk.Entry(master=frame_ppFormat, width=8)
if graph_parameters["axesValuesSize"]:  # Value of False means the box is empty
    axesValuesEntry.insert(tk.END, graph_parameters["axesValuesSize"])
# Label, entry box, and entry box insertion for legend font size
legendFontLabel = tk.Label(master=frame_ppFormat, text="Legend font size")
legendFontEntry = tk.Entry(master=frame_ppFormat, width=8)
if graph_parameters["legendFontSize"]:  # Value of False means the box is empty
    legendFontEntry.insert(tk.END, graph_parameters["legendFontSize"])

titleFontLabel.grid(row=1, column=0, sticky=tk.W, padx=2)
titleFontEntry.grid(row=1, column=1, columnspan=3, sticky=tk.E, padx=2)
figureDimsLabel.grid(row=2, column=0, sticky=tk.W, padx=2)
figureDimsEntryX.grid(row=2, column=1)
figureDimsSeparator.grid(row=2, column=2)
figureDimsEntryY.grid(row=2, column=3, padx=2)
axesTitleLabel.grid(row=3, column=0, sticky=tk.W, padx=2)
axesTitleEntry.grid(row=3, column=1, columnspan=3, sticky=tk.E, padx=2)
axesValuesLabel.grid(row=4, column=0, sticky=tk.W, padx=2)
axesValuesEntry.grid(row=4, column=1, columnspan=3, sticky=tk.E, padx=2)
legendFontLabel.grid(row=5, column=0, sticky=tk.W, padx=2)
legendFontEntry.grid(row=5, column=1, columnspan=3, sticky=tk.E, pady=2, padx=2)

frame_ppGraphLabels.grid(row=0, column=0, padx=8, sticky=tk.N)
frame_ppFormat.grid(row=0, column=1, padx=8, sticky=tk.N)
graph_button.grid(row=4, sticky="nsew")
ppsave_button.grid(row=5, pady=5)

#-----------------------------------------------------------------------------
# Multivariate Curve Resolution Frame
#-----------------------------------------------------------------------------

component_a_label = tk.Label(filesdecomp_frame, text="Substrate + Water Spectrum") 
component_a_label.grid(row=0, column=0, sticky=tk.W, padx=5)  
a_file_entry = tk.Entry(master=filesdecomp_frame)
a_file_entry.grid(row=1, column=0, columnspan=2, sticky=tk.EW, padx=50)
a_file_button = tk.Button(master=filesdecomp_frame, text="Select file", command=lambda:select_component(a_file_entry))
a_file_button.grid(row=1, column=2)

component_b_label = tk.Label(filesdecomp_frame, text="Substrate + Water + Analyte Spectrum") 
component_b_label.grid(row=2, column=0, sticky=tk.W, padx=5)  
b_file_entry = tk.Entry(master=filesdecomp_frame)
b_file_entry.grid(row=3, column=0, columnspan=2, sticky=tk.EW, padx=50)
b_file_button = tk.Button(master=filesdecomp_frame, text="Select file", command=lambda:select_component(b_file_entry))
b_file_button.grid(row=3, column=2)

units_label = tk.Label(filesdecomp_frame, text="Units (e.g. %)")
units_label.grid(row=4, column=0, padx=5, sticky=tk.W)
units_entry = tk.Entry(filesdecomp_frame, textvariable=units_var, width=4)
units_entry.grid(row=4, column=0, sticky=tk.E)

reference_quantity_label = tk.Label(filesdecomp_frame, text="Reference quantity:")
reference_quantity_label.grid(row=5, column=0, padx=5, sticky=tk.W)
reference_quantity_entry = tk.Entry(master=filesdecomp_frame, width=4)
reference_quantity_entry.grid(row=5, column=0, sticky=tk.E)
reference_units_label = tk.Label(filesdecomp_frame, textvariable=units_var)
reference_units_label.grid(row=5, column=1, sticky=tk.W)

# Create header label for the frame of decomposition parameters
decomp_method_label = tk.Label(filesdecomp_frame, text="Decomposition method")
decomp_method_label.grid(row=10, column=0, sticky=tk.W, padx=5, pady=20)

decomp_dropdown_list = list(decompositions.keys()) # MCR, NNEN etc.
decomp_dropdown_default = tk.StringVar()
decomp_dropdown_default.set(decomp_dropdown_list[0])  # i.e. MCR
# Can we just not have a preselected method and have it start blank?
decomp_option_menu = tk.OptionMenu(filesdecomp_frame, decomp_dropdown_default, *decomp_dropdown_list, command=lambda chosen_decomp : change_thu_decomp(chosen_decomp, decomp_option_dict))
decomp_option_menu.grid(row=10, column=1, pady=20)

frame_decomp_graph = ttk.Frame(tab2, relief="solid")

decomp_labelled_check = tk.Checkbutton(tab2, text="Labelled data", variable=labelledDataVar, command= lambda: labelled_true())
decomp_labelled_check.grid(row=2, column=0, pady=5)
decomp_plot_var = tk.IntVar()
decomp_plot_var.set(True)
decomp_plot_check = tk.Checkbutton(frame_decomp_graph, variable=decomp_plot_var, text="Plot Graph")
decomp_plot_check.grid(row=0, column=0, padx=3, pady=3, sticky=tk.W)
log_dropdown_label = tk.Label(frame_decomp_graph, text="Log / linear graph")
log_dropdown_label.grid(row=1, column=0, padx=3)
scale_options = ["Linear", "Log"]
scale_var = tk.StringVar()
scale_var.set(scale_options[0])
log_dropdown = tk.OptionMenu(frame_decomp_graph, scale_var, *scale_options)
log_dropdown.grid(row=1, column=1, padx=3)
xy_var = tk.IntVar()
xy_line = tk.Checkbutton(frame_decomp_graph, variable=xy_var, text="Show y=x line")
xy_line.grid(row=2, column=0, padx=3, sticky=tk.W)
least_sq_var = tk.IntVar()
least_sq_var.set(True)
least_sq_line = tk.Checkbutton(frame_decomp_graph, variable=least_sq_var, text="Show least squares line")
least_sq_line.grid(row=3, column=0, padx=3, sticky=tk.W)
stats_var = tk.IntVar()
stats_var.set(True)
output_stats = tk.Checkbutton(frame_decomp_graph, variable=stats_var, text="Output statistics")
output_stats.grid(row=4, column=0, padx=3, pady=3, sticky=tk.W)

find_elements_button = tk.Button(master=tab2, text="Calculate", command=calculate_elements)
find_elements_button.grid(row=4, column=0, pady=25)

save_parameters_button = tk.Button(master=tab2, text="Save processing setup", command=save_parameters)
save_parameters_button.grid(row=5, column=0, sticky=tk.W, padx=50, pady=20)

save_concentration_button = tk.Button(master=tab2, text="Save results", command=save_concentrations)
save_concentration_button.grid(row=5, column=0, sticky=tk.E, padx=50, pady=20)

#-----------------------------------------------------------------------------
# Window display
#-----------------------------------------------------------------------------

# Assembling the second level frames onto the root grid
frame_files.grid(row=0, column=0)
frame_processing.grid(row=2, column=0)
frame_display.grid(row=3, column=0, pady=25)

filesdecomp_frame.grid(row=0, column=0)
decomp_frame.grid(row=1, column=0)
decomp_method_frame.grid(row=0, column=0)

# Display the window
root.mainloop()