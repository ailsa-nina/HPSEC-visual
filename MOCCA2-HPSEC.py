#!/usr/bin/env python
# coding: utf-8

# In[26]:


import glob
import os
import pandas as pd
import mocca2 as mc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import re


# ## Part 1: Process file names to extract samples metadata

# In[28]:


def process_filename(filepath):
    """ Given the name of an HPLC file, extract the sample name, ID and wavelength.
    For example, given "37-BOOSTER-FD-32_S_214.CSV", return {'name': 'BOOSTER-FD', 'id': 32', 
    'wavelength': 214}. 
    
    If the filename does not match the expected format, return None.

    """
    name = os.path.basename(filepath)
    # Regular expression to extract the main components
    pattern = r"^(?:\d+-)?([A-Z0-9\-]+)_S_(\d+)\.CSV$"

    match = re.match(pattern, name)
    if match:
        full_sample_name = match.group(1)
        wavelength = match.group(2)

        # Split sample name into name and ID if it ends with numbers
        name_match = re.match(r"^(.*?)(-\d{2}-\d{2}|\-\d+)?$", full_sample_name)
        if name_match:
            sample_name = name_match.group(1)  # Main name
            sample_id = name_match.group(2)   # Optional ID (e.g., "-21-22" or "-30")
            if sample_id:
                sample_id = sample_id[1:]  # Remove leading hyphen
        return({'name': sample_name, 'id': sample_id, 'wavelength': wavelength, 'path': filepath})
    else:
        return None


# In[29]:


HPLC_filepaths = list(glob.glob("./data/*.CSV"))
processed_names = list(map(process_filename, HPLC_filepaths))
df_names = pd.DataFrame([x for x in processed_names if x is not None])


# In[30]:


# Sanity check
df_names.head()


# ## Separate the baselines in their own dataframe

# In[32]:


baseline_ids = df_names["name"].str.contains("PBS")
df_baseline_names = df_names[baseline_ids]  # Rows where "name" contains "PBS"

# Remove these rows from the original DataFrame
df_hplc_names = df_names[~baseline_ids]


# ## Part 2: For each sample, process peaks

# In[34]:


def integrate_peaks(chromatogram, MIN_PEAK_HEIGHT=2):
    """ Given a chromatogram, find the peaks, plot them and return the integrated peak values.

    """
    chromatogram.correct_baseline()
    chromatogram.find_peaks(min_height=MIN_PEAK_HEIGHT)
    chromatogram.deconvolve_peaks(
        model="FraserSuzuki",
        min_r2=0.95,
        relaxe_concs=False,
        max_comps=7)
    # Plot for verification.
    chromatogram.plot()
    plt.show()
    data = pd.DataFrame([{
        "Elution Time": chromatogram.time[component.elution_time],
        "Integral": component.integral,}
        for component in chromatogram.all_components()])
    try:
        data["Area %"] = np.round(data["Integral"] / data["Integral"].sum() * 100, decimals=2)
    except: # if no peaks
        data["Area %"] = 0.0
    return data


# In[35]:


def process_file_peaks(filepath, freq, blank_path, MIN_PEAK_HEIGHT=2):
    df = pd.read_csv(filepath, sep='\t', encoding='utf-16le', names=['t', 'A'])
    mc_data = mc.classes.Data2D(df["t"].to_numpy(), np.array([freq]), df["A"].to_numpy().reshape(1, -1))

    try:
        df_baseline = pd.read_csv(blank_path, sep='\t', encoding='utf-16le', names=['t', 'A'])
        mc_data_baseline = mc.classes.Data2D(
            df_baseline["t"].to_numpy(), np.array([freq]),
            df_baseline["A"].to_numpy().reshape(1, -1))
        chromatogram = mc.Chromatogram(sample = mc_data, blank = mc_data_baseline,
                                   interpolate_blank=True)
    except:
        print("Error") # ugly but Ok for today.
    peak_data = integrate_peaks(chromatogram, MIN_PEAK_HEIGHT)
    peak_data["filename"] = os.path.basename(filepath) # add path information for latter use
    return peak_data


# In[36]:


def _apply(row):
    w = int(row['wavelength'])
    if w == 214:
        baseline_path = './data/2PBS_S_214.CSV'
    elif w == 280:
        baseline_path = './data/2PBS_S_280.CSV'

    print(f"{row['name']} / {row['id']} wavelength: {w}")
    res = process_file_peaks(row['path'], w, baseline_path)
    return res
res = df_hplc_names.apply(_apply, axis = 1)


# # Part 3: Postprocess results
# We add peak numbering and concatenate everything in a single dataframe.

# In[38]:


df_results = pd.concat([x for x in res])
df_results['peak nr'] = df_results.index + 1
df_results.reset_index(inplace=True, drop=True)
df_results.head()


# In[39]:


# Save
df_results.to_csv("HPLC_peaks_processed_Ailsa_Jan2025.csv")


# In[ ]:




