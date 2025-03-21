#!/usr/bin/env python
# coding: utf-8

# # Extract the data from file 

# In[2]:


# Extract data 

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import re

# Load the CSV file
file_path = "/Users/ailsamoody/Desktop/HPSEC/060125/HPLC_data_visual.csv"
df = pd.read_csv(file_path) 


# In[6]:


# Filter out standard files (which don't follow fraction_number-sample_name format)
df_filtered = df[df['filename'].str.contains(r"^\d+-[A-Za-z0-9]+", regex=True)].copy()

print(df_filtered)


# In[8]:


df_filtered['filename'].unique()


# # Clear name so there is no spaces 

# In[11]:


df_filtered['processed names'] = (
    df_filtered['filename']
    .str.replace('NF-R', 'NFR')
    .str.replace('BOOSTER-FD', 'BOOSTERFD')
    .str.replace('PPC-QA', 'PPCQA')
)
df_filtered['processed names'].unique()


# # Extract the information from the fixed filenames

# In[16]:


def extract_sample_fraction(s):
    # Find first and second hyphens
    first_hyphen = s.find('-')
    second_hyphen = s.find('-', first_hyphen + 1)
    
    # Extract the sample name
    sample_name = s[first_hyphen + 1 : second_hyphen]
    
    # Extract everything from the second hyphen up to '_S_'
    s_s_index = s.find('_S_')
    fraction_part = s[second_hyphen + 1 : s_s_index] if s_s_index != -1 else ""
    
    return sample_name, fraction_part

df_filtered['sample_name'], df_filtered['fraction_part'] = zip(*df_filtered['processed names'].map(extract_sample_fraction))


# In[18]:


df_filtered.head()


# # Extract the fraction numbers as fraction low and fraction high for samples with fraction range

# In[21]:


def parse_fraction_part(x):
    if '-' in x:
        low, high = x.split('-')
        return low, high
    return x, x

df_filtered['fraction_low'], df_filtered['fraction_high'] = zip(*df_filtered['fraction_part'].map(parse_fraction_part))
df_filtered.head()


# In[23]:


df_plot = df_filtered[['Area %', 'kDa', 'sample_name', 'fraction_low', 'fraction_high']].copy()
df_plot['fraction_low'] = pd.to_numeric(df_plot['fraction_low'], errors='coerce')
df_plot.head()


# # Save cleaned data as csv

# In[28]:


# Save the cleaned data
df_plot.to_csv("cleaned_HPLC_data.csv", index=False)


# # Create interactive plot 3D

# In[60]:


import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import plotly.subplots as psub
import re


# In[162]:


# Generate a 2x2 grid of 3D plots
sample_names = df_plot['sample_name'].unique()
fig = psub.make_subplots(rows=2, cols=2, specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}], [{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
                         subplot_titles=sample_names, horizontal_spacing=0.02, vertical_spacing=0.02)


# Manually adjust annotation y-positions to move titles lower
for annotation in fig['layout']['annotations']:
    annotation['y'] -= 0.03

positions = [(1, 1), (1, 2), (2, 1), (2, 2)]
scenes = {}
traces = []

# Set axis limits manually
x_min, x_max = 0, 33
y_min, y_max = 0, 1500
z_min = df_filtered['Area %'].min()
z_max = df_filtered['Area %'].max()

for i, (sample, (row, col)) in enumerate(zip(sample_names, positions)):
    df_sample = df_plot[df_plot['sample_name'] == sample]
    scatter3d = go.Scatter3d(
        x=df_sample['fraction_low'],
        y=df_sample['kDa'],
        z=df_sample['Area %'],
        mode='markers',
        marker=dict(
            size=6,
            color=df_sample['Area %'],  # Color by Area %
            colorscale=[(0, 'blue'), (1, 'red')],
            opacity=0.6,
            showscale=False  # Disable multiple colorbars
        ),
        name=sample,
        text=df_sample['sample_name'],
        showlegend=False  # Remove legend entries
    )
    fig.add_trace(scatter3d, row=row, col=col)
    scenes[f'scene{i+1}'] = dict(
        xaxis=dict(title='Fraction Number', titlefont=dict(size=10), range=[x_min, x_max]),
        yaxis=dict(title='Molecular Weight (kDa)', titlefont=dict(size=10), range=[y_min, y_max]),
        zaxis=dict(title='Area (%)', titlefont=dict(size=10))
    )

# Add a single color bar with correct scale
fig.add_trace(go.Scatter3d(
    x=[None], y=[None], z=[None],
    mode='markers',
    marker=dict(
        size=8,
        color=[df_plot['Area %'].min(), df_plot['Area %'].max()],
        colorscale=[(0, 'blue'), (1, 'red')],
        cmin=df_plot['Area %'].min(),
        cmax=df_plot['Area %'].max(),
        showscale=True,
        colorbar=dict(title='Area (%)', x=1.1)
    ),
    showlegend=False  # Ensure this does not appear in the legend
))

# Layout configuration
fig.update_layout(
    title='',
    margin=dict(l=10, r=10, t=10, b=10),
    showlegend=False,
    **scenes
)




# In[164]:


# Save as interactive HTML file
html_file = "interactive_chromatography.html"
pio.write_html(fig, html_file)

print(f"Filtered data saved as cleaned_HPLC_data.csv")
print(f"Interactive visualization saved as {html_file}. Open it in a browser to explore!")


# In[ ]:




