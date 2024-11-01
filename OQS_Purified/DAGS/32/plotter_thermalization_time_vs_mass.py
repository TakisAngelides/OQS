import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.pyplot import xticks, yticks
import re

# Define the directory and file names
file_names = [
    'Plots/Thermalization_Heatmap/data/0.1.pickle', 
    'Plots/Thermalization_Heatmap/data/0.25.pickle', 
    'Plots/Thermalization_Heatmap/data/0.5.pickle', 
    'Plots/Thermalization_Heatmap/data/0.75.pickle', 
    'Plots/Thermalization_Heatmap/data/1.0.pickle'
]

# Load the data from the pickle files
plots = []
masses = [0.1, 0.25, 0.5, 0.75, 1.0]

for file_name in file_names:
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        plots.append(data)
        
# Prepare to collect data for plotting
thermalization_data = {}

# Iterate through each DataFrame and each l_0_1 value
for plot, mass in zip(plots, masses):
    if isinstance(plot, pd.DataFrame):
        for i, index in enumerate(plot.index):
            for j, column in enumerate(plot.columns):
                
                thermalization_data[(mass, i, j)] = plot.iloc[i, j]
                
for i, index in enumerate(plot.index):
    for j, column in enumerate(plot.columns):
        l = []
        for mass in masses:
            l.append(thermalization_data[(mass, i, j)])
        plt.plot(l)
        plt.title(f'D = {index}, l_0 = {column}')
        plt.savefig(f'Plots/thermalization_vs_mass/d_{index}_l0_{column}.png')
        plt.close()
            
