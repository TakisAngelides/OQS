import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib import rc
rc('text', usetex=True)

plt.tight_layout()

# Define the models
def exponential_saturation(aT, T_inf, k):
    return T_inf * (1 - np.exp(-k * aT))

def hyperbolic_saturation(aT, T_inf, K):
    return (T_inf * aT) / (K + aT)

def logarithmic_growth(aT, T_0, k):
    return T_0 + k * np.log(1 + aT)

def power_law(aT, a, b):
    return a * aT**b

# Directory containing the pickle files
data_dir = "Thermalization_Heatmap_50"

# Data structure to store the thermalization times
data = {}

# Step 1: Reading data from pickle files
for filename in os.listdir(data_dir):
    if filename.endswith(".pickle"):
        # Extract parameter values from the filename
        try:
            parts = filename.replace(".pickle", "").split("_")
            ma = float(parts[0])
            l_0_1 = float(parts[1])
            aD = float(parts[2])
            aT = float(parts[3])
        except ValueError:
            print(f"Filename {filename} does not match the expected format.")
            continue

        # Read the thermalization time from the pickle file
        with open(os.path.join(data_dir, filename), "rb") as f:
            thermalization_time = pickle.load(f)

        # Organize data
        key = (ma, l_0_1, aD)
        if key not in data:
            data[key] = []
        data[key].append((aT, thermalization_time))

# Define a list of markers for differentiation
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'H', '+', 'x']
marker_cycle = iter(markers)
color_cycle = iter(colors)

# Step 2: Fit data for each group and plot
# Adjustments for better aesthetics
plt.figure(figsize=(12, 8))

# Set global font size for ticks and axis labels
plt.rcParams.update({
    'axes.labelsize': 20,      # Axes label font size
    'xtick.labelsize': 20,     # X-axis tick font size
    'ytick.labelsize': 20,     # Y-axis tick font size
})

# Set tick direction
plt.tick_params(axis='both', direction='in')
models = {
    "Exponential Saturation": exponential_saturation,
    "Hyperbolic Saturation": hyperbolic_saturation,
    "Logarithmic Growth": logarithmic_growth,
    "Power Law": power_law
}

# Sort the keys of the data dictionary
sorted_keys = sorted(data.keys(), key=lambda x: (x[0], x[1], x[2]))

# Loop over each group of data in sorted order
for (ma, l_0_1, aD) in sorted_keys:
    values = data[(ma, l_0_1, aD)]
    
    # Sort values by aT
    values.sort(key=lambda x: x[0])
    aT_values, thermalization_times = zip(*values)
    aT_values = np.array(aT_values)
    thermalization_times = np.array(thermalization_times)
    
    # Fit data to each model
    # fit_results = {}
    # for label, model in models.items():
    #     try:
    #         popt, _ = curve_fit(model, aT_values, thermalization_times, maxfev=10000)
    #         fitted_values = model(aT_values, *popt)
    #         r2 = r2_score(thermalization_times, fitted_values)

    #         # Store results
    #         fit_results[label] = {
    #             "parameters": popt,
    #             "r2": r2,
    #             "fitted_values": fitted_values,
    #         }
    #     except RuntimeError:
    #         print(f"{label} model fitting failed for ma={ma}, l_0_1={l_0_1}, aD={aD}.")

    # Determine the best fit
    # if fit_results:
    #     best_fit_label = max(fit_results, key=lambda k: fit_results[k]["r2"])
    #     best_fit = fit_results[best_fit_label]

    # Assign a unique marker for this group
    marker = next(marker_cycle)
    
    color = next(color_cycle)
    
    # Plot the original data and the best fit
    plt.errorbar(aT_values, thermalization_times, yerr = 0.1, label=f"$m$ = {ma}, $l_0$ = {l_0_1}, $D$ = {aD}", fmt = marker, zorder=5, color = color, markerfacecolor='none', markersize=10, linestyle = '-')
        
    # plt.plot(aT_values, best_fit["fitted_values"], linestyle="--", linewidth=2, color = color)
    # print(f"Best Fit for ma={ma}, l_0_1={l_0_1}, aD={aD}: {best_fit_label}")
    # print(f"Parameters: {best_fit['parameters']}")
    # print(f"R^2: {best_fit['r2']:.4f}")
    
        
# Step 3: Finalize the plot
plt.xlabel("$T$")
plt.ylabel("$\mathcal{T}$")
plt.legend(loc="upper center", fontsize=12, ncol=4, bbox_to_anchor=(0.5, 1.125))
plt.savefig("thermalization_time_vs_aT.pdf", bbox_inches='tight', dpi = 1200)
