import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
import pickle
from collections import defaultdict
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('text', usetex=True)
plt.tight_layout()

# Set the directory where files are located
directory_path = '/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_50'

# Load data for the first plot
with open(directory_path + '/PN_1.0_0.5_5_7.0.pickle', "rb") as f:
    pn7 = np.real(pickle.load(f))
with open(directory_path + '/PN_1.0_0.5_5_100.0.pickle', "rb") as f:
    pn100 = np.real(pickle.load(f))

# Load data for the second plot
with open(directory_path + '/KE_1.0_0.5_5_7.0.pickle', "rb") as f:
    KE7 = np.real(pickle.load(f))
with open(directory_path + '/KE_1.0_0.5_5_100.0.pickle', "rb") as f:
    KE100 = np.real(pickle.load(f))

# Load imshow data for plots (c) and (d)
with open(directory_path + '/EF_1.0_0.5_5_7.0.pickle', "rb") as f:
    EF7 = np.real(pickle.load(f))
with open(directory_path + '/EF_1.0_0.5_5_100.0.pickle', "rb") as f:
    EF100 = np.real(pickle.load(f))

# Create a figure with four subplots (2 rows, 2 columns)
fig, ax = plt.subplots(2, 2, figsize=(12, 10), gridspec_kw={'wspace': 0.3, 'hspace': 0.3})
t_vals = np.linspace(1 * 0.05, 2000 * 0.05, 2000)

# Plot (a) in the top-left subplot
ax[0, 0].plot(t_vals, KE7, label='$T = 7$, $\Delta K_{max} = $' + f' ${np.max(KE7):.3f}$')
ax[0, 0].plot(t_vals, KE100, label='$T = 100$, $\Delta K_{max} = $' + f' ${np.max(KE100):.3f}$', linestyle='--')
ax[0, 0].legend(fontsize=16, loc='upper right')
ax[0, 0].set_xlabel('$t$', fontsize=20)
ax[0, 0].set_ylabel('$\Delta K$', fontsize=20)
ax[0, 0].tick_params(axis='both', direction='in', labelsize=20)
ax[0, 0].text(0.5, 1.07, r'$(a)$', fontsize=20, transform=ax[0, 0].transAxes, ha='center', va='center')

# Plot (b) in the top-right subplot
ax[0, 1].plot(t_vals, pn7, label='$T = 7$, $\mathcal{P}(t = 100) = $' + f' ${pn7[-1]:.3f}$')
ax[0, 1].plot(t_vals, pn100, label='$T = 100$, $\mathcal{P}(t = 100) = $' + f' ${pn100[-1]:.3f}$', linestyle='--')
ax[0, 1].legend(fontsize=16, loc='lower right')
ax[0, 1].set_xlabel('$t$', fontsize=20)
ax[0, 1].set_ylabel('$\mathcal{P}$', fontsize=20)
ax[0, 1].tick_params(axis='both', direction='in', labelsize=20)
ax[0, 1].text(0.5, 1.07, r'$(b)$', fontsize=20, transform=ax[0, 1].transAxes, ha='center', va='center')

# Plot (c) in the bottom-left subplot with time values on the x-axis
im_c = ax[1, 0].imshow(EF7, aspect='auto', origin='lower', cmap='jet',
                       extent=[t_vals[0], t_vals[-1], 0, EF7.shape[0]])
ax[1, 0].set_xlabel('$t$', fontsize=20)
ax[1, 0].set_ylabel('$L$', fontsize=20)
ax[1, 0].tick_params(axis='both', labelsize=20)
ax[1, 0].text(0.5, 1.07, r'$(c)$', fontsize=20, transform=ax[1, 0].transAxes, ha='center', va='center')
# Add colorbar for (c)
# cbar_c = fig.colorbar(im_c, ax=ax[1, 0], fraction=0.046, pad=0.04)
# cbar_c.ax.tick_params(labelsize=16)

# Plot (d) in the bottom-right subplot with time values on the x-axis
im_d = ax[1, 1].imshow(EF100, aspect='auto', origin='lower', cmap='jet',
                       extent=[t_vals[0], t_vals[-1], 0, EF100.shape[0]])
ax[1, 1].set_xlabel('$t$', fontsize=20)
ax[1, 1].set_ylabel('$L$', fontsize=20)
ax[1, 1].tick_params(axis='both', labelsize=20)
ax[1, 1].text(0.5, 1.07, r'$(d)$', fontsize=20, transform=ax[1, 1].transAxes, ha='center', va='center')
# Add colorbar for (d)
cbar_d = fig.colorbar(im_d, ax=ax[1, 1], fraction=0.046, pad=0.04)
cbar_d.ax.tick_params(labelsize=20)
cbar_d.set_label(r'$\Delta F$', fontsize=20)

# Adjust layout for better spacing
plt.tight_layout()
plt.savefig('data_50.pdf', bbox_inches='tight', dpi=1200)
# plt.show()
