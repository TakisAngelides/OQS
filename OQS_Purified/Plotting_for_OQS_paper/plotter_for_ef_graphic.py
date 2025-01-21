import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc 
# rc('text', usetex=True)

# Parameters
base_paths = {
    12: "/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_47"
}
masses = [0.5]
l_0_values = [0.0]
aD_values = [2.0]
tau = 0.01
steps = 10000
t_over_a_list = [0] + list(tau * np.arange(1, steps))

# Collect all EF and Q values to determine the min and max for color scale
all_ef_values = []
all_q_values = []

for N in base_paths:
    for mass in masses:
        for l_0 in l_0_values:
            for aD in aD_values:
                
                ef_file_name = f"EF_{mass}_{l_0}_{aD}.pickle"
                q_file_name = f"Q_{mass}_{l_0}_{aD}.pickle"  # Assuming Q data is stored similarly
                
                ef_file_path = os.path.join(base_paths[N], ef_file_name)
                q_file_path = os.path.join(base_paths[N], q_file_name)
                
                # Load EF data
                try:
                    with open(ef_file_path, "rb") as f:
                        ef_values = pickle.load(f)
                        all_ef_values.append(ef_values)
                except FileNotFoundError:
                    print(f"File not found: {ef_file_path}")
                    continue
                
                # Load Q data
                try:
                    with open(q_file_path, "rb") as f:
                        q_values = pickle.load(f)
                        all_q_values.append(q_values)
                except FileNotFoundError:
                    print(f"File not found: {q_file_path}")
                    continue

# Plot setup with constrained layout
fig, axs = plt.subplots(1, 2, figsize=(30, 12), constrained_layout=True)

for i, N in enumerate([12]):
    data_path = base_paths[N]
    for j, mass in enumerate(masses):
        for k, l_0 in enumerate(l_0_values):
            for l, aD in enumerate(aD_values):
                ef_file_name = f"EF_{mass}_{l_0}_{aD}.pickle"
                q_file_name = f"Q_{mass}_{l_0}_{aD}.pickle"
                
                ef_file_path = os.path.join(data_path, ef_file_name)
                q_file_path = os.path.join(data_path, q_file_name)
                
                try:
                    with open(ef_file_path, "rb") as f:
                        ef_values = pickle.load(f)
                except FileNotFoundError:
                    continue
                
                try:
                    with open(q_file_path, "rb") as f:
                        q_values = pickle.load(f)
                except FileNotFoundError:
                    continue
                
                # EF Heatmap
                tt = 3000
                im1 = axs[0].imshow(ef_values[:, :tt], aspect='auto', origin='lower',
                                    extent=[t_over_a_list[0], t_over_a_list[tt], 0, ef_values.shape[0]],
                                    vmin=np.min(ef_values[:, :tt]), vmax=np.max(ef_values[:, :tt]), cmap='jet')
                
                axs[0].set_ylabel(r'$n$', fontsize=24)
                axs[0].tick_params(axis='y', labelsize=24)
                axs[0].set_yticks(range(0, 11))
                
                axs[0].set_xlabel(r"$t$", fontsize=24)
                axs[0].tick_params(axis='x', labelsize=24)
                
                # Add colorbar for EF
                cbar1 = fig.colorbar(im1, ax=axs[0], orientation='horizontal', location='top', fraction=0.03)
                cbar1.ax.tick_params(labelsize=24)
                cbar1.set_label(r'$(a): $' + r'$\Delta F(n)$', fontsize=24)
                
                # Q Heatmap
                tt = 700
                im2 = axs[1].imshow(q_values[:, :tt], aspect='auto', origin='lower',
                                    extent=[t_over_a_list[0], t_over_a_list[tt], 0, q_values.shape[0]],
                                    vmin=np.min(q_values[:, :tt]), vmax=np.max(q_values[:, :tt]), cmap='jet')
                
                axs[1].set_ylabel(r'$n$', fontsize=24)
                axs[1].tick_params(axis='y', labelsize=24)
                axs[1].set_yticks(range(0, 12))
                
                axs[1].set_xlabel(r"$t$", fontsize=24)
                axs[1].tick_params(axis='x', labelsize=24)
                
                # Add colorbar for Q
                cbar2 = fig.colorbar(im2, ax=axs[1], orientation='horizontal', location='top', fraction=0.03)
                cbar2.ax.tick_params(labelsize=24)
                cbar2.set_label(r'$(b): $' + r'$\Delta Q(n)$', fontsize=24)
                
                for spine in cbar1.ax.spines.values():
                    spine.set_visible(False)
                for spine in cbar2.ax.spines.values():
                    spine.set_visible(False)
                for spine in axs[0].spines.values():
                    spine.set_visible(False)
                for spine in axs[1].spines.values():
                    spine.set_visible(False)


plt.savefig('ef_q_N_12_graphic.pdf', bbox_inches='tight', dpi=1200)
