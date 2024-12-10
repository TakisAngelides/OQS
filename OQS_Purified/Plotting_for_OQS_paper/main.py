import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
import seaborn as sns
import pickle
from matplotlib.ticker import LogFormatterMathtext
from matplotlib import rc 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('text', usetex=True)

def plot_ef_symmetry():
            
    f_inputs = h5py.File("/Users/takisangelides/Downloads/Plotting_for_OQS_Paper/inputs.h5", 'r')

    for file in os.listdir("/Users/takisangelides/Downloads/Plotting_for_OQS_Paper/HDF5"):
        
        # try:
        
        g = f_inputs[file[:-3]]
        attributes_dict = {attr_name: attr_value for attr_name, attr_value in g.attrs.items()}
        
        if attributes_dict['wis'] == 'dirac_vacuum_with_string':

            # With string
            time_step_limit = -1
            N = attributes_dict['N']
            if N < 50:
                continue
            l_0_1 = attributes_dict['l_0_1']
            ma, aD, aT, cqns, cutoff, l_0_1, teo, waf, x_val = attributes_dict['ma'], attributes_dict['aD'], attributes_dict['aT'], attributes_dict['cqns'], attributes_dict['cutoff'], attributes_dict['l_0_1'], attributes_dict['teo'], attributes_dict['waf'], np.round(attributes_dict['x'], decimals = 3)
            tau = attributes_dict['tau']
            if tau > 0.002:
                continue
            staggering = np.array([(-1)**n for n in range(N)])
            f = h5py.File(f'{"/Users/takisangelides/Downloads/Plotting_for_OQS_Paper/HDF5"}/{file}', 'r')
            try:
                z_configs = np.asarray(f['z_configs'])
            except:
                continue
            q_configs = np.array([0.5*(np.real(z_configs[:, i]) + staggering) for i in range(z_configs.shape[1])])
            ef_configs = np.transpose(np.array([np.array([l_0_1 + sum(q_configs[i][0:j + 1]) for j in range(q_configs.shape[1] - 1)]) for i in range(q_configs.shape[0])]))
            f.close()
        
            # Without string
            file_without_string = f'{int(file[:-3])-1}'
            attributes_dict_without_string = {attr_name: attr_value for attr_name, attr_value in f_inputs[file_without_string].attrs.items()}

            N = attributes_dict_without_string['N']
            l_0_1 = attributes_dict_without_string['l_0_1']
            staggering = np.array([(-1)**n for n in range(N)])
            f_without_string = h5py.File(f'{"/Users/takisangelides/Downloads/Plotting_for_OQS_Paper/HDF5"}/{file_without_string}.h5', 'r')
            try:
                z_configs_without_string = np.asarray(f_without_string['z_configs'])
            except:
                continue
            q_configs_without_string = np.array([0.5*(np.real(z_configs_without_string[:, i]) + staggering) for i in range(z_configs_without_string.shape[1])])
            ef_configs_without_string = np.transpose(np.array([np.array([l_0_1 + sum(q_configs_without_string[i][0:j + 1]) for j in range(q_configs_without_string.shape[1] - 1)]) for i in range(q_configs_without_string.shape[0])]))
            f_without_string.close()
            
            z = ef_configs - ef_configs_without_string
            t_over_a_list = [0] + list(tau*(np.arange(1, z.shape[1])))
            x = np.round(t_over_a_list, decimals = 3)
            y = list(np.arange(1, N))
            
            # Assuming t_over_a_list, z, and N are defined
            fig = plt.figure(figsize=(10, 8))

            # Top subplot (full data)
            ax_top = fig.add_subplot(3, 1, 1)
            im = ax_top.imshow(z, aspect='auto', origin='lower', extent=[t_over_a_list[0], t_over_a_list[-1], 0, z.shape[0]], cmap='jet')
            ax_top.set_ylabel(r"$L$", fontsize=24)
            ax_top.tick_params(axis='both', labelsize=24)
            ax_top.set_xticklabels([])
            for spine in ax_top.spines.values():
                spine.set_visible(False)

            # Adding custom ticks for `ax_top`
            full_links = [0, z.shape[0]//2, z.shape[0] - 1]  # Start, middle, end links
            ax_top.set_yticks(full_links)
            # Use LaTeX formatting for the labels
            ax_top.set_yticklabels([r"${}$".format(tick) for tick in full_links], fontsize=24)

            # Middle subplot (middle 20 links)
            middle_links_start = z.shape[0] // 2 - 5
            middle_links_end = z.shape[0] // 2 + 5
            ax_middle = fig.add_subplot(3, 1, 2)
            im_middle = ax_middle.imshow(z[middle_links_start:middle_links_end, :], aspect='auto', origin='lower',
                                        extent=[t_over_a_list[0], t_over_a_list[-1], middle_links_start, middle_links_end], cmap='jet')
            ax_middle.set_ylabel(r"$L$", fontsize=24)
            ax_middle.tick_params(axis='both', labelsize=24)
            for spine in ax_middle.spines.values():
                spine.set_visible(False)

            # Adding custom ticks for `ax_middle`
            middle_links = [middle_links_start, (middle_links_start + middle_links_end) // 2, middle_links_end]  # Start, middle, end links
            ax_middle.set_yticks(middle_links)
            ax_middle.set_xticklabels([])
            # Use LaTeX formatting for the labels
            ax_middle.set_yticklabels([r"${}$".format(tick) for tick in middle_links], fontsize=24)

            # Bottom subplot (plot data as per the original code)
            ax_bottom = fig.add_subplot(3, 1, 3)
            alpha_list = [(i / (N//2 - 2))**4 for i in range(N//2 - 1)]
            for i in range(N//2 - 1):
                ax_bottom.plot(t_over_a_list, np.abs(z[i, :] - z[N-1-i-1, :]), label=f'$i$ = {i}', alpha=alpha_list[i], color='k')
            ax_bottom.set_xlabel(r"$t$", fontsize=24)
            ax_bottom.set_xlim(0, 3.5)
            ax_bottom.set_ylabel(r"$P$", fontsize=24)
            ax_bottom.set_yscale('log')
            ax_bottom.tick_params(axis='both', labelsize=24)

            # Adding custom y-ticks for `ax_bottom`
            # Choose ticks based on data range; here, we add two more ticks for better readability
            y_ticks_bottom = [10**-14, 10**-9, 10**-4]  # Example tick values; adjust according to your data range
            ax_bottom.set_yticks(y_ticks_bottom)
            ax_bottom.set_yticklabels([f"{tick:.1e}" for tick in y_ticks_bottom], fontsize=24)
            ax_bottom.yaxis.set_major_formatter(LogFormatterMathtext())

            plt.tight_layout()

            # Single color bar for the top two subplots
            # fig.subplots_adjust(hspace=0.4)
            cbar = fig.colorbar(im, ax=[ax_top, ax_middle], location='top', fraction=0.05)
            cbar.ax.tick_params(labelsize=20)
            cbar.set_label(r"$\Delta F$", fontsize=24)
            for spine in cbar.ax.spines.values():
                spine.set_visible(False)

            # Adjust layout and annotations
            fig.text(0.99, 0.79, r'$(a)$', ha='center', fontsize=20)
            fig.text(0.99, 0.52, r'$(b)$', ha='center', fontsize=20)
            fig.text(0.99, 0.23, r'$(c)$', ha='center', fontsize=20)

            # Save figure
            plt.savefig(f'ef_symmetry_{tau}.pdf', bbox_inches='tight', dpi=1200)
            # plt.savefig(f'ef_symmetry_{tau}.png', bbox_inches='tight')
            
            # with open(f'Plots/data/EF_{N}_{ma}_{l_0_1}_{aD}_{tau}.pickle', 'wb') as f:
            #     pickle.dump(z, f)
                
            
                            
        #     else:
                
        #         continue
            
        # except:
            
        #     print(file)

plot_ef_symmetry()
