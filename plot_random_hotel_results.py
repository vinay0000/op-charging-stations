""" This script plots the results of the execution of the large number of cases with randomized hotel placement.
    While many cases are simulated, some result in infeasible solutions. 
    The sample size which is considered is given by the parameter 'sample_size'.
"""
import os, shutil
import pickle
import pandas as pd
import numpy as np

import utils

file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results_random_hotel')
r_unif_result_dir = os.path.join(result_dir, f'random_uniform')
kmeans_result_dir = os.path.join(result_dir, f'kmeans')

r_unif_plot_dir = os.path.join(r_unif_result_dir, f'plots')
utils.create_directory(r_unif_plot_dir)

kmeans_plot_dir = os.path.join(kmeans_result_dir, f'plots')
utils.create_directory(kmeans_plot_dir)

combined_plot_dir = os.path.join(result_dir, f'plots') # for stacked histogram plots
utils.create_directory(combined_plot_dir)

######## define the nominal parameters of the scenario ########
sci_case_name = 'Plume_24_scival'   # Bioassessment_24_scival, Fracking_24_scival, Plume_24_scival, RivNetworkContinuity_24_scival

N = 24
H = 3
D = 4
T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = -1 # timeout in seconds. Enter negative number if optimal value is desired.
num_of_iterations = 300

sample_size = 200

def save_histogram_plot(optimal_value_list, plot_dir, max_objective):

    #print(optimal_value_list)
    utils.plt.figure(figsize=(4, 3))
    utils.plt.hist(optimal_value_list, bins=30, edgecolor='black', alpha=0.7)

    #utils.plt.title(f'{science_case_label}')
    utils.plt.xlim(0,max_objective)
    utils.plt.ylim(0,50)
    utils.plt.xlabel('Optimal Value')
    utils.plt.ylabel('Frequency')

    utils.plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{sci_case_name}_optval_dist.{utils.figure_format}')
    utils.plt.savefig(plot_fp , format=utils.figure_format)


def collate_results(result_parent_folder_path, plot_dir):
    
    optimal_value_list = []
    Si = None
    counter = 0
    for iter_n in range(num_of_iterations):
        
        result_folder_path = os.path.join(result_parent_folder_path, f'{iter_n}') # Absolute folder path where the results are present

        result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl')
        try:
            (_H, _N, _D, _T_Max, _T_CH, c_pos, Si, score, transitions, halt_times, 
                                    max_flight_times, flight_times, subtour_u, optimal_value, gap, 
                                    nconss, nvars, optimal_sol, process_time) = utils.load_results(result_fp)
            optimal_value_list.append(optimal_value)
            counter = counter + 1
            if counter == sample_size:
                break
        except:
            pass
    
    max_objective = sum(Si)
    save_histogram_plot(optimal_value_list, plot_dir, max_objective)

    return (optimal_value_list, counter, max_objective)

# kmeans
(kmeans_optimal_val, kmeans_counter, max_objective) = collate_results(os.path.join(kmeans_result_dir, f'{sci_case_name}'), kmeans_plot_dir)
print(f'Implemented kmeans sample size: {kmeans_counter}')

# Random uniform
(r_unif_optimal_val, r_unif_counter, max_objective) = collate_results(os.path.join(r_unif_result_dir, f'{sci_case_name}'), r_unif_plot_dir)
print(f'Implemented random-uniform sample size: {r_unif_counter}')

# Plot both histograms on same plot area
fig, (ax1, ax2) = utils.plt.subplots(2, 1, sharex=True, figsize=(4, 3))
# Plot the first histogram on the first subplot
ax1.hist(r_unif_optimal_val, bins=30, color='skyblue', edgecolor='black')
ax1.set_ylabel('Frequency')
ax1.text(0.1, 0.85, 'Random-uniform', transform=ax1.transAxes, fontsize=10, color='blue')

#ax1.set_title('Separate Histograms with Shared X-Axis')

# Plot the second histogram on the second subplot
ax2.hist(kmeans_optimal_val, bins=30, color='lightgreen', edgecolor='black')
ax2.set_xlabel('Optimal Value')
ax2.set_ylabel('Frequency')
ax2.text(0.1, 0.85, 'K-Means', transform=ax2.transAxes, fontsize=10, color='green')

# Set the x-axis and y-axis limits
ax1.set_xlim(0, max_objective)

# Adjust layout to prevent overlap
utils.plt.tight_layout()

# Show the plot
plot_fp = os.path.join(combined_plot_dir, f'{sci_case_name}__H{H}_D{D}_optval_dist.{utils.figure_format}')
utils.plt.savefig(plot_fp , format=utils.figure_format)
