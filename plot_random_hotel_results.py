""" This script plots the results of the execution of the large number of cases with randomized hotel placement.
    While many cases are simulated, some result infeasible solutions. 
    The sample size which is considered is given by the parameter 'sample_size'.
"""
import os, shutil
import pickle
import pandas as pd
import numpy as np

import utils

file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results_random_hotel/random_uniform')
plot_dir = os.path.join(result_dir, f'plots')
utils.create_directory(plot_dir)

# define the nominal parameters of the scenario
sci_case_name = 'Bioassessment_24_scival'   # Bioassessment_24_scival, Fracking_24_scival, Plume_24_scival, RivNetworkContinuity_24_scival

N = 24
H = 3
D = 4
T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = -1 # timeout in seconds. Enter negative number if optimal value is desired.
num_of_iterations = 250

sample_size = 200

result_parent_folder_path = os.path.join(result_dir, f'{sci_case_name}')

optimal_value_list = []
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

print(f'Implemented sample size: {counter}')
#print(optimal_value_list)
utils.plt.figure(figsize=(4, 3))
utils.plt.hist(optimal_value_list, bins=30, edgecolor='black', alpha=0.7)

#utils.plt.title(f'{science_case_label}')
#utils.plt.xlim(0,1)
utils.plt.xlabel('Score')
utils.plt.ylabel('Frequency')

utils.plt.tight_layout() # Adjust layout to prevent cutting off labels

plot_fp = os.path.join(plot_dir, f'{sci_case_name}_scores.{utils.figure_format}')
utils.plt.savefig(plot_fp , format=utils.figure_format)
