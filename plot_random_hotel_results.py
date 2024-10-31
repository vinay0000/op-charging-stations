""" This script plots the results of the execution of the large number of cases with randomized hotel placement.
"""
import os, shutil
import pickle
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 10,           # Set font size to 10
    'font.family': 'serif',    # Use serif fonts
    'font.serif': ['Times New Roman']  # Set font to Times New Roman
})
figure_format = 'pdf'

def load_results(mip_results_fp):
    """
    Load result data from a pickle file.
    
    Parameters:
    -----------
    mip_results_fp : str
        Absolute file path to the pickle file which contains the results of the optimization.
    
    Returns:
    --------
    tuple : All loaded results from the pickle file.
    """
    # Ensure file exists
    if not os.path.exists(mip_results_fp):
        raise FileNotFoundError(f"Results file '{mip_results_fp}' not found.")
    
    with open(mip_results_fp, 'rb') as file:
        H = pickle.load(file)
        N = pickle.load(file)
        D = pickle.load(file)
        T_Max = pickle.load(file)
        T_CH = pickle.load(file)
        c_pos = pickle.load(file)
        Si = pickle.load(file)
        score = pickle.load(file)
        transitions = pickle.load(file)
        halt_times = pickle.load(file)
        max_flight_times = pickle.load(file)
        flight_times = pickle.load(file)
        subtour_u = pickle.load(file)
        subtour_u = np.array(subtour_u) - min(np.array(subtour_u))
        optimal_value = pickle.load(file)
        gap = pickle.load(file)
        nconss = pickle.load(file)
        nvars = pickle.load(file)
        optimal_sol = pickle.load(file)
        process_time = pickle.load(file)

    return (H, N, D, T_Max, T_CH, c_pos, Si, score, transitions, halt_times, 
            max_flight_times, flight_times, subtour_u, optimal_value, gap, 
            nconss, nvars, optimal_sol, process_time)

file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results_random_hotel/random_uniform')

# define the nominal parameters of the scenario
sci_case_name = 'RivNetworkContinuity_24_scival'   # Bioassessment_24_scival, Fracking_24_scival, Plume_24_scival, RivNetworkContinuity_24_scival

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

result_parent_folder_path = os.path.join(result_dir, f'{sci_case_name}')

optimal_value_list = []
counter = 0
for iter_n in range(num_of_iterations):
    
    result_folder_path = os.path.join(result_parent_folder_path, f'{iter_n}') # Absolute folder path where the results are present

    result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl')
    try:
        (_H, _N, _D, _T_Max, _T_CH, c_pos, Si, score, transitions, halt_times, 
                                max_flight_times, flight_times, subtour_u, optimal_value, gap, 
                                nconss, nvars, optimal_sol, process_time) = load_results(result_fp)
        optimal_value_list.append(optimal_value)
        counter = counter + 1
    except:
        pass

print(f'Valid number of data points: {counter}')
#print(optimal_value_list)
plt.figure(figsize=(4, 3))
plt.hist(optimal_value_list, bins=30, edgecolor='black', alpha=0.7)

#plt.title(f'{science_case_label}')
#plt.xlim(0,1)
plt.xlabel('Score')
plt.ylabel('Frequency')

plt.tight_layout() # Adjust layout to prevent cutting off labels

#plot_fp = os.path.join(plot_dir, f'{sci_case_name}_scores.{figure_format}')
#plt.savefig(plot_fp , format=figure_format)
plt.show()