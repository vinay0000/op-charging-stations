""" This script collects the results of execution of number of different scenarios and generates insights.
"""
import os
import pickle
import pandas as pd
import numpy as np

file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results')

T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = 2*3600 # timeout in seconds. Enter negative number if optimal value is desired.

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

#### Parameter values ####
science_case_list = [('Bioassessment_24_scival', 24), ('Bioassessment_50_scival', 50), ('Bioassessment_103_scival', 103), 
                        ('Fracking_24_scival',24), ('Fracking_50_scival',50), ('Fracking_103_scival',103),
                        ('Plume_30_scival',30), ('Plume_58_scival',58), ('Plume_110_scival',110),
                        ('RivNetworkContinuity_24_scival',24), ('RivNetworkContinuity_50_scival',50), ('RivNetworkContinuity_103_scival',103)]

D_range = range(1,14) # number of trips
H_range = range(2,4) # number of hotels

# search the 'results' folder for available results and consolidate all of them into a pandas dataframe
df = pd.DataFrame()
for sci_idx, science_case in enumerate(science_case_list):
    
    sci_case_name =  science_case[0]
    N = science_case[1]
    sci_raw_data_fp = os.path.join(data_dir, f'{sci_case_name}.csv') # Absolute path to the raw science data file
    result_folder_path = os.path.join(result_dir, f'{sci_case_name}') # Absolute folder path where the results are to be written

    for D in D_range:
        for H in H_range:
            result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl') 
            if os.path.isfile(result_fp):
                print(result_fp)
                (_H, _N, _D, _T_Max, _T_CH, _c_pos, _Si, _score, transitions, halt_times, 
                    max_flight_times, flight_times, subtour_u, optimal_value, gap, 
                    nconss, nvars, optimal_sol, process_time) = load_results(result_fp)
                
                if not (_H==H and _N==N):
                    raise RuntimeError(f'unexpected results in pickle file. H{H} vs _H{_H}')
                
                # Creating a dictionary to hold all the variables as columns
                data = {
                    'scence_case': [sci_case_name],
                    'H': [_H],
                    'N': [_N],
                    'D': [_D],
                    'T_Max': [_T_Max],
                    'T_CH': [_T_CH],
                    #'c_pos': [c_pos],
                    #'Si': [Si],
                    #'score': [score],
                    #'transitions': [transitions],
                    #'halt_times': [halt_times],
                    #'max_flight_times': [max_flight_times],
                    #'flight_times': [flight_times],
                    #'subtour_u': [subtour_u],
                    'optimal_value': [optimal_value],
                    'gap': [gap],
                    'nconss': [nconss],
                    'nvars': [nvars],
                    'optimal_sol': [optimal_sol],
                    'process_time': [process_time]
                }

                # Creating the pandas DataFrame
                new_row = pd.DataFrame(data)                

                # Adding the new row using pd.concat
                df = pd.concat([df, new_row], ignore_index=True)
                
print(df)