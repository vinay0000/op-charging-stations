""" This script collects the results of execution of number of different scenarios and generates insights.
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

def create_directory(dir_name):
    # Check if the directory already exists
    if os.path.exists(dir_name):
        # Remove the existing directory
        shutil.rmtree(dir_name)
        print(f"Removed existing directory: {dir_name}")

    # Create a new directory
    os.makedirs(dir_name)
    print(f"Created new directory: {dir_name}")

file_dir = os.path.dirname(os.path.abspath(__file__))
#data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results_for_paper_set2')

plot_dir = os.path.join(result_dir, f'plots')
create_directory(plot_dir)
# (re)make plot directory


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
                        ('Plume_24_scival',24), ('Plume_50_scival',50), ('Plume_103_scival',103),
                        ('RivNetworkContinuity_24_scival',24), ('RivNetworkContinuity_50_scival',50), ('RivNetworkContinuity_103_scival',103)]


D_range = range(1,14) # number of trips
H_range = range(2,5) # number of hotels


# search the 'results' folder for available results and consolidate all of them into a pandas dataframe
df = pd.DataFrame()
for sci_idx, science_case in enumerate(science_case_list):
    
    sci_case_name =  science_case[0]
    N = science_case[1]
    result_folder_path = os.path.join(result_dir, f'{sci_case_name}') # Absolute folder path where the results are to be written

    for D in D_range:
        for H in H_range:
            result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl') 
            if os.path.isfile(result_fp):
                #print(result_fp)
                (_H, _N, _D, _T_Max, _T_CH, c_pos, Si, score, transitions, halt_times, 
                    max_flight_times, flight_times, subtour_u, optimal_value, gap, 
                    nconss, nvars, optimal_sol, process_time) = load_results(result_fp)
                
                if not (_H==H and _N==N and _D==D and _T_Max==T_Max and _T_CH==T_CH):
                    raise RuntimeError(f'Unexpected data in pickle file. Check: H{H} vs _H{_H}, N{N} vs _N{_N}, _D{_D} vs D{D} and _T_Max{_T_Max} vs T_Max{T_Max} and _T_CH{_T_CH} vs T_CH{T_CH}')
                
                try:
                    # Creating a dictionary to hold all the variables as columns
                    data = {
                        'science_case': [sci_case_name],
                        'H': [_H],
                        'N': [_N],
                        'D': [_D],
                        'T_Max': [_T_Max],
                        'T_CH': [_T_CH],
                        'c_pos': [c_pos],
                        'Si': [Si],
                        'score': [score],
                        'transitions': [transitions],
                        'halt_times': [halt_times],
                        'max_flight_times': [max_flight_times],
                        'flight_times': [flight_times],
                        'subtour_u': [subtour_u],
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

                except Exception as e:
                    print(f"Error loading file {result_fp}: {e}")
                
#print(df)
# Create plots for each science case
unique_cases = df['science_case'].unique()

for case in unique_cases:

    case_df = df[df['science_case'] == case]

    if 'Bioassessment' in case:
        case_label = 'Bioassessment'
    elif 'Fracking' in case:
        case_label = 'Fracking'
    elif 'Plume' in case:
        case_label = 'Low flow estimation'
    elif 'RivNetworkContinuity' in case:
        case_label = 'River network continuity'
    else:
        raise RuntimeError(f'Unknown case {case}')

    plt.figure(figsize=(4, 3))
    
    # Create a plot for each unique value of H
    for h_value in case_df['H'].unique():
        h_df = case_df[case_df['H'] == h_value]
        plt.plot(h_df['D'], h_df['optimal_value'], marker='o', label=f'H={h_value}')
    
    plt.title(f'{case_label}')
    plt.xlabel('D (Number of Trips)')
    plt.ylabel('Optimal Value')
    plt.legend()
    plt.grid()
    plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{case}_optimal_value_vs_D.{figure_format}')
    plt.savefig(plot_fp, format=figure_format)
    plt.close()

    # Create a plot for the sum of halt_times vs D
    plt.figure(figsize=(4, 3))
    for h_value in case_df['H'].unique():
        h_df = case_df[case_df['H'] == h_value]
        # Sum the arrays in halt_times for each row and store the results
        halt_time_sums = 1/60 * h_df['halt_times'].apply(lambda x: np.sum(x))  # Sum each array
        plt.plot(h_df['D'], halt_time_sums, label=f'H={h_value}', marker='x')
    
    plt.title(f'{case_label}')
    plt.xlabel('D (Number of Trips)')
    plt.ylabel('Total Charging Time [minutes]')
    plt.legend()
    plt.grid()
    plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{case}_charging_time_vs_D.{figure_format}')
    plt.savefig(plot_fp, format=figure_format)
    plt.close()


    # Create a plot for the total tour time vs D
    plt.figure(figsize=(4, 3))
    for h_value in case_df['H'].unique():
        h_df = case_df[case_df['H'] == h_value]
        # Sum the arrays in halt_times for each row and store the results
        tour_time = 1/60 * (h_df['flight_times'].apply(lambda x: np.sum(x)) +  h_df['halt_times'].apply(lambda x: np.sum(x)))  # Sum each array
        plt.plot(h_df['D'], tour_time, label=f'H={h_value}', marker='x')  # Use scatter for clearer visualization
    
    plt.title(f'{case_label}')
    plt.xlabel('D (Number of Trips)')
    plt.ylabel('Tour Time [minutes]')
    plt.legend()
    plt.grid()
    plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{case}_tour_time_vs_D.{figure_format}')
    plt.savefig(plot_fp, format=figure_format)
    plt.close()


#### DEBUG: To verify that the hotels generated by K-means is the same and in the same order across all cases. Simultaneously the vertex position data is also verified. ####

def read_data_file(input_fp):
        """
        Reads input data file and extracts vertex locations and scores.
        Assumes the file contains 'x y Si' data where:
            x  = x coordinate
            y  = y coordinate
            Si = score

        Returns:
            location (list): List of (x, y) coordinates.
            Si (list): List of corresponding scores.
        """
        if not os.path.exists(input_fp):
            raise FileNotFoundError(f"Input file '{input_fp}' not found.")

        location = []
        Si = []

        with open(input_fp, 'r') as file:
            next(file)  # Skip the first line (header)
            for line in file:
                values = list(map(float, line.split()))
                x, y, score = values[-3], values[-2], values[-1]
                location.append([x, y])
                Si.append(score)

        return location, Si


for sci_idx, science_case in enumerate(science_case_list):
    
    sci_case_name =  science_case[0]
    N = science_case[1]
    result_folder_path = os.path.join(result_dir, f'{sci_case_name}') # Absolute folder path where the results are to be written

    for D in D_range:
        for H in H_range:
            result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl') 
            if os.path.isfile(result_fp):
                #print(result_fp)
                (_H, _N, _D, _T_Max, _T_CH, c_pos, Si, score, transitions, halt_times, 
                    max_flight_times, flight_times, subtour_u, optimal_value, gap, 
                    nconss, nvars, optimal_sol, process_time) = load_results(result_fp)
                
                # load the vertex data
                vertex_data_fp = os.path.join(result_folder_path, 'vertex_data.opc')

                # load the hotel position data
                hotel_fp = os.path.join(result_folder_path, f'hotels_H{H}.opc')

                c_pos_hotel_verify, Si_hotel_verify = read_data_file(hotel_fp)
                c_pos_verify, Si_verify = read_data_file(vertex_data_fp)

                c_pos_verify = c_pos_hotel_verify + c_pos_verify
                Si_verify = Si_hotel_verify + Si_verify

                if c_pos_verify == c_pos:
                    print("Hotel and vertex verification succeeded.")
                else:
                    raise RuntimeError("Hotel and vertex verification failed.")






                

