""" This script collects the results of execution of number of different scenarios and generates insights.
"""
import os, shutil
import pickle
import pandas as pd
import numpy as np

import utils

# List of colors to cycle through
colors = ['blue', 'orange', 'green', 'purple', 'pink']

file_dir = os.path.dirname(os.path.abspath(__file__))
result_dir = os.path.join(file_dir, f'results_set0')

# (re)make plot directory
plot_dir = os.path.join(result_dir, f'plots')
utils.create_directory(plot_dir)


T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = 2*3600 # timeout in seconds. Enter negative number if optimal value is desired.


#### Parameter values ####
science_case_list = [('Bioassessment_24_scival', 24), ('Bioassessment_50_scival', 50), ('Bioassessment_103_scival', 103), 
                        ('Fracking_24_scival',24), ('Fracking_50_scival',50), ('Fracking_103_scival',103),
                        ('Plume_24_scival',24), ('Plume_50_scival',50), ('Plume_103_scival',103),
                        ('RivNetworkContinuity_24_scival',24), ('RivNetworkContinuity_50_scival',50), ('RivNetworkContinuity_103_scival',103)]


D_range = range(1,27) # number of trips
H_range = range(2,7) # number of hotels

#science_case_list = [('Bioassessment_24_scival', 24)]

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
                    nconss, nvars, optimal_sol, process_time) = utils.load_results(result_fp)
                
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

    fig, ax1 = utils.plt.subplots(figsize=(4.5, 3.5))
    ax2 = ax1.twinx()

    # Create a plot for each unique value of H
    for i, h_value in enumerate(case_df['H'].unique()):
        h_df = case_df[case_df['H'] == h_value]
         # Plot optimal values on the primary y-axis (left)
        ax1.plot(h_df['D'], h_df['optimal_value'], color=colors[i % len(colors)], marker='o', markersize=2, linewidth=1, label=f'H={h_value}')
        
        # Plot process time on the secondary y-axis (right)
        ax2.plot(h_df['D'], h_df['process_time']/60, ':', color=colors[i % len(colors)], marker='x', markersize=3, linewidth=1, label=f'Process Time (H={h_value})')
            
    # Label axes and title
    #ax1.set_title(f'{case_label}')
    sum_of_scores = sum(np.array(h_df['Si'].iloc[0]))
    ax1.axhline(y=sum_of_scores, color='r', linestyle='-', linewidth=1) # Add a horizontal line at the maximum possible Score
    ax1.text(x=0.15, y=sum_of_scores - 1,  # Position slightly below the line
         s=f'{sum_of_scores:.2f}',  # Text content
         color='r', fontsize=10, ha='center', va='bottom', 
         transform=ax1.get_yaxis_transform())  # Relative to y-axis position
    
    ax1.set_xlabel('D (Number of Trips)')
    ax1.set_ylabel('Optimal Value')
    ax2.set_ylabel('Process Time [minutes]')

    # Adjust tick color for clarity
    ax1.tick_params(axis="y")
    ax2.tick_params(axis="y")

    # Handle legends for each axis separately
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1, labels_1, loc='center right')


    # Set grid and layout adjustments
    utils.plt.grid(True)
    utils.plt.tight_layout()  # Adjust layout to prevent cutting off labels

    # Save the plot
    plot_fp = os.path.join(plot_dir, f'{case}_optimal_value_vs_D.{utils.figure_format}')
    utils.plt.savefig(plot_fp, format=utils.figure_format)
    utils.plt.close()

    # Create a plot for the sum of halt_times vs D
    utils.plt.figure(figsize=(4, 3))
    for i, h_value in enumerate(case_df['H'].unique()):
        h_df = case_df[case_df['H'] == h_value]
        # Sum the arrays in halt_times for each row and store the results
        halt_time_sums = 1/60 * h_df['halt_times'].apply(lambda x: np.sum(x))  # Sum each array
        utils.plt.plot(h_df['D'], halt_time_sums, label=f'H={h_value}', color=colors[i % len(colors)], linewidth=1, marker='x', markersize=3)
    
    #utils.plt.title(f'{case_label}')
    #utils.plt.ylim([210, 230])
    #utils.plt.xlim([8, 24])
    utils.plt.xlabel('D (Number of Trips)')
    utils.plt.ylabel('Total Charging Time [minutes]')
    utils.plt.legend()
    utils.plt.grid()
    utils.plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{case}_charging_time_vs_D.{utils.figure_format}')
    utils.plt.savefig(plot_fp, format=utils.figure_format)
    utils.plt.close()

    # Create a plot for the total tour time vs D
    utils.plt.figure(figsize=(5, 3))
    for i, h_value in enumerate(case_df['H'].unique()):
        h_df = case_df[case_df['H'] == h_value]
        # Sum the arrays in halt_times for each row and store the results
        tour_time = 1/60 * (h_df['flight_times'].apply(lambda x: np.sum(x)) +  h_df['halt_times'].apply(lambda x: np.sum(x)))  # Sum each array
        utils.plt.plot(h_df['D'], tour_time, label=f'H={h_value}', color=colors[i % len(colors)], linewidth=1, marker='x', markersize=3)  # Use scatter for clearer visualization
    
    max_tour_time = h_df['T_Max'].iloc[0]/60 # maximum tour time in minutes
    utils.plt.axhline(y=max_tour_time, color='r', linestyle='--', linewidth=1) # Add a horizontal line at maximum possible tour time
    utils.plt.text(x=2, y=max_tour_time-35,  # Position below line
            s=f'{max_tour_time:.0f} mins',  # Text content
            color='r', fontsize=10, ha='center', va='bottom')  # Anchor relative to axis
    
    #utils.plt.title(f'{case_label}')
    #utils.plt.ylim([460, 490])
    #utils.plt.xlim([8, 24])
    utils.plt.xlabel('D (Number of Trips)')
    utils.plt.ylabel('Tour Time [minutes]')
    utils.plt.legend()
    utils.plt.grid()
    utils.plt.tight_layout() # Adjust layout to prevent cutting off labels

    plot_fp = os.path.join(plot_dir, f'{case}_tour_time_vs_D.{utils.figure_format}')
    utils.plt.savefig(plot_fp, format=utils.figure_format)
    utils.plt.close()

    


#### DEBUG: To verify that the hotels generated by K-means is the same and in the same order across all cases. Simultaneously the vertex position data is also verified. ####


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
                    nconss, nvars, optimal_sol, process_time) = utils.load_results(result_fp)
                
                # load the vertex data
                vertex_data_fp = os.path.join(result_folder_path, 'vertex_data.opc')

                # load the hotel position data
                hotel_fp = os.path.join(result_folder_path, f'hotels_H{H}.opc')

                c_pos_hotel_verify, Si_hotel_verify = utils.read_data_file(hotel_fp)
                c_pos_verify, Si_verify = utils.read_data_file(vertex_data_fp)

                c_pos_verify = c_pos_hotel_verify + c_pos_verify
                Si_verify = Si_hotel_verify + Si_verify

                if c_pos_verify == c_pos:
                    print("Hotel and vertex verification succeeded.")
                else:
                    raise RuntimeError("Hotel and vertex verification failed.")






                

