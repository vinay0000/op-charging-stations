""" This script is meant to be used for getting results of random hotel placement.
Prior to running this script, please set the ``algo`` and ``random_state`` parameters of the ``make_hotels.py`` script appropriately.
A separate folder is created for each science case (with different number of vertices).
The filenames of the resulting pickle objects contain meta infromation about the scenario.
"""

import subprocess
import os
import pickle
import numpy as np

def create_directory_if_not_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

def run_script(command):
    """
    Executes a subprocess command and streams output in real-time.

    Args:
        command (list): The command to execute as a list of strings.
    """
    try:
        # Run the script and stream output in real-time
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
            for stdout_line in iter(process.stdout.readline, ''):
                print(stdout_line, end='')  # Stream standard output
            for stderr_line in iter(process.stderr.readline, ''):
                print(stderr_line, end='')  # Stream error output

            process.stdout.close()
            process.stderr.close()
            process.wait()

            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command)

    except subprocess.CalledProcessError as e:
        print(f"Error: Command {e.cmd} failed with return code {e.returncode}")

def script_call(sci_raw_data_fp, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout, result_folder_path):
    command = [
        'python', 'main_driver.py',
        sci_raw_data_fp, str(N), str(H), str(D), str(T_Max), str(T_CH),
        str(uav_s), str(k_ch), str(k_dis), str(timeout), str(result_folder_path)
    ]
    print(f"Executing command: {' '.join(command)}")
    run_script(command)

file_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(file_dir, f'data/science_data/') 
result_dir = os.path.join(file_dir, f'results_random_hotel/random_uniform')

#### Define nominal parameter values ####
#science_case_list = [ ('Bioassessment_24_scival',24), ('Fracking_24_scival',24), ('Plume_24_scival',24), ('RivNetworkContinuity_24_scival',24)                        ]
sci_case_name = 'Bioassessment_24_scival'
N = 24
H = 3
D = 4
T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = -1 # timeout in seconds. Enter negative number if optimal value is desired.


sci_raw_data_fp = os.path.join(data_dir, f'{sci_case_name}.csv') # Absolute path to the raw science data file
result_parent_folder_path = os.path.join(result_dir, f'{sci_case_name}')
create_directory_if_not_exists(result_parent_folder_path)

num_of_iterations = 250

print(f'Science case: {sci_case_name}')
print(f'Case: # Hotels H = {H}, # Trips D = {D}')

optimal_value_list = []
for iter_n in range(num_of_iterations):
    print(f'Iteration number {iter_n}')
    
    result_folder_path = os.path.join(result_parent_folder_path, f'{iter_n}') # Absolute folder path where the results are to be written
    create_directory_if_not_exists(result_folder_path)

    script_call(sci_raw_data_fp, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout, result_folder_path)
    try:
        result_fp = os.path.join(result_folder_path, f'MIPver5__N{N}_H{H}_D{D}_Tmax{T_Max}_Tch{T_CH}_UAVsp{uav_s}_kch{k_ch}_kdis{k_dis}.pkl') 
        
    except:
        pass

print(optimal_value_list)


