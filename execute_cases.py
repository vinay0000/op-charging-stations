""" This script can be used to run batches of different scenarios.
A separate folder is created for each science case (with different number of vertices).
The filenames of the resulting pickle objects contain meta infromation about the scenario.
"""

import subprocess
import os

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
result_dir = os.path.join(file_dir, f'results')

#### Define nominal parameter values ####
science_case_list = [('Bioassessment_24_scival', 24), ('Bioassessment_50_scival', 50), ('Bioassessment_103_scival', 103), 
                        ('Fracking_24_scival',24), ('Fracking_50_scival',50), ('Fracking_103_scival',103),
                        ('RivNetworkContinuity_24_scival',24), ('RivNetworkContinuity_50_scival',50), ('RivNetworkContinuity_103_scival',103)]


#science_case_list = [
#                    ('Plume_30_scival',30)]

science_case_list = [ ('Bioassessment_50_scival', 50), 
                        ('Fracking_50_scival',50), 
                        ('RivNetworkContinuity_50_scival',50), 
                        ('Bioassessment_103_scival', 103), ('Fracking_103_scival',103), ('RivNetworkContinuity_103_scival',103)]


T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = -1 # timeout in seconds. Enter negative number if optimal value is desired.

D_range = [1,3] #range(1,14) # number of trips
H_range = [2,3] #range(2,4) # number of hotels

for sci_idx, science_case in enumerate(science_case_list):

    sci_case_name =  science_case[0]
    N = science_case[1]
    sci_raw_data_fp = os.path.join(data_dir, f'{sci_case_name}.csv') # Absolute path to the raw science data file
    result_folder_path = os.path.join(result_dir, f'{sci_case_name}') # Absolute folder path where the results are to be written

    create_directory_if_not_exists(result_folder_path)


    for D in D_range:
        for H in H_range:
            print(f'Science case: {science_case}')
            print(f'Case: # Hotels H = {H}, # Trips D = {D}')
            script_call(sci_raw_data_fp, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout, result_folder_path)
