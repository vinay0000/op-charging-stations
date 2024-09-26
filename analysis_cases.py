""" This script can be used to run batches of different scenarios.
A separate folder is created for each science case (with different number of vertices).
The filenames of the resulting pickle objects contain meta infromation about the scenario.
"""

import subprocess
import os

#### Define nominal parameter values ####
science_case = 'Bioassessment_25_scival'
N = 20
H = 4  # number of hotels
D = 3 # number of trips
T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = 2*3600 # timeout in seconds. Enter negative number if optimal value is desired.

current_dir = os.path.dirname(os.path.abspath(__file__))
sci_raw_data_fp = os.path.join(current_dir, f'data/science_data/{science_case}.csv') # Absolute path to the raw science data file

def create_directory_if_not_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

result_folder_path = os.path.join(current_dir, f'results/{science_case}') # Absolute folder path where the results are to be written
create_directory_if_not_exists(result_folder_path)

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

#### Analysis cases ####
D = 3
print("Case D = 3")
script_call(sci_raw_data_fp, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, -1, result_folder_path)
