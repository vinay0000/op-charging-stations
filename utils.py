import os, shutil
import pickle
import pandas as pd
import numpy as np
import subprocess

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