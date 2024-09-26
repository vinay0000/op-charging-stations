import subprocess
import argparse

# Helper function to run subprocess commands and handle outputs
def run_command(command):
    """
    Runs a subprocess command and streams the output in real-time.
    Args:
        command (list): The command to execute as a list of strings.
    """
    try:
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
            # Stream the command output in real-time
            for stdout_line in iter(process.stdout.readline, ''):
                print(stdout_line, end='') # Print each line of output as it is received
            for stderr_line in iter(process.stderr.readline, ''):
                print(stderr_line, end='') # Print each line of error as it is received
            
            # Ensure the subprocess has finished
            process.stdout.close()
            process.stderr.close()
            process.wait()

            # Check if the process exited with an error
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, command)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command {e.cmd} failed with return code {e.returncode}")

def main(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout):
    """
    Main function that orchestrates the MIP planner flow. 
    The input science data file is expected to be in the 'data' folder.
    The results are written as a pickle file in the 'results' folder.
    Args:
        sci_raw_fn (str): Input science filename (e.g., 'Bioassessment_50_scival.csv')
        N (int): Number of vertices.
        H (int): Number of hotels.
        D (int): Number of trips.
        T_Max (float): Maximum tour length in seconds.
        T_CH (float): Maximum flight time on a full charge.
        uav_s (float): UAV speed in m/s.
        k_ch (float): Charging factor.
        k_dis (float): Discharge factor.
        timeout (float): Timeout value in seconds. Enter a negative value for an optimal solution without timeout.
    """
    # MIP formulation version
    ver = 5

    # Generate filenames based on input parameters
    sci_fn = sci_raw_fn.replace(".csv", ".opc")
    hotel_fn = f'hotels_H{H}.opc'
    result_fn = (f'ver{ver}_{sci_fn}_N{int(N)}_H{int(H)}_D{int(D)}_Tmax{int(T_Max)}_'
                 f'Tch{int(T_CH)}_UAVsp{int(uav_s)}_kch{int(k_ch)}_kdis{int(k_dis)}.pkl')

    # Step 1: Run the data preprocessing script
    print(f"Running data preprocessing for {sci_raw_fn}...")
    run_command(['python', 'process_input_data_ver2.py', sci_raw_fn, sci_fn])

    # Step 2: Run the make_hotels script
    print(f"Generating hotel data for {H} hotels...")
    run_command(['python', 'make_hotels.py', sci_fn, str(H), hotel_fn])

    # Step 3: Run the MIP planner script
    print("Running the MIP planner...")
    mip_command = [
        'python', f'scip/uav_charging_ver{ver}.py',
        sci_fn, hotel_fn, str(N), str(H), str(D), str(T_Max), str(T_CH), 
        str(uav_s), str(k_ch), str(k_dis), str(timeout), result_fn
    ]
    run_command(mip_command)

    # Step 4: Display the results
    print(f"Displaying the results from {result_fn}...")
    run_command(['python', 'print_results_ver3.py', '--cartesian_plot', result_fn])

if __name__ == "__main__":
    # Command-line argument parser
    parser = argparse.ArgumentParser(description='MIP planner script.')

    parser.add_argument('sci_raw_fn', type=str, 
                        help='Input CSV science filename in the data/science_data directory.')
    parser.add_argument('N', type=int, help='Number of vertices.')
    parser.add_argument('H', type=int, help='Number of hotels.')
    parser.add_argument('D', type=int, help='Number of trips.')
    parser.add_argument('T_Max', type=float, help='Total tour length in seconds (e.g., 28800 for 8 hours).')
    parser.add_argument('T_CH', type=float, help='Maximum flight time on a full charge in seconds (e.g., 1800 for 30 minutes).')
    parser.add_argument('uav_s', type=float, help='Speed of UAV in meters per second.')
    parser.add_argument('k_ch', type=float, help='Charging factor.')
    parser.add_argument('k_dis', type=float, help='Discharge factor.')
    parser.add_argument('timeout', type=float, help='Timeout in seconds. Enter negative for no timeout (i.e. try to find optimal solution).')

    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.sci_raw_fn, args.N, args.H, args.D, args.T_Max, args.T_CH, args.uav_s, args.k_ch, args.k_dis, args.timeout)