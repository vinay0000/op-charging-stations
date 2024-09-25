import subprocess
import argparse

def main(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout):
    ################ user Parameters ################
    #sci_raw_fn = 'Bioassessment_50_scival.csv' # science file name
    #N = 44
    #H = 5  # number of hotels
    #D = 3 # number of trips
    #T_Max = 28800  # [seconds] = 8 hrs total tour length
    #T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
    #uav_s = 7 # speed of UAV [m/s]
    #k_ch = 1 # charging factor
    #k_dis = 1 # discharge factor
    #timeout = 3600 # timeout in seconds. Enter negative number if optimal value is desired.
    ver = 5 # MIP formulation version

    # name of the file containing the processed science data
    sci_fn = sci_raw_fn.replace(".csv", ".ophs")
    # name of the file containing the hotel data (in the data/science_data location)
    hotel_fn = f'hotels_H{H}.ophs'
    # define a (planner) result filename based on the above parameters
    result_fn = f'ver{ver}_' + sci_fn + f'_N{int(N)}' + f'_H{int(H)}' + f'_D{int(D)}' + f'_Tmax{int(T_Max)}' + f'_Tch{int(T_CH)}' + f'_UAVsp{int(uav_s)}' + f'_kch{int(k_ch)}' + f'_kdis{int(k_dis)}' + '.pkl'


    ################ Call the data preprocessing script ################
    command = [
        'python', 'process_input_data_ver2.py',
        sci_raw_fn, sci_fn
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    # Call the script and capture output in real-time
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        # Iterate over the output lines
        for stdout_line in iter(process.stdout.readline, ''):
            print(stdout_line, end='')  # Print each line of output as it is received
        for stderr_line in iter(process.stderr.readline, ''):
            print(stderr_line, end='')  # Print each line of error as it is received

        process.stdout.close()
        process.stderr.close()
        process.wait()  # Wait for the subprocess to finish

    ################ Call the make_hotels script ################
    command = [
        'python', 'make_hotels.py',
        sci_fn, str(H), hotel_fn
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    # Call the script and capture output in real-time
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        # Iterate over the output lines
        for stdout_line in iter(process.stdout.readline, ''):
            print(stdout_line, end='')  # Print each line of output as it is received
        for stderr_line in iter(process.stderr.readline, ''):
            print(stderr_line, end='')  # Print each line of error as it is received

        process.stdout.close()
        process.stderr.close()
        process.wait()  # Wait for the subprocess to finish
    ################################################################

    ################ Call the MIP planner script ################
    command = [
        'python', f'scip/uav_charging_ver{ver}.py',
        sci_fn, hotel_fn, str(N), str(H), str(D), str(T_Max), str(T_CH), str(uav_s), str(k_ch), str(k_dis), str(timeout), result_fn
    ]
    # Call the script and capture output in real-time
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        # Iterate over the output lines
        for stdout_line in iter(process.stdout.readline, ''):
            print(stdout_line, end='')  # Print each line of output as it is received
        for stderr_line in iter(process.stderr.readline, ''):
            print(stderr_line, end='')  # Print each line of error as it is received

        process.stdout.close()
        process.stderr.close()
        process.wait()  # Wait for the subprocess to finish
    ################################################################
    
    ################ Display the results ################
    command = [
        'python', 'print_results_ver2.py', '--disable_plot',
        result_fn
    ]
    # Call the script and capture output in real-time
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        # Iterate over the output lines
        for stdout_line in iter(process.stdout.readline, ''):
            print(stdout_line, end='')  # Print each line of output as it is received
        for stderr_line in iter(process.stderr.readline, ''):
            print(stderr_line, end='')  # Print each line of error as it is received

        process.stdout.close()
        process.stderr.close()
        process.wait()  # Wait for the subprocess to finish
    ################################################################
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MIP planner script.')

    parser.add_argument('sci_raw_fn', type=str, help='Input science filename. File must be present in the data/science_data directory with .ophs extension. e.g., Fracking_25_scival.ophs') 
    parser.add_argument('N', type=int, help='Number of vertices')
    parser.add_argument('H', type=int, help='Number of hotels')
    parser.add_argument('D', type=int, help='Number of trips')
    parser.add_argument('T_Max', type=float, help='Total tour length in seconds (e.g., 28800 for 8 hours)')
    parser.add_argument('T_CH', type=float, help='Maximum flight time on full-charge in seconds (e.g., 1800 for 30 minutes)')
    parser.add_argument('uav_s', type=float, help='speed of UAV [m/s]')
    parser.add_argument('k_ch', type=float, help='charging factor')
    parser.add_argument('k_dis', type=float, help='discharge factor')
    parser.add_argument('timeout', type=float, help='Timeout value in seconds. Enter negative # if optimal solution without timeout is desired.')

    args = parser.parse_args()
    
    main(args.sci_raw_fn, args.N, args.H, args.D, args.T_Max, args.T_CH, args.uav_s, args.k_ch, args.k_dis, args.timeout)