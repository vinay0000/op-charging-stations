import subprocess

#### Define nominal parameter values ####
sci_raw_fn = 'Bioassessment_25_scival.csv' # science file name
N = 20
H = 4  # number of hotels
D = 3 # number of trips
T_Max = 28800  # [seconds] = 8 hrs total tour length
T_CH = 1800  # [seconds] = 30 min maximum flight time on full-charge
uav_s = 7 # speed of UAV [m/s]
k_ch = 1 # charging factor
k_dis = 1 # discharge factor
timeout = 2*3600 # timeout in seconds. Enter negative number if optimal value is desired.


def script_call(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout):
    command = [
    'python', 'main_driver.py',
    sci_raw_fn, str(N), str(H), str(D), str(T_Max), str(T_CH), str(uav_s), str(k_ch), str(k_dis), str(timeout)
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

#### Analysis cases ####


D = 5
print("Case 5")
script_call(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, -1)
