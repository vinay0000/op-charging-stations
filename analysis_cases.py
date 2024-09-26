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

def script_call(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout):
    command = [
        'python', 'main_driver.py',
        sci_raw_fn, str(N), str(H), str(D), str(T_Max), str(T_CH),
        str(uav_s), str(k_ch), str(k_dis), str(timeout)
    ]
    print(f"Executing command: {' '.join(command)}")
    run_script(command)

#### Analysis cases ####
D = 3
print("Case D = 3")
script_call(sci_raw_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, -1)
