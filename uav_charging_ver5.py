"""
 Divsalar, A., Vansteenwegen, P., & Cattrysse, D. (2013). A variable neighborhood search method for the orienteering problem 
    with hotel selection. International Journal of Production Economics, 145(1), 150-160.

    In ver-5, scores are assigned to nodes, and not to arcs. 
    In ver-5, like in ver-4 hotel-to-hotel transitions are given a Score = 0

    Some times sub-tours are given as the solution. However note that while SCIP states it as a optimal solution, it also says it is a in-feasible solution.
"""

from pyscipopt import Model
import pyscipopt as scip
import numpy as np
import os
import pickle
import argparse
import time

def main(input_data_fn, hotel_fn, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout, result_fn):
    """ Example usage: 
            python uav_charging_ver4.py Fracking_25.ophs 22 3 2 28800 1800 5 5 1 Fracking_25_scival_N22_H3_D2_Tmax28800_Tch1800_UAVs5_kch5_kdis1.pkl

            The results are stored as a pickle object in the 'results' directory, with the same filename as the input filename.
    """

    #### User inputs ####
    # input_data_fn = 'Fracking_25.ophs'
    # N = 22 # number of vertices 
    # H = 3 # number of hotels >=2
    # D = 3 # number of trips
    # Tmax = 28800 # [seconds] = 8 hrs Max tour length 
    # T_CH = 1800 # [seconds] = 30 min Maximum flight time on full-charge 
    # uav_s = 5 # speed of UAV [m/s]
    # k_ch = 5
    # k_dis = 1
    # result_fn = 'Fracking_25_scival_N24_H1_D10_Tmax28800_Tch1800.pkl'
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_data_dir = os.path.join(current_dir,'data/science_data')
    input_data_fp = os.path.join(input_data_dir, input_data_fn)
    hotel_fp = os.path.join(input_data_dir, hotel_fn)

    # Formulate the result file path
    results_dir = os.path.join(current_dir, 'results')
    result_fp = os.path.join(results_dir, result_fn)
    print(result_fp)


    # Start timing here
    start_time = time.time()

    def read_data_file(input_data_fp):
        with open(input_data_fp, 'r') as file:
            """
                location and score:  x y Si
                    x	= x coordinate 
                    y	= y coordinate
                    Si	= score

                    the first line is the starting hotel
                    the second line is the ending hotel
                    the next "H" lines are the extra hotels
                    the remaining lines (N) are the vertices

            """
            # Extract location and score data
            location = []
            Si = []
            next(file)  # Skip the first line (header)
            for line in file:
                    values = list(map(float, line.split()))
                    x, y, _Si = values[-3], values[-2], values[-1]  # Assuming the first two values are x and y, and third is Si
                    location.append([x, y])
                    Si.append(_Si)

        return location, Si

    # read in the input data 
    c_pos_hotel, Si_hotel = read_data_file(hotel_fp)
    c_pos, Si = read_data_file(input_data_fp)

    c_pos = c_pos_hotel + c_pos
    Si = Si_hotel + Si

    print(c_pos)
    print(Si)


    # Print the extracted values
    print(f"number of vertices is N: {N}")
    print(f"number of hotels is H: {H}")
    print(f"number for trips is D: {D}")
    print(f"maximum tour time T_Max: {T_Max}")
    print(f"Max flight time on full charge T_CH: {T_CH}")
    #print("Location  Data:", c_pos)
    #print(len(c_pos))

    # calculate the time of flight
    t = np.empty((H+N, H+N))
    print("calculate time of flight")
    for i in range(0, H+N):
        for j in range(0, H+N):
            t[i][j] = np.sqrt((c_pos[i][0] - c_pos[j][0])**2 + (c_pos[i][1] - c_pos[j][1])**2)/uav_s
            #print(t[i][j]) 


    ##### Create a new model #####
    m = Model()

    ##### Below settings are enabled to prevent SCIP from giving infesible solutio #####
    # Disable presolving
    m.setParam('presolving/maxrounds', 0)  # Disables presolving rounds
    # Tighten feasibility tolerance
    # https://stackoverflow.com/questions/22452332/obtain-best-feasible-solution-with-scip
    m.setParam('numerics/feastol', 1e-8)  # Default is 1e-6, use a smaller value for tighter feasibility. Change to 1e-10 if required.
    m.setEmphasis(3)  #  SCIP_PARAMEMPHASIS_FEASIBILITY https://www.scipopt.org/doc/html/type__paramset_8h_source.php 
    #################################################################

    ##### Predefined variables #####
    # x_{i,j,d} is a binary variable that indicates whether edge (i,j) is traversed in trip d.
    x = [[[m.addVar(vtype="B") for d in range(D)] for j in range(H + N)] for i in range(H + N)]


    # Maximum flight time for trip d 
    t_M = [m.addVar(vtype="C", lb=0, ub=T_Max) for d in range(0, D)]

    #Actual flight time for trip d 
    #Has a strong equality constraint with rest of the variables.
    t_F = [m.addVar(vtype="C", lb=0, ub=T_Max) for d in range(0, D)]

    # Halt time at the end of trip d
    t_H = [m.addVar(vtype="C", lb=0, ub=T_Max) for d in range(0, D)]

    # subtour elimination variable
    u = [m.addVar(vtype="I", lb=0) for i in range(0, N)]


    ##### Define the objective function #####
    m.setObjective(
        scip.quicksum(x[i][j][d]*Si[i] for i in range(0, H+N) for j in range(0, H+N) for d in range(0, D)),
        sense="maximize"
    )

    # Add constraints
    # The zeroth-node (hotel) is the start point (first trip).
    m.addCons(scip.quicksum(x[0][l][0] for l in range(1, H+N)) == 1) 

    # The first-node (hotel) is the end point (last trip).
    m.addCons(scip.quicksum(x[k][1][D-1] for k in range(0, H+N)) == 1) 

    # Disallow all same node-to-node (including hotels) transitions. (No self-tours.)
    for d in range(0, D):
        for k in range(0, H+N):
            m.addCons(x[k][k][d] == 0)

    # Ensure that the trip starts at a hotel.
    # Does not restrict hotel-to-hotel transitions. 
    for d in range(0, D):
        m.addCons(scip.quicksum(x[h][l][d] for h in range(0, H) for l in range(0, H+N)) == 1) 

    # ensure that the trip ends in a hotel
    # does not restrict hotel-to-hotel transitions
    for d in range(0, D):
        m.addCons(scip.quicksum(x[k][h][d] for h in range(0, H) for k in range(0, H+N)) == 1) 

    # ensure trip begins from the same hotel where it last ended. (Connectivity between trips.)
    for d in range(0, D-1):
        for h in range(0, H):
            m.addCons(scip.quicksum(x[k][h][d] for k in range(0, H+N)) - scip.quicksum(x[h][l][d+1] for l in range(0, H+N)) == 0) 

    # ensure connectivity within a trip
    for d in range(0, D):
        for k in range(H, H+N):
            m.addCons(scip.quicksum(x[i][k][d] for i in range(0, H+N)) - scip.quicksum(x[k][j][d] for j in range(0, H+N)) == 0) 

    # limit the same vertex-to-hotel and same vertex-to-vertex transitions to atmost 1.
    # The same hotel-to-vertex and same hotel-to-hotel transitions are *not* limited. But self-tours (same hotel-to-hotel transitions) are disallowed by a previous constraint.
    for i in range(H, H+N):
        m.addCons(scip.quicksum(x[i][j][d]  for d in range(0, D) for j in range(0, H+N)) <=1) 

    # new condition to help convergence
    # restrict same hotel-to-vertex transitions to atmost 1.
    for j in range(H, H+N):
        m.addCons(scip.quicksum(x[i][j][d]  for d in range(0, D) for i in range(0, H)) <=1)     


    '''
    # 	Limit each trip to the max time allotted per trip (T_d)
    for d in range(0, D):
        m.addCons(scip.quicksum(t[i][j]*x[i][j][d] for i in range(0, H+N) for j in range(0, H+N)) <= T_d)
    '''
    # time of flight for trip d
    for d in range(0, D):
        m.addCons(t_F[d] == scip.quicksum(t[i][j]*x[i][j][d] for i in range(H+N) for j in range(H+N))) 

    # maximum available flight time for trip 'd'
    m.addCons(t_M[0] == T_CH)
    for d in range(1, D):
        m.addCons(t_M[d] == t_M[d-1] - t_F[d-1]  + k_ch*t_H[d-1])

    # Flight time must be less than the maximum allowed flight time for the respective trip d.
    for d in range(0, D):
        m.addCons(t_F[d] <= t_M[d]) 
    
    # constraint on the charging/halt time at the end of trip 'd'
    for d in range(0, D):
        m.addCons(k_ch*t_H[d] <= T_CH - (t_M[d] - k_dis*t_F[d]))

    # new condition to help convergence
    #M= 10000000
    M = N*(N-1)/2
    for d in range(0, D):
        m.addCons(M*(1 - scip.quicksum(x[i][j][d] for i in range(0, H) for j in range(0, H))) >= scip.quicksum(x[i][j][d] for i in range(H, H+N) for j in range(H, H+N))) 


    # limit max tour time
    m.addCons(scip.quicksum(t_H[d] for d in range(D)) + scip.quicksum(t[i][j]*x[i][j][d] for i in range(0, H+N) for j in range(0, H+N) for d in range(0, D))<=T_Max)


    # subtour elimination
    for i in range(H, H+N):
        for j in range(H, H+N):
            m.addCons( u[i-H] - u[j-H] + 1 <= (N-1)* (1 - scip.quicksum(x[i][j][d] for d in range(0, D))) )

    # Set a time limit (e.g., 60 seconds)
    if timeout > 0:
        m.setParam('limits/time', timeout)

    # Solve the problem
    m.optimize()

    # End timing here
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"execution time: {elapsed_time} seconds")


    # Write the model to an LP file
    m.writeProblem("output_model.lp")

    def get_solution(timeout=None):

        solution = m.getBestSol() # best solution which could be either an optimal solution or the solution within a timeout period
        
        optimal_sol = 0
        # check if solution is optimal
        if m.getStatus() == "optimal":
            print("Optimal solution found:")
            optimal_sol = 1
        elif(timeout < 0):
            raise Exception("Optimal solution is requested but is not found.")
            
        
        if solution is not None:
            print("Best solution:")
            m.printBestSol()

            gap = m.getGap()
            nconss = m.getNConss() # Retrieve number of constraints
            nvars = m.getNVars() # Retrieve number of variables in the problems

            transitions = [[[round(m.getVal(x[i][j][d])) for d in range(D)] for j in range(H+N)] for i in range(H+N)]
            halt_times = [float(m.getVal(t_H[d])) for d in range(D)]
            max_flight_times = [float(m.getVal(t_M[d])) for d in range(D)]
            flight_times = [float(m.getVal(t_F[d])) for d in range(D)]
            subtour_u = [float(m.getVal(u[i])) for i in range(N)]
            
            # Print the values of the decision variables x
            '''
            for i in range(H + N):
                for j in range(H + N):
                    for d in range(D):
                        print(solution[i][j][d])
            '''

            # Iterate through the 'd' dimension and print 2D slices
            matrix_3d = np.array(transitions)
            # Iterate through the 'd' dimension and print 2D slices
            for d in range(D):
                print(f"Slice {d}:")
                for matrix_2d in matrix_3d:
                    print(' '.join(map(str, [row[d] for row in matrix_2d])))
                print()

            optimal_value = m.getObjVal()

        return transitions, halt_times, max_flight_times, flight_times, subtour_u, optimal_sol, optimal_value, gap, nconss, nvars


    #### Print the solution ####
    transitions, halt_times, max_flight_times, flight_times, subtour_u, optimal_sol, optimal_value, gap, nconss, nvars = get_solution(timeout)

    # Save data objects to a file using pickle
    score = Si
    with open(result_fp, 'wb') as file:
        pickle.dump(H, file)
        pickle.dump(N, file)
        pickle.dump(D, file)
        pickle.dump(T_Max, file)
        pickle.dump(T_CH, file)
        pickle.dump(c_pos, file)
        pickle.dump(Si, file)
        pickle.dump(score, file)
        pickle.dump(transitions, file)
        pickle.dump(halt_times, file)
        pickle.dump(max_flight_times, file)
        pickle.dump(flight_times, file)
        pickle.dump(subtour_u, file)
        pickle.dump(optimal_value, file)
        pickle.dump(gap, file)
        pickle.dump(nconss, file)
        pickle.dump(nvars, file)
        pickle.dump(optimal_sol, file)
        pickle.dump(elapsed_time, file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MIP planner script.')

    parser.add_argument('input_data_fn', type=str, help='Input data filename. File must be present in the data/science_data directory with .ophs extension. e.g., Fracking_25_scival.ophs') 
    parser.add_argument('hotel_fn', type=str, help='Hotel filename')
    parser.add_argument('N', type=int, help='Number of vertices')
    parser.add_argument('H', type=int, help='Number of hotels')
    parser.add_argument('D', type=int, help='Number of trips')
    parser.add_argument('T_Max', type=float, help='Total tour length in seconds (e.g., 28800 for 8 hours)')
    parser.add_argument('T_CH', type=float, help='Maximum flight time on full-charge in seconds (e.g., 1800 for 30 minutes)')
    parser.add_argument('uav_s', type=float, help='speed of UAV [m/s]')
    parser.add_argument('k_ch', type=float, help='charging factor')
    parser.add_argument('k_dis', type=float, help='discharge factor')
    parser.add_argument('timeout', type=float, help='Timeout value in seconds. Enter negative # if optimal solution without timeout is desired.')
    parser.add_argument('result_fn', type=str, help='Name of the pickle object in which the results are to be stored.')

    args = parser.parse_args()
    
    main(args.input_data_fn, args.hotel_fn, args.N, args.H, args.D, args.T_Max, args.T_CH, args.uav_s, args.k_ch, args.k_dis, args.timeout, args.result_fn)