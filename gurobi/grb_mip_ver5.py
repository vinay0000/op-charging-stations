"""
 Divsalar, A., Vansteenwegen, P., & Cattrysse, D. (2013). A variable neighborhood search method for the orienteering problem 
    with hotel selection. International Journal of Production Economics, 145(1), 150-160.

    This file has been developed off the SCIP version `scip/scip_mip_ver5.py`.

"""

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import os
import pickle
import argparse
import time

import importlib.util

import importlib.util
import os

# Define the path to the file one level up
current_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(current_dir, "../utils.py")

# Create a module specification from the file path
spec = importlib.util.spec_from_file_location("utils", file_path)
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)

def main(vertex_data_fp, hotel_fp, N, H, D, T_Max, T_CH, uav_s, k_ch, k_dis, timeout, result_fp):
    """
    MIP planner script for planning UAV flight paths in the presence of multiple charging stations.
    
    Parameters:
    -----------
    vertex_data_fp : str
        Absolute file path to the file containing the vertex data (i.e. list of vertices and their scores).
    hotel_fp : str
        Absolute file path containing the strategic hotel locations.
    N : int
        Number of vertices.
    H : int
        Number of hotels.
    D : int
        Number of trips.
    T_Max : float
        Maximum tour length in seconds (e.g., 28800 for 8 hours operations).
    T_CH : float
        Maximum flight time on full charge in seconds (e.g., 1800 for 30 minutes UAV flight).
    uav_s : float
        Speed of UAV in meters per second.
    k_ch : float
        Charging factor.
    k_dis : float
        Discharge factor.
    timeout : float
        Timeout value in seconds. Enter a negative value if an optimal solution without timeout is desired.
    result_fp : str
        Absolute file path where the results are to be stored as a pickle file.
    
    Example:
    --------
    python uav_charging_ver5.py Fracking_25.opc 22 3 2 28800 1800 5 1 1 results
    """

    # Start timing here
    start_time = time.time()


    # read in the input data 
    c_pos_hotel, Si_hotel = utils.read_data_file(hotel_fp)
    c_pos, Si = utils.read_data_file(vertex_data_fp)

    c_pos = c_pos_hotel + c_pos
    Si = Si_hotel + Si

    print("c_pos: ", c_pos)
    print(f'len(c_pos) is {len(c_pos)}')
    #print(Si)

    # Print the extracted values
    print(f"number of vertices is N: {N}")
    print(f"number of hotels is H: {H}")
    print(f"number for trips is D: {D}")
    print(f"maximum tour time T_Max: {T_Max}")
    print(f"Max flight time on full charge T_CH: {T_CH}")
    #print("Location  Data:", c_pos)
    
    # calculate the time of flight
    t = np.empty((H+N, H+N))
    print("calculate time of flight")
    for i in range(0, H+N):
        for j in range(0, H+N):
            t[i][j] = np.sqrt((c_pos[i][0] - c_pos[j][0])**2 + (c_pos[i][1] - c_pos[j][1])**2)/uav_s
            #print(t[i][j]) 


    ##### Create a new model #####
    m = gp.Model("opc5")

    ##### Predefined variables #####
    # x_{i,j,d} is a binary variable that indicates whether edge (i,j) is traversed in trip d.
    x = [[[m.addVar(vtype=GRB.BINARY) for d in range(D)] for j in range(H + N)] for i in range(H + N)]


    # Maximum flight time for trip d 
    t_M = [m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=T_Max) for d in range(0, D)]

    #Actual flight time for trip d 
    #Has a strong equality constraint with rest of the variables.
    t_F = [m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=T_Max) for d in range(0, D)]

    # Halt time at the end of trip d
    t_H = [m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=T_Max) for d in range(0, D)]

    # subtour elimination variable TODO: Restrict the upper rbound of u
    u = [m.addVar(vtype=GRB.INTEGER, lb=0) for i in range(0, N)]


    ##### Define the objective function #####
    m.setObjective(
        gp.quicksum(x[i][j][d]*Si[i] for i in range(0, H+N) for j in range(0, H+N) for d in range(0, D)),
        GRB.MAXIMIZE
    )

    # Add constraints
    # The zeroth-node (hotel) is the start point (first trip).
    m.addConstr(gp.quicksum(x[0][l][0] for l in range(1, H+N)) == 1) 

    # The first-node (hotel) is the end point (last trip).
    m.addConstr(gp.quicksum(x[k][1][D-1] for k in range(0, H+N)) == 1) 

    # Disallow all same node-to-node (including hotels) transitions. (No self-tours.)
    for d in range(0, D):
        for k in range(0, H+N):
            m.addConstr(x[k][k][d] == 0)

    # Ensure that the trip starts at a hotel.
    # Does not restrict hotel-to-hotel transitions. 
    for d in range(0, D):
        m.addConstr(gp.quicksum(x[h][l][d] for h in range(0, H) for l in range(0, H+N)) == 1) 

    # ensure that the trip ends in a hotel
    # does not restrict hotel-to-hotel transitions
    for d in range(0, D):
        m.addConstr(gp.quicksum(x[k][h][d] for h in range(0, H) for k in range(0, H+N)) == 1) 

    # ensure trip begins from the same hotel where it last ended. (Connectivity between trips.)
    for d in range(0, D-1):
        for h in range(0, H):
            m.addConstr(gp.quicksum(x[k][h][d] for k in range(0, H+N)) - gp.quicksum(x[h][l][d+1] for l in range(0, H+N)) == 0) 

    # ensure connectivity within a trip
    for d in range(0, D):
        for k in range(H, H+N):
            m.addConstr(gp.quicksum(x[i][k][d] for i in range(0, H+N)) - gp.quicksum(x[k][j][d] for j in range(0, H+N)) == 0) 

    # limit the vertex-to-hotel and vertex-to-vertex transitions to atmost 1.
    # The hotel-to-vertex and hotel-to-hotel transitions are *not* limited. But self-tours (same hotel-to-hotel transitions) are disallowed by a previous constraint.
    for i in range(H, H+N):
        m.addConstr(gp.quicksum(x[i][j][d]  for d in range(0, D) for j in range(0, H+N)) <=1) 

    # new condition to help convergence
    # restrict hotel-to-vertex transitions to atmost 1.
    for j in range(H, H+N):
        m.addConstr(gp.quicksum(x[i][j][d]  for d in range(0, D) for i in range(0, H)) <=1)     

    '''
    # 	Limit each trip to the max time allotted per trip (T_d)
    for d in range(0, D):
        m.addConstr(gp.quicksum(t[i][j]*x[i][j][d] for i in range(0, H+N) for j in range(0, H+N)) <= T_d)
    '''
    # time of flight for trip d
    for d in range(0, D):
        m.addConstr(t_F[d] == gp.quicksum(t[i][j]*x[i][j][d] for i in range(H+N) for j in range(H+N))) 

    # maximum available flight time for trip 'd'
    m.addConstr(t_M[0] == T_CH)
    for d in range(1, D):
        m.addConstr(t_M[d] == t_M[d-1] - t_F[d-1]  + k_ch*t_H[d-1])

    # Flight time must be less than the maximum allowed flight time for the respective trip d.
    for d in range(0, D):
        m.addConstr(t_F[d] <= t_M[d]) 
    
    # constraint on the charging/halt time at the end of trip 'd'
    for d in range(0, D):
        m.addConstr(k_ch*t_H[d] <= T_CH - (t_M[d] - k_dis*t_F[d]))

    # new condition to help convergence
    # if the trip is a hotel-to-hotel transition, there are no vertex-to-vertex transitions
    #M= 10000000
    M = N*(N-1)/2
    for d in range(0, D):
        m.addConstr(M*(1 - gp.quicksum(x[i][j][d] for i in range(0, H) for j in range(0, H))) >= gp.quicksum(x[i][j][d] for i in range(H, H+N) for j in range(H, H+N))) 


    # limit max tour time
    m.addConstr(gp.quicksum(t_H[d] for d in range(D)) + gp.quicksum(t[i][j]*x[i][j][d] for i in range(0, H+N) for j in range(0, H+N) for d in range(0, D))<=T_Max)


    # subtour elimination
    for i in range(H, H+N):
        for j in range(H, H+N):
            m.addConstr( u[i-H] - u[j-H] + 1 <= (N-1)* (1 - gp.quicksum(x[i][j][d] for d in range(0, D))) )

    # Set a time limit (e.g., 60 seconds)
    if timeout > 0:
        m.setParam(GRB.Param.TimeLimit, timeout)

    # Solve the problem
    m.optimize()

    # End timing here
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"execution time: {elapsed_time} seconds")

    # Write the model to an LP file
    #m.write("output_model.lp")

    def get_solution(timeout=None):

        # Status checking
        solver_status = m.Status

        # Check for infeasibility, unboundedness, or infeasibility/unboundedness ambiguity
        if solver_status in (GRB.INF_OR_UNBD, GRB.INFEASIBLE, GRB.UNBOUNDED):
            raise Exception("The model cannot be solved because it is infeasible or unbounded.")

        # Initialize optimal solution flag
        optimal_sol = 0

        # Check for specific statuses
        if solver_status == GRB.OPTIMAL:
            optimal_sol = 1
        elif solver_status == GRB.TIME_LIMIT:
            print('Reached time limit, returning best feasible solution.')
        elif solver_status == GRB.ITERATION_LIMIT:
            print('Iteration limit reached, returning current best solution.')
        elif solver_status == GRB.NODE_LIMIT:
            print('Node limit reached, returning current best solution.')
        elif solver_status == GRB.SUBOPTIMAL:
            print('Suboptimal solution found due to limit or other interruption.')
        elif solver_status == GRB.NUMERIC:
            raise Exception("Solver stopped due to numerical issues.")
        elif solver_status == GRB.INTERRUPTED:
            raise Exception("Optimization was interrupted by the user.")
        else:
            # Catch-all for unexpected statuses
            print(f"Optimization was stopped with unknown status {solver_status}")
            raise RuntimeError("Unknown solver status encountered.")

        # Additional handling if optimal solution is not found and timeout is used
        if optimal_sol == 0 and timeout < 0:
            raise Exception("Optimal solution is requested but not found.")
        
        m.printQuality()
        m.printStats()

        gap = m.getAttr("MIPGap")

        # Retrieve the number of linear constraints
        nconss = m.getAttr("NumConstrs")
        print(f"Number of linear constraints: {nconss}")

        # Retreive miscellaneous constraints (quadratic, SOS, general) and check that they are '0'
        nmisc_conss = m.getAttr("NumQConstrs") + m.getAttr("NumSOS") + m.getAttr("NumGenConstrs")
        
        if nmisc_conss != 0:
            raise RuntimeError("Unexpected constraints!")
        
        # Retrieve the total number of variables
        nvars = m.getAttr("NumVars")
        print(f"Number of variables: {nvars}")

        transitions = [[[round(x[i][j][d].X) for d in range(D)] for j in range(H+N)] for i in range(H+N)]
        halt_times = [float(t_H[d].X) for d in range(D)]
        max_flight_times = [float(t_M[d].X) for d in range(D)]
        flight_times = [float(t_F[d].X) for d in range(D)]
        subtour_u = [float(u[i].X) for i in range(N)]
        
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

        optimal_value = m.getAttr("ObjVal")

        return transitions, halt_times, max_flight_times, flight_times, subtour_u, optimal_sol, optimal_value, gap, nconss, nvars


    # Get the solution
    transitions, halt_times, max_flight_times, flight_times, subtour_u, optimal_sol, optimal_value, gap, nconss, nvars = get_solution(timeout)

    # Store results in a pickle file
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
    parser = argparse.ArgumentParser(description='MIP planner script for planning UAV flight paths in the presence of multiple charging stations.')

    parser.add_argument('vertex_data_fp', type=str, help='Absolute file path to the file containing the vertex data (i.e. list of vertices and their scores).') 
    parser.add_argument('hotel_fp', type=str, help='Absolute file path containing the strategic hotel locations.')
    parser.add_argument('N', type=int, help='Number of vertices')
    parser.add_argument('H', type=int, help='Number of hotels')
    parser.add_argument('D', type=int, help='Number of trips')
    parser.add_argument('T_Max', type=float, help='Total tour length in seconds (e.g., 28800 for 8 hours operations)')
    parser.add_argument('T_CH', type=float, help='Maximum flight time on full-charge in seconds (e.g., 1800 for 30 minutes UAV flight)')
    parser.add_argument('uav_s', type=float, help='speed of UAV [m/s]')
    parser.add_argument('k_ch', type=float, help='charging factor')
    parser.add_argument('k_dis', type=float, help='discharge factor')
    parser.add_argument('timeout', type=float, help='Timeout value in seconds. Enter negative # if no timeout is desired (i.e. try to find optimal solution).')
    parser.add_argument('result_fp', type=str, help='Absolute file path where the results are to be stored as a pickle file.')

    args = parser.parse_args()
    
    main(args.vertex_data_fp, args.hotel_fp, args.N, args.H, args.D, args.T_Max, args.T_CH, args.uav_s, args.k_ch, args.k_dis, args.timeout, args.result_fp)