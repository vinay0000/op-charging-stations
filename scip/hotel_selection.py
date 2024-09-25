"""
 Divsalar, A., Vansteenwegen, P., & Cattrysse, D. (2013). A variable neighborhood search method for the orienteering problem 
    with hotel selection. International Journal of Production Economics, 145(1), 150-160.
"""

from pyscipopt import Model
import pyscipopt as scip
import matplotlib.pyplot as plt
import numpy as np


##### Predefined constants #####

H = 2 # number of hotels
N = 14 # number of 'nodes'
D = 4 # number of trips

T_d = 5 #  limit to each trip

uav_s = 1 # speed of UAV

##### Read in data file containing the position of the vertices (hotels and nodes) #####
file_name = "hotel_sel"
with open(f"./data/{file_name}") as f:
    lines = f.readlines()

num_of_vertices = len(lines)
print(f"num_of_vertices: {num_of_vertices}")

c_pos = [None] * (H+N)

for i in range(0, H+N):
    x_str, y_str = lines[i].split()
    c_pos[i] = [float(x_str), float(y_str)]
print("c_pos: ", c_pos)

# calculate the time of flight
t = np.empty((H+N, H+N))
print("time of flight")
for i in range(0, H+N):
    for j in range(0, H+N):
        t[i][j] = np.sqrt((c_pos[i][0] - c_pos[j][0])**2 + (c_pos[i][1] - c_pos[j][1])**2)/uav_s
        print(t[i][j]) 

##### Create a new model #####
m = Model()

##### Predefined variables #####

# x_{i,j,d} is a binary variable that indicates whether edge (i,j) is traversed in trip d.
x = [[[m.addVar(vtype="B") for d in range(D)] for j in range(H + N)] for i in range(H + N)]

# subtour elimination variable
u = [m.addVar(vtype="I") for i in range(0, N)]

##### Define the objective function #####

m.setObjective(
    scip.quicksum(x[i][j][d] for i in range(0, H+N) for j in range(0, H+N) for d in range(0, D)),
    sense="maximize"
)

# Add constraints
# The zeroth-vertex (hotel) is the start point (first trip).
m.addCons(scip.quicksum(x[0][l][0] for l in range(1, H+N)) == 1) 

# The first-vertex (hotel) is the end point (last trip).
m.addCons(scip.quicksum(x[k][1][D-1] for k in range(0, H+N)) == 1) 

# Disallow all vertex-to-vertex (including hotels) transitions. (No self-tours.)
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

# limit the vertex-to-hotel and vertex-to-vertex transitions.
# The hotel-to-vertex and hotel-to-hotel transitions are *not* limited. 
for i in range(H, H+N):
    m.addCons(scip.quicksum(x[i][j][d]  for d in range(0, D) for j in range(0, H+N)) <=1) 

# 	Limit each trip to the max time allotted per trip (T_d)
for d in range(0, D):
    m.addCons(scip.quicksum(t[i][j]*x[i][j][d] for i in range(0, H+N) for j in range(0, H+N)) <= T_d)

# subtour elimination
for i in range(H, H+N):
    for j in range(H, H+N):
        m.addCons( u[i-H] - u[j-H] + 1 <= (N-1)* (1 - scip.quicksum(x[i][j][d] for d in range(0, D))) )



# Solve the problem
m.optimize()

def  get_solution():
    # Get the optimal solution
    if m.getStatus() == "optimal":
        print("Optimal solution found:")

        solution = [[[round(m.getVal(x[i][j][d])) for d in range(D)] for j in range(H+N)] for i in range(H+N)]
        
        # Print the values of the decision variables x
        '''
        for i in range(H + N):
            for j in range(H + N):
                for d in range(D):
                    print(solution[i][j][d])
        '''

        # Iterate through the 'd' dimension and print 2D slices
        matrix_3d = np.array(solution)
        for d in range(D):
            for matrix_2d in matrix_3d:
                for row in matrix_2d:
                    print(row[d], end=" ")  # Print the corresponding element in each row
                print() # Add an empty line to separate 2D slices
            print()
        

        # ... (print other solution information as needed) ...

    else:
        raise Exception("Optimal solution not found.")

    return solution


### Print the solution ###
solution = get_solution()


### Plot the solution ###
# Initialize the figure
plt.figure(figsize=(8, 6))

# Separate the first H_count 'H' coordinates from the rest
hotels = c_pos[:H]
nodes = c_pos[H:]

# Extract the 'x' and 'y' values for both 'H' and the rest
x_h, y_h = zip(*hotels)
x_n, y_n = zip(*nodes)

# Create the scatter plot with different marker styles and colors
plt.scatter(x_h, y_h, color='orange', s=100, marker='o', label='Hotel')
plt.scatter(x_n, y_n, color='black', s=100, label='Node')


# Draw lines connecting cities
colors = ['black', 'cyan', 'magenta', 'grey','yellow'] # corresponding to different trips
for i in range(H+N):
    for j in range(H+N):
        for d in range(D):
            if solution[i][j][d] == 1:
                [x1, y1] = c_pos[i]
                [x2, y2] = c_pos[j]
                plt.plot([x1, x2], [y1, y2], color=colors[d], lw=2)

# Add a colorbar
plt.colorbar()



# Show the plot
plt.show()

