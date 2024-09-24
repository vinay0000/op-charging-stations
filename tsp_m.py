# https://opensourc.es/blog/mip-tsp/

import pyscipopt as scip
import matplotlib.pyplot as plt
import numpy as np

file_name = "tsp_25_m"
with open(f"./data/{file_name}") as f:
    lines = f.readlines()

N = len(lines)
print(f"N: {N}")

c_pos = [None] * N

for i in range(N):
    x_str, y_str = lines[i].split()
    c_pos[i] = [float(x_str), float(y_str)]
print("c_pos: ", c_pos)

# Create a new model
m = scip.Model()

# Define the number of points N
N = len(c_pos)

# Define the total distance which can be traversed.
D = 1000

# Define the variables
x = [[m.addVar(vtype="B") for _ in range(N)] for _ in range(N)]

# Define the objective function
m.setObjective(
    scip.quicksum(x[i][j] for i in range(N) for j in range(N)),
    "maximize"
)

# Add constraints
# The first-city is the base-city

# force the base-city to always be one of the toured cities.
m.addCons(scip.quicksum(x[0][j] for j in range(0, N)) == 1) 
m.addCons(scip.quicksum(x[i][0] for i in range(0, N)) == 1) 

m.addCons(x[0][0] == 0) # no self-tours to base-city

for i in range(1, N):
    m.addCons(x[i][i] == 0) # no self-tours
    m.addCons(scip.quicksum(x[i][j] for j in range(0, N)) <= 1) # a city is entered only once 

for j in range(1, N):
    m.addCons(scip.quicksum(x[i][j] for i in range(0, N)) <= 1) # a city, is exited only once

for f in range(0,N):
    for t in range(0,N):
        m.addCons(x[f][t] + x[t][f] <= 1)


for j in range(1,N):
    m.addCons(scip.quicksum(x[j][k] for k in range(0, N)) - scip.quicksum(x[i][j] for i in range(0, N)) == 0)


m.addCons(scip.quicksum(x[i][j] * scip.quicksum((c_pos[i][k] - c_pos[j][k]) ** 2 for k in range(2)) for i in range(N) for j in range(N)) <= D)

def print_solution(m):    
    # Get the optimal solution
    solution = [[int(m.getVal(x[i][j])) for j in range(N)] for i in range(N)]

    # Print the solution
    for i in range(N):
        print(solution[i])
    
    print("first solution total distance contraint")
    total_distance = scip.quicksum(solution[i][j] * scip.quicksum((c_pos[i][k] - c_pos[j][k]) ** 2 for k in range(2)) for i in range(N) for j in range(N))
    print("total_distance: ", total_distance)
    
    return solution

# Solve the problem
m.optimize()
#  print solution
print("first solution")
solution = print_solution(m)


def is_one_cycle(m ,x):

    N = len(x)
    x_val = [m.getVal(var) for row in x for var in row]
    print("length of x_val:", len(x_val))

    # Find cycle which encompasses the base-city (the first city)
    cycle_idx = [0] # base-station
    cycle_idx.append(max(range(N), key=lambda i: x_val[cycle_idx[-1]*N + i]))

    print("Logged cycle through the first-city.")

    while True:
        v, idx = max((x_val[cycle_idx[-1]*N + i], i) for i in range(N))
        if idx == cycle_idx[0]:
            break
        else:
            cycle_idx.append(idx)

    print(f"cycle_idx: {cycle_idx}")
    print(f"Length: {len(cycle_idx)}")

    number_of_cities =  sum(x_val[i] for i in range(N*N)) # total number of cities traversed, counting across all the smaller loops
    print("number_of_cities traversed, counting across all the smaller loops : ", number_of_cities)
    if len(cycle_idx) < number_of_cities:

        m.freeTransform()
        m.addCons(scip.quicksum(x[i][j] for i in cycle_idx for j in cycle_idx) <= len(cycle_idx) - 1) # r

        m.optimize()
        return False
    else:
        return True

    
isOneCycle = False
while not isOneCycle:
    isOneCycle = is_one_cycle(m,x)

print("final solution")
solution = print_solution(m)

### Plot the solution ###

# Assuming `matrix` is your N x N matrix
# Initialize the figure
plt.figure(figsize=(8, 6))

x, y = zip(*c_pos)
plt.scatter(x, y, color='black', s=100)  # Adjust 's' for marker size



# Draw lines connecting cities
for i in range(N):
    for j in range(N):
        if solution[i][j] == 1:
            [x1, y1] = c_pos[i]
            [x2, y2] = c_pos[j]
            plt.plot([x1, x2], [y1, y2], color='red', lw=2)

# Add a colorbar
plt.colorbar()



# Show the plot
plt.show()






