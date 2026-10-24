import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
from dedalus.extras.plot_tools import plot_bot_2d
figkw = {'figsize':(6,4), 'dpi':100}

# %%%
from matplotlib import animation

# Bases
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=np.float64)
xbasis = d3.RealFourier(xcoord, 1024, bounds=(0, 300), dealias=2)

# Fields
u = dist.Field(name='u', bases=xbasis)

# Problem
problem = d3.IVP([u], namespace=locals())

# Substitutions
dx = lambda A: d3.Differentiate(A, xcoord)
c = 1
nu = 0.02

# Add main equation, with linear terms on the LHS and nonlinear terms on the RHS
problem.add_equation("dt(u) - nu*dx(dx(u))= - u*dx(u)")

# Build solver
solver = problem.build_solver(d3.RK222)
# Stopping criteria
solver.stop_sim_time = 150
# Setup a sine wave
x = dist.local_grid(xbasis)
u['g'] =  np.exp(- (x-100)**2 /100)


# Setup storage
u.change_scales(1)
u_list = [np.copy(u['g'])]
t_list = [solver.sim_time]

# Main loop
timestep = 0.05
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 10 == 0:
        u.change_scales(1)
        u_list.append(np.copy(u['g']))
        t_list.append(solver.sim_time)
    if solver.iteration % 1000 == 0:
        print('Completed iteration {}'.format(solver.iteration))

# Convert storage lists to arrays
u_array = np.array(u_list)
t_array = np.array(t_list)

# Plot solution
# %%%
fig = plt.figure(figsize=(12, 4), dpi=100)
def update(n):
    if n%1==0:
        plt.cla()
        plt.plot(x, u_array[n])
        plt.xlabel('x')
        plt.ylabel('u')
        plt.title(f'1-D Convection evolution (t={n*timestep:.2f}): u')

fig.tight_layout()
ani = animation.FuncAnimation(
    fig, update, blit=False, frames=solver.stop_sim_time*2, interval=30, repeat=False
)
plt.show()
# %%%

plt.figure(figsize=(8, 9), dpi=100)
plt.pcolormesh(x, t_array, np.abs(u_array), shading='nearest')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('1-D Convection: u')
plt.tight_layout()
plt.show()
