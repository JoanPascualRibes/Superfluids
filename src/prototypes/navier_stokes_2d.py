import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
from dedalus.extras.plot_tools import plot_bot_2d
import matplotlib.animation as animation
figkw = {'figsize':(6,4), 'dpi':100}

#Basis
coords = d3.CartesianCoordinates('x', 'y')
xbasis = d3.RealFourier(coords['x'], 64, bounds=(-np.pi, np.pi))
ybasis = d3.RealFourier(coords['y'], 64, bounds=(-np.pi, np.pi))
dist = d3.Distributor(coords, dtype= np.float64 )

#Fields
u = dist.VectorField(coords, name='u', bases = (xbasis, ybasis))

# Problem
problem = d3.IVP([u], namespace=locals())

#Substitutions
dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])
lap =lambda A: d3.Laplacian(A)

#Equations
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
nu = 0.5
problem.add_equation('dt(u) -nu*lap(u) = -u@grad(u)')

#Initial value
u['g'][0] = 1
u['g'][1] = x
#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 5

# Main loop

timestep = 0.01
u_list = [np.copy(u['g'])]
t_list = [solver.sim_time]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 1 == 0:
        u_list.append(np.copy(u['g']))
        t_list.append(solver.sim_time)
    if solver.iteration % 1000 == 0:
        print('Completed iteration {}'.format(solver.iteration))

u_array = np.array(u_list)
t_array = np.array(t_list)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
def update(n):
    plt.cla()
    speed = np.sqrt(u_array[n,0]**2 + u_array[n,1]**2)
    lw = 2*speed / speed.max()
    plt.streamplot(X, Y, u_array[n,0], u_array[n,1], density=0.6, color='k' , linewidth=lw)
    #plt.streamplot(X, Y, u_array[n,0], u_array[n,1], color=u_array[n,0], linewidth=2, cmap='autumn')
    plt.title('Varying Line Width')
    plt.ylabel('x')
    plt.xlabel('y')

ani = animation.FuncAnimation(fig, update, blit=False, frames=len(u_array), interval = 2, repeat = False )
plt.show()