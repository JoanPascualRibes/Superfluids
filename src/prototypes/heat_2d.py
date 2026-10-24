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
u = dist.Field(name='u', bases = (xbasis, ybasis))


# Problem
problem = d3.IVP([u], namespace=locals())

#Substitutions
dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])
lap =lambda A: d3.Laplacian(A)

#Equations
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
problem.add_equation('dt(u)  -lap(u) = 0')

#Initial value
u['g'] = np.sin(x) -np.sin(y)

#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 5

# Main loop
fig, ax = plt.subplots()
timestep = 0.01
artists = []
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 1 == 0:
        u.change_scales(1)
        container = ax.pcolormesh(np.transpose(np.copy(u['g'])), vmin=-1, vmax=1, cmap='RdBu_r')
        artists.append([container])
    if solver.iteration % 1000 == 0:
        print('Completed iteration {}'.format(solver.iteration))
ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=20, repeat=False)
fig.colorbar(container)
plt.show()