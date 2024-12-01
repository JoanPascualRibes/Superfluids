import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
from dedalus.extras.plot_tools import plot_bot_2d
import matplotlib.animation as animation

#Basis
coords = d3.CartesianCoordinates('x', 'y')
xbasis = d3.RealFourier(coords['x'], 64, bounds=(-np.pi, np.pi))
ybasis = d3.RealFourier(coords['y'], 64, bounds=(-np.pi, np.pi))
dist = d3.Distributor(coords, dtype= np.float64 )

#Fields
v = dist.VectorField(coords, name='v', bases=(xbasis, ybasis))
vs = dist.VectorField(coords, name='vs', bases=(xbasis, ybasis))
s = dist.Field(name='s', bases=(xbasis, ybasis))
rho = dist.Field(name='rho', bases=(xbasis, ybasis))

#Variables
xis0 = 1
xin0 = 1
xip = 1
s0 = 1
u1 = 1
rho0 = 1
T0 = 1
C = 1

# Problem
problem = d3.IVP([v, vs,s,rho], namespace=locals())

#Substitutions
dx = lambda A: d3.Differentiate(A, coords['x'])
dy = lambda A: d3.Differentiate(A, coords['y'])
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
lap =lambda A: d3.Laplacian(A)
xin = xin0 + xip*(s-s0)
xis = xis0 - xip*(s-s0)
vn = (v - xis*vs)/xin
vns = vn - vs
idn = dist.TensorField( (coords,coords), name='idn', bases=(xbasis,ybasis))
idn['g'][0,0] = 1
idn['g'][0,1] = 0
idn['g'][1,0] = 0
idn['g'][1,1] = 1
P = u1**2*(rho-rho0)
grad_vn = d3.Gradient(vn)
grad_vnt = d3.TransposeComponents(grad_vn)
T = T0*(1+(s-s0)/C)

#Equations

problem.add_equation('dt(rho)  = - v@grad(rho) -rho*div(v) ')
problem.add_equation('dt(s) = -v@grad(s) -(1/rho)*div(rho*s*xis*vns)')
problem.add_equation('dt(v) = -v@grad(v) -(1/rho)*div(rho*xin*xis*vns*vns + idn*P )')
problem.add_equation('dt(vs) = -v@grad(vs) +xin*grad_vnt@vns - grad(P)/rho + s*grad(T)')

#Initial value
rho['g'] = 1.e-2*(np.sin(x) - np.sin(y))


#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 1

# Main loop
timestep = 1.e-3
rho_list = [np.copy(rho['g'])]
s_list = [np.copy(s['g'])]
v_list = [np.copy(v['g'])]
vs_list = [np.copy(vs['g'])]
t_list = [solver.sim_time]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 10 == 0:
        rho_list.append(np.copy(rho['g']))
        s_list.append(np.copy(s['g']))
        v_list.append(np.copy(v['g']))
        vs_list.append(np.copy(vs['g']))
        t_list.append(solver.sim_time)
    if solver.iteration % 100 == 0:
        print('Completed iteration {}'.format(solver.iteration))

rho_array = np.array(rho_list)
s_array = np.array(s_list)
v_array = np.array(v_list)
vs_array = np.array(vs_list)
t_array = np.array(t_list)
X, Y = np.meshgrid(x, y)

#Plot
fig, axs = plt.subplots(1, 4, figsize=(11, 4))
axs[0].set_title('rho')
axs[1].set_title('s')
axs[2].set_title('v')
axs[3].set_title('vs')

def update(n):
    plt.cla()
    axs[0].pcolormesh(rho_array[n],  vmin=-1, vmax=1, cmap='RdBu_r')
    axs[1].pcolormesh(s_array[n],  vmin=-1, vmax=1, cmap='RdBu_r')
    speed = np.sqrt(v_array[n,0]**2 + v_array[n,1]**2)
    speeds = np.sqrt(vs_array[n,0]**2 + vs_array[n,1]**2)
    lw = 2*speed / speed.max()
    lws = 2*speeds / speeds.max()
    axs[2].streamplot(X, Y, v_array[n,0], v_array[n,1], density=0.3, color='k' , linewidth=lw, broken_streamlines=False)
    axs[3].streamplot(X, Y, vs_array[n,0], vs_array[n,1], density=0.3, color='k' , linewidth=lws, broken_streamlines=False)


ani = animation.FuncAnimation(fig, update, blit=False, frames=len(s_array), interval = 2, repeat = False )
plt.show()
#ani.save(filename="/tmp/pillow_example.gif", writer="pillow")