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
xin0 = 0.4195
xis0 = 1-xin0
xip = 5.697*1.e-4
s0 = 725.5
u1 = 40
rho0 = 145.5
T0 = 1.9
C = 3902

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
P = (u1**2)*(rho-rho0)
grad_vn = d3.Gradient(vn)
grad_vnt = d3.TransposeComponents(grad_vn)
T = T0*(1+(s-s0)/C)

#Equations

problem.add_equation('dt(rho)  = - v@grad(rho) -rho*div(v) ')
problem.add_equation('dt(s) = -v@grad(s) -(1/rho)*div(rho*s*xis*vns)')
problem.add_equation('dt(v) = -v@grad(v) -(1/rho)*div(rho*xin*xis*vns*vns + idn*P )')
problem.add_equation('dt(vs) = -v@grad(vs) +xin*grad_vnt@vns - grad(P)/rho + s*grad(T)')


#Initial value
rho['g'] = 1000 + np.exp(-(y**2)) 
s['g'] = 1000 
v['g'][0] = 0.0001
v['g'][1] = 0
vs['g'][0] = 0.0001
vs['g'][1] = 0

#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 1.e-2

# Main loop
timestep = 1.e-5
rho_list = [np.copy(rho['g'])]
s_list = [np.copy(s['g'])]
v_list = [np.copy(v['g'])]
vs_list = [np.copy(vs['g'])]
t_list = [solver.sim_time]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 25 == 0:
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
fig, axs = plt.subplots(1, 4, figsize=(16, 4))


def update(n):
    axs[0].cla()
    axs[1].cla()
    axs[2].cla()
    axs[3].cla()
    #axs[4].cla()
    axs[0].set_title('rho')
    axs[1].set_title('s')
    axs[2].set_title('v')
    axs[3].set_title('vs')
    #axs[4].set_title('speed')
    axs[0].imshow(np.transpose(rho_array[n]),  cmap='RdBu_r')
    axs[1].pcolormesh(np.transpose(s_array[n]),   cmap='RdBu_r')
    speed = np.sqrt(v_array[n,0]**2 + v_array[n,1]**2)
    speeds = np.sqrt(vs_array[n,0]**2 + vs_array[n,1]**2)
    if np.max(speed) > 0:
        lw = 2*speed / np.max(speed)
    else:
        lw = speed
    if np.max(speeds) > 0:
        lws = 2*speeds / np.max(speeds)
    else:
        lws = speeds
    axs[2].streamplot(X, Y, np.transpose(v_array[n,0]), np.transpose(v_array[n,1]), density=0.5, color='k' , linewidth=lw, broken_streamlines=True)
    axs[3].streamplot(X, Y, np.transpose(vs_array[n,0]), np.transpose(vs_array[n,1]), density=0.5, color='k' , linewidth=lws, broken_streamlines=True)
    #axs[4].pcolormesh(speed, cmap='RdBu_r')

ani = animation.FuncAnimation(fig, update, blit=False, frames=len(rho_list), interval = 200, repeat = True )
plt.show()

#ani.save(filename="Superfluid.gif", writer="pillow")