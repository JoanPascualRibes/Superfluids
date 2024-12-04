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

#Tracers
tr = dist.Field(name = 'tr', bases = (xbasis,ybasis))
trs = dist.Field(name = 'trs', bases = (xbasis,ybasis))

#Variables
xin0 = 0.4195
xis0 = 1-xin0
xip = 5.697*1.e-4
s0 = 725.5
u1 = 40
rho0 = 145.5
T0 = 1.9
C = 3902
beta = 1.e+3
mu = 1.e3

# Problem
problem = d3.IVP([v, vs,s,rho,tr,trs], namespace=locals())

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
Dn = 1/2*(d3.Gradient(vn) + d3.TransposeComponents(d3.Gradient(vn)))
zeta = beta*d3.Gradient(T)@d3.Gradient(T) + 2*mu*d3.Trace(Dn@d3.TransposeComponents(Dn))

#Equations

problem.add_equation('dt(rho)  = - v@grad(rho) -rho*div(v) ')
problem.add_equation('dt(s) = -v@grad(s) -(1/rho)*div(rho*s*xis*vns) +(beta/rho)*lap(T) + zeta/(rho*T)')
problem.add_equation('dt(v) = -v@grad(v) -(1/rho)*div(rho*xin*xis*vns*vns + idn*P ) + 2*mu/rho*div(Dn)')
problem.add_equation('dt(vs) = -v@grad(vs) +xin*grad_vnt@vns - grad(P)/rho + s*grad(T)')
problem.add_equation('dt(tr)  = -1*v@grad(tr)')
problem.add_equation('dt(trs)  = -vs@grad(trs)')


#Initial value
rho['g'] = 1000 + 100*np.sin(x + y)
s['g'] = 1000 #0.001*(x**2+y**2)
v['g'][0] = 0.01
v['g'][1] =  0
vs['g'][0] = 0.01
vs['g'][1] = 0

#tracer
tr['g'] = np.sin(5*x) + np.sin(5*y)
trs['g'] = np.sin(5*x) + np.sin(5*y)

#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 1*1.e-1

# Main loop
timestep = 1.e-4
rho_list = [np.copy(rho['g'])]
s_list = [np.copy(s['g'])]
v_list = [np.copy(v['g'])]
vs_list = [np.copy(vs['g'])]
tr_list = [np.copy(tr['g'])]
trs_list = [np.copy(trs['g'])]
t_list = [solver.sim_time]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 25 == 0:
        rho_list.append(np.copy(rho['g']))
        s_list.append(np.copy(s['g']))
        v_list.append(np.copy(v['g']))
        vs_list.append(np.copy(vs['g']))
        tr_list.append(np.copy(tr['g']))
        trs_list.append(np.copy(trs['g']))
        t_list.append(solver.sim_time)
    if solver.iteration % 100 == 0:
        print('Completed iteration {}'.format(solver.iteration))

rho_array = np.array(rho_list)
s_array = np.array(s_list)
v_array = np.array(v_list)
vs_array = np.array(vs_list)
tr_array = np.array(tr_list)
trs_array = np.array(trs_list)
t_array = np.array(t_list)
X, Y = np.meshgrid(x, y)

#Plot
fig, axs = plt.subplots(2, 3, figsize=(12, 8))
rho_min = np.min(rho_array)
rho_max = np.max(rho_array)
s_min = np.min(s_array)
s_max = np.max(s_array)
im_rho = axs[0,0].imshow(np.transpose(rho_array[0]), vmin=rho_min , vmax=rho_max, cmap='jet')
im_s = axs[1,0].imshow(np.transpose(s_array[0]),vmin=s_min, vmax=s_max,   cmap='jet')
fig.colorbar(im_rho, ax= axs[0,0])
fig.colorbar(im_rho, ax= axs[1,0])


def update(n):
    fig.suptitle(f"t={t_array[n]:.4f}")
    for i in range(2):
        for j  in range(3):
            axs[i,j].cla()
    axs[0,0].set_title('$ \\rho $')
    axs[1,0].set_title('$s$')
    axs[0,1].set_title('$v$')
    axs[1,1].set_title('$v_s$')
    axs[0,2].set_title('Tracer')
    axs[1,2].set_title('Superfluid Tracer')
    axs[0,0].imshow(np.transpose(rho_array[n]), vmin=rho_min , vmax=rho_max, cmap='jet')
    axs[1,0].imshow(np.transpose(s_array[n]),vmin=s_min, vmax=s_max,   cmap='jet')
    axs[0,2].imshow(np.transpose(tr_array[n]),   cmap='RdBu_r')
    axs[1,2].imshow(np.transpose(trs_array[n]),   cmap='RdBu_r')
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
    axs[0,1].streamplot(X, Y, np.transpose(v_array[n,0]), np.transpose(v_array[n,1]), density=0.5, color=speed ,  broken_streamlines=True, cmap ='jet')
    axs[1,1].streamplot(X, Y, np.transpose(vs_array[n,0]), np.transpose(vs_array[n,1]), density=0.5, color=speeds ,  broken_streamlines=True, cmap='jet')

ani = animation.FuncAnimation(fig, update, blit=False, frames=len(rho_list), interval = 200, repeat = True )
plt.show()

#ani.save(filename="Superfluid.gif", writer="pillow")