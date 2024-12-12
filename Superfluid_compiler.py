import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
from dedalus.extras.plot_tools import plot_bot_2d
import matplotlib.animation as animation

#Basis
coords = d3.CartesianCoordinates('x', 'y')
xbasis = d3.RealFourier(coords['x'], 96, bounds=(-np.pi, np.pi))
ybasis = d3.RealFourier(coords['y'], 96, bounds=(-np.pi, np.pi))
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
beta = 1.e1
mu = 1.e1

u2 = np.sqrt( xis0*T0*s0**2 /(xin0*C) )
print(u2)

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
problem.add_equation('dt(tr)  = -25*v@grad(tr)')
problem.add_equation('dt(trs)  = -vs@grad(trs)')


#Initial value
rho['g'] = rho0 + 5*np.exp(-(x-1)**2-y*y) + 5*np.exp(-(x+1)**2-y*y)
s['g'] = s0

# rho['g'] = 1000 + 100*np.exp(-x*x-y*y)
# rho['g'] = 1000 #+ 10*np.sin(x + y)
# rho['g'] = 1000 + 100*np.exp(-4*y*y)
# s['g'] = 1000 #+ 10*np.exp(-x*x-y*y)
# s['g'] = 1000+1*np.exp(-4*x*x-4*y*y)*np.exp(-x*x-y*y)
v['g'][0] = 0.0
v['g'][1] =  0.0
vs['g'][0] = 0.0
vs['g'][1] = 0.0

#tracer
tr['g'] = np.sin(8*x) + np.sin(8*y)
trs['g'] = np.sin(8*x) + np.sin(8*y)

#Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 1.e-1

print(np.shape(v['g']), np.shape(vs['g']), np.shape(vn['g']))
# 85% ; 10:16 @ 1.e-0 , 5.e-5 -->> 71% ; 10:33 completed simulation -->>  11:32 process killed, why???
# Main loop
timestep = 1.e-4
rho_list = [np.copy(rho['g'])]
s_list = [np.copy(s['g'])]
v_list = [np.copy(v['g'])]
vs_list = [np.copy(vs['g'])]
vn_list = [np.copy(vn['g'])]
T_list = [np.copy(T['g'])]
P_list = [np.copy(P['g'])]
tr_list = [np.copy(tr['g'])]
t_list = [solver.sim_time]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 25 == 0:
        rho_list.append(np.copy(rho['g']))
        s_list.append(np.copy(s['g']))
        v_list.append(np.copy(v['g']))
        vs_list.append(np.copy(vs['g']))
        vn_list.append(np.copy(vn['g']))
        T_list.append(np.copy(T['g']))
        P_list.append(np.copy(P['g']))
        tr_list.append(np.copy(tr['g']))
        t_list.append(solver.sim_time)
    if solver.iteration % 100 == 0:
        print('Completed iteration {}'.format(solver.iteration))

rho_array = np.array(rho_list)
s_array = np.array(s_list)
v_array = np.array(v_list)
vs_array = np.array(vs_list)
vn_array = np.array(vn_list)
vns_array = vn_array - vs_array
T_array = np.array(T_list)
P_array = np.array(P_list)
tr_array = np.array(tr_list)
t_array = np.array(t_list)
X, Y = np.meshgrid(x, y)

#Plot
fig, axs = plt.subplots(3, 3, figsize=(16, 13))

rho_min, rho_max = max(np.min(rho_array),0), np.max(rho_array)
s_min, s_max = max(np.min(s_array),0), np.max(s_array)
T_min, T_max = max(np.min(T_array),0), np.max(T_array)
P_min, P_max = np.min(P_array), np.max(P_array)

im_rho = axs[0,0].imshow(np.transpose(rho_array[0])[::-1], vmin=rho_min , vmax=rho_max, cmap='jet')
im_s = axs[1,0].imshow(np.transpose(s_array[0])[::-1],vmin=s_min, vmax=s_max,   cmap='plasma')
im_P = axs[0,1].imshow(np.transpose(P_array[0])[::-1],vmin=P_min, vmax=P_max,   cmap='BuPu')
im_T = axs[1,1].imshow(np.transpose(T_array[0])[::-1],vmin=T_min, vmax=T_max,   cmap='OrRd')
fig.colorbar(im_rho, ax= axs[0,0], fraction=0.046, pad=0.04)
fig.colorbar(im_s, ax= axs[1,0], fraction=0.046, pad=0.04)
fig.colorbar(im_P, ax= axs[0,1], fraction=0.046, pad=0.04)
fig.colorbar(im_T, ax= axs[1,1], fraction=0.046, pad=0.04)

def update(n):
    fig.suptitle(f"t={t_array[n]:.4f}")
    for i in range(3):
        for j  in range(3):
            axs[i,j].cla()

    axs[0,0].set_title('$\\text{Density field: } \\rho $')
    axs[1,0].set_title('$\\text{Entropy field: } s$')
    axs[0,1].set_title('$\\text{Pressure difference field: } P$')
    axs[1,1].set_title('$\\text{Temperature field: } T$')
    axs[0,2].set_title('Tracer')
    axs[1,2].set_title('$\\text{Total fluid velocity field: } v$')
    axs[2,0].set_title('$\\text{Superfluid velocity field: } v_s$')
    axs[2,1].set_title('$\\text{Normal fluid velocity field: } v_n$')
    axs[2,2].set_title('$\\text{Conterflow velocity field: } v_n - v_s$')
    

    axs[0,0].imshow(np.transpose(rho_array[n])[::-1], vmin=rho_min, vmax=rho_max, cmap='jet', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    axs[1,0].imshow(np.transpose(s_array[n])[::-1], vmin=s_min, vmax=s_max, cmap='plasma', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    axs[0,1].imshow(np.transpose(P_array[n])[::-1], vmin=P_min, vmax=P_max, cmap='BuPu', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    axs[1,1].imshow(np.transpose(T_array[n])[::-1], vmin=T_min, vmax=T_max, cmap='OrRd', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    axs[0,2].imshow(np.transpose(tr_array[n])[::-1], cmap='RdBu_r', extent=[np.min(x), np.max(x), np.min(y), np.max(y)])

    speed = np.sqrt(v_array[n,0]**2 + v_array[n,1]**2)
    speeds = np.sqrt(vs_array[n,0]**2 + vs_array[n,1]**2)
    speedn = np.sqrt(vn_array[n,0]**2 + vn_array[n,1]**2)
    speedns = np.sqrt(vns_array[n,0]**2 + vns_array[n,1]**2)

    axs[1,2].streamplot(X, Y, np.transpose(v_array[n,0]), np.transpose(v_array[n,1]), density=0.6, color=speed ,  broken_streamlines=True, cmap ='jet')
    axs[2,0].streamplot(X, Y, np.transpose(vs_array[n,0]), np.transpose(vs_array[n,1]), density=0.6, color=speeds ,  broken_streamlines=True, cmap='YlGnBu')
    axs[2,1].streamplot(X, Y, np.transpose(vn_array[n,0]), np.transpose(vn_array[n,1]), density=0.6, color=speedn ,  broken_streamlines=True, cmap='YlOrBr')
    axs[2,2].streamplot(X, Y, np.transpose(vns_array[n,0]), np.transpose(vns_array[n,1]), density=0.6, color=speedn ,  broken_streamlines=True, cmap='RdPu')

    # fig.savefig(f"./snapshots/SF01_{n}.png") # Descomentar para guardar los frames

ani = animation.FuncAnimation(fig, update, blit=False, frames=len(rho_list), interval = 20, repeat = True, repeat_delay = 100)
plt.show()
# ani.save(filename="Superfluid01.gif", writer="pillow")
