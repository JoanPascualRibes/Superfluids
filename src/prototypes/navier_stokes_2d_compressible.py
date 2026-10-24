# Compressible Navier-Stokes where p = rho ** gamma
# With tracer
# Pero solo me funciona si la densidad es constante
# Si lo pongo variable (linea 62) peta despues de unos iteraciones

import numpy as np
import dedalus.public as d3
import logging
import matplotlib.pyplot as plt
from matplotlib import animation

logger = logging.getLogger(__name__)


## Parameters
Lx, Lz = 2, 2
Nx, Nz = 64, 64
eta = 1 / 5e4  # Viscosity
D = 1e-2
gamma = 1.4
# Reynolds = 5e4
# Schmidt = 1
dealias = 3 / 2
stop_sim_time = 2
timestepper = d3.RK222
timestep = 1e-2
dtype = np.float64

## Bases
coords = d3.CartesianCoordinates("x", "z")
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords["x"], size=Nx, bounds=(0, Lx), dealias=dealias)
zbasis = d3.RealFourier(coords["z"], size=Nz, bounds=(-Lz / 2, Lz / 2), dealias=dealias)

## Fields
rho = dist.Field(name="rho", bases=(xbasis, zbasis))
u = dist.VectorField(coords, name="u", bases=(xbasis, zbasis))
s = dist.Field(
    name="s", bases=(xbasis, zbasis)
)  # Tracer: sigue el fluido para poder visualizar su movimiento
vor = -d3.div(d3.skew(u))  # Vorticity

## Substitutions
x, z = dist.local_grids(xbasis, zbasis)
ex, ez = coords.unit_vector_fields(dist)

## Equations
problem = d3.IVP([u, rho, s], namespace=locals())
problem.add_equation("dt(u) =  -u@grad(u) - grad(rho**gamma)/rho + eta/rho*lap(u)")
problem.add_equation("dt(rho) = -div(rho*u)")  # Continuity
problem.add_equation("dt(s) - D * lap(s)= - u @ grad(s)")

## Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

## Initial conditions
u["g"][0] = 1 / 2 + 1 / 2 * (np.tanh((z - 0.5) / 0.1) - np.tanh((z + 0.5) / 0.1))
# u["g"][1] += 0.1 * np.sin(2 * np.pi * x / Lx) * np.exp(-((z - 0.5) ** 2) / 0.01)
# u["g"][1] += 0.1 * np.sin(2 * np.pi * x / Lx) * np.exp(-((z + 0.5) ** 2) / 0.01)
rho["g"] = 1
# rho["g"] = 1 + 0.001 * np.cos(4 * np.pi * x)
s["g"] = np.cos(4 * np.pi * x)

## Main loop
s_list = []
rho_list = []
vor_list = []
t_list = []
try:
    logger.info("Starting main loop")
    while solver.proceed:
        solver.step(timestep)
        if (solver.iteration) % 1 == 0:
            s_list.append(s["g"].copy())
            rho_list.append(rho["g"].copy())
            vor_list.append(vor["g"].copy())
            t_list.append(solver.sim_time)
        if (solver.iteration - 1) % 10 == 0:
            logger.info(
                "Iteration=%i, Time=%e, dt=%e"
                % (solver.iteration, solver.sim_time, timestep)
            )
except:
    logger.error("Exception raised, triggering end of main loop.")
    raise
finally:
    solver.log_stats()


## Plot
fig, axs = plt.subplots(1, 3, figsize=(11, 4))
plot_arrays = [s_list, rho_list, vor_list]
axs[0].set_title("Tracer")
axs[1].set_title("Density")
axs[2].set_title("Vorticity")

cmap = "plasma"

ims = [
    axs[0].imshow(plot_arrays[0][0], cmap=cmap),
    axs[1].imshow(plot_arrays[1][0], cmap=cmap),
    axs[2].imshow(plot_arrays[2][0], cmap=cmap),
]
for i in range(len(plot_arrays)):
    fig.colorbar(ims[i], ax=axs[i], location="bottom")
    # axs[i].set_ylabel("x")
    # axs[i].set_xlabel("y")


def update(n):
    fig.suptitle(f"t={t_list[n]:.2f}")
    for i in range(len(plot_arrays)):
        ims[i].set_data(plot_arrays[i][n])


ani = animation.FuncAnimation(
    fig, update, blit=False, frames=len(s_list), interval=2, repeat=False
)
plt.show()
