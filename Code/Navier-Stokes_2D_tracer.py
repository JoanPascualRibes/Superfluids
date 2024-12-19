import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import matplotlib.animation as animation

figkw = {"figsize": (6, 4), "dpi": 100}

# Parameters
nu = 0.5
D = 0.01

# Basis
coords = d3.CartesianCoordinates("x", "y")
xbasis = d3.RealFourier(coords["x"], 64, bounds=(-np.pi, np.pi))
ybasis = d3.RealFourier(coords["y"], 64, bounds=(-np.pi, np.pi))
dist = d3.Distributor(coords, dtype=np.float64)

# Fields
u = dist.VectorField(coords, name="u", bases=(xbasis, ybasis))
# Tracer: sigue el fluido para poder visualizar su movimiento
s = dist.Field(name="s", bases=(xbasis, ybasis))

# Problem
problem = d3.IVP([u, s], namespace=locals())

# Substitutions
dx = lambda A: d3.Differentiate(A, coords["x"])
dy = lambda A: d3.Differentiate(A, coords["y"])
lap = lambda A: d3.Laplacian(A)

# Equations
x = dist.local_grid(xbasis)
y = dist.local_grid(ybasis)
problem.add_equation("dt(u) -nu*lap(u) = -u@grad(u)")
problem.add_equation("dt(s) -D*lap(s)= - u@grad(s)")

# Initial value
# u["g"][0] = 1
u["g"][0] = (1 / 2 + 1 / 2 * (np.tanh((y - 0.5) / 0.1) - np.tanh((y + 0.5) / 0.1))) * 5
u["g"][1] += 1 * np.sin(2 * np.pi * x) * np.exp(-(y**2) / 0.01)
u["g"][1] += 1 * np.sin(2 * np.pi * x) * np.exp(-(y**2) / 0.01)
# u["g"][1] = np.sin(x) + np.sin(y)
s["g"] = np.sin(4 * x) + np.sin(4 * y)


# Solver
solver = problem.build_solver(d3.RK222)
solver.stop_sim_time = 5

# Main loop
vor = -d3.div(d3.skew(u))
u_list = [np.copy(u["g"])]
s_list = [np.copy(s["g"])]
vor_list = [np.copy(vor["g"])]
t_list = [solver.sim_time]

timestep = 0.01
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 1 == 0:
        u_list.append(np.copy(u["g"]))
        s_list.append(np.copy(s["g"]))
        vor_list.append(np.copy(vor["g"]))
        t_list.append(solver.sim_time)
    if solver.iteration % 100 == 0:
        print("Completed iteration {}".format(solver.iteration))

u_list = np.array(u_list)
vel_list = np.sqrt(u_list[:, 0] ** 2 + u_list[:, 1] ** 2)

# Plot
import matplotlib.colors as colors

fig, axs = plt.subplots(1, 3, figsize=(11, 4))
plot_arrays = [s_list, vel_list, vor_list]
axs[0].set_title("Tracer")
axs[1].set_title("Velocity")
axs[2].set_title("Vorticity")

cmap = "plasma"
norm = colors.SymLogNorm(
    linthresh=0.015
)  # La viscosidad decrece rápidamente así que con una escala logarítmica se ve mejor

ims = [
    axs[0].imshow(plot_arrays[0][0], cmap=cmap),
    axs[1].imshow(plot_arrays[1][0], cmap=cmap),
    axs[2].imshow(plot_arrays[2][0], cmap=cmap, norm=norm),
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
