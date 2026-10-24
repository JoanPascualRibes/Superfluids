# Superfluid helium-II simulation

![Simulation snapshot](https://github.com/user-attachments/assets/54fb43a3-50d1-410d-8752-765b1c3a5ef6)

Team project on superfluidity in helium-II. The [paper](paper/superfluids.pdf) introduces
Bose-Einstein condensation and the two-fluid picture of He-II, and then tests some of its
properties with our own 2D simulation of the superfluid, written with spectral methods in
[Dedalus](https://dedalus-project.org/).

## Model

We solve the two-fluid model of Kincl, Schmoranzer and Pavelka (Physics of Fluids 35, 047124, 2023)
for the density $\rho$, specific entropy $s$, coflow velocity $\mathbf{v}$ and superfluid velocity
$\mathbf{v}_s$:

$$
\begin{aligned}
\frac{D\rho}{Dt} &= -\rho\,\nabla\cdot\mathbf{v}, \\
\frac{Ds}{Dt} &= -\frac{1}{\rho}\nabla\cdot(\rho s \chi_s \mathbf{v}_{ns}) + \frac{\beta}{\rho}\Delta T + \frac{\zeta}{\rho T}, \\
\frac{D\mathbf{v}}{Dt} &= -\frac{1}{\rho}\nabla\cdot(\rho\chi_n\chi_s\,\mathbf{v}_{ns}\otimes\mathbf{v}_{ns} + p I) + \frac{2\mu}{\rho}\nabla\cdot D_n, \\
\frac{D\mathbf{v}_s}{Dt} &= \chi_n (\nabla\mathbf{v}_n)^T\mathbf{v}_{ns} - \frac{\nabla p}{\rho} + s\nabla T,
\end{aligned}
$$

with $\mathbf{v}_{ns} = \mathbf{v}_n - \mathbf{v}_s$ and $p$, $T$, $\chi_n$ linearised around
$T_0 = 1.9$ K (section 5 of the paper). The domain is periodic in 2D, discretised with a real
Fourier basis and integrated with a second-order Runge-Kutta scheme (RK222).

## Results

| Simulation | Initial condition | What we see | Video |
|---|---|---|---|
| 1. First and second sound | velocity pulse $v_x = e^{-10x^2}$ | density waves at ~40 m/s (first sound) and entropy waves at ~19 m/s, against the predicted 18.8 m/s for second sound | [mp4](media/videos/simulation_1_first_second_sound.mp4) |
| 2. Peak in density | $\rho = \rho_0 + e^{-10x^2}$ | the density peak splits into two sound waves and creates entropy waves that travel at the speeds of first and second sound | [mp4](media/videos/simulation_2_density_peak.mp4) |
| 3. Peak in entropy | $s = s_0 + e^{-10x^2}$ | two density waves: one follows second sound, one travels ahead at the speed of sound | [mp4](media/videos/simulation_3_entropy_peak.mp4) |
| 4. Standing waves | $\mathbf{v} = (0.1\cos y, -0.1\cos x)$ | density and entropy oscillate at ~8.8 Hz and ~4 Hz, i.e. waves at 39 m/s and 17.7 m/s; the superfluid component stays irrotational | [mp4](media/videos/simulation_4_standing_waves.mp4) |
| 5. Fountain effect | hot spot $s = s_0 + 100\,e^{-(x^2+y^2)}$ | superfluid flows towards the hot spot, raising pressure and density there | [mp4](media/videos/simulation_5_fountain_effect.mp4) |

Snapshots of every run are in the appendix of the paper. Boundaries, superfluid creep and
quantum vortices are outside the scope of this model (see the conclusion of the paper).

## Running it

Dedalus is easiest to install with conda (see the
[installation guide](https://dedalus-project.readthedocs.io/en/latest/pages/installation.html)):

```
conda create -n superfluids -c conda-forge dedalus matplotlib
conda activate superfluids
python src/superfluid_2d.py
```

`src/superfluid_2d.py` runs the full model and animates the density, entropy, pressure,
temperature and velocity fields. The initial conditions are set in the `#Initial value` block;
the commented lines there are the ones we used for other runs. To save frames, uncomment the
`savefig` line in `update`.

## Repository layout

```
paper/superfluids.pdf   paper (theory, derivations, results)
src/superfluid_2d.py    two-fluid He-II solver
src/prototypes/         earlier steps: heat equation, 1D Euler, 2D (compressible)
                        Navier-Stokes, first version of the solver, a C++ attempt
media/videos/           videos of the simulations
media/snapshots/        frames from two runs
```

## Authors

Eloi Estèvez, Eki González, Joan Pascual and Timothy Skipper (January 2025).

Thanks to Ondřej Kincl (University of Trento) for helping us understand his paper.
