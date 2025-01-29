# Superfluid Simulation
![image](https://github.com/user-attachments/assets/54fb43a3-50d1-410d-8752-765b1c3a5ef6)

## Overview
We make an introduction to the superfluidity phenomenon, focusing on the particular case
of helium-II, dealing with the basics of Bose-Einstein condensation for a better understanding.
We comment and justify the unique properties of He-II, by performing our own simulation of
the superfluid in which we observe some of these behaviors and drawing conclusions. 


This repository contains the code, videos, and paper associated with our research on the simulation of superfluid phenomena, particularly in helium-II. 

The project explores Bose-Einstein condensates, the two-fluid model, and key superfluid behaviors, using numerical simulations with spectral methods with Dedalus.

## Paper

The full paper can be found [here](./Superfluids.pdf). It provides a theoretical background, derivations of the governing equations, and an analysis of the simulation results.

## Simulations

The code simulates various superfluid phenomena, including:

1. **First and Second Sound** - Demonstrates the propagation of density and entropy waves.
2. **Peak in Density** - Examines how density perturbations influence entropy.
3. **Peak in Entropy** - Studies the effects of entropy perturbations on density.
4. **Standing Waves** - Shows the oscillatory behavior of density and entropy.
5. **Fountain Effect** - Simulates the thermomechanical effect in superfluid helium.

Snapshots from the simulations are included in the paper's appendix, and the corresponding videos can be found in this repository.

## Code Structure

```
/                     # Root directory
|-- Code/              # Source code for simulations
|-- Simulation X.mp4     # Video recordings of simulation results
|-- Superfluids.pdf   # Research paper
|-- README.md         # This document
```

## Requirements

To run the simulations, you need:

- **Dedalus framework** for spectral simulations (Python-based).
- Python and basic libraries.

### Installation

To install the latest version of Dedalus Project, we refer you to the [official website](https://dedalus-project.org/).


### Results
The generated data can be visualized using Python scripts included in the repository. You can also see the conclusions with the uploaded videos. Feel free to experiment with the initial conditions of the superfluid.

## Contributors

- Eloi Estèvez
- Eki González
- Joan Pascual
- Timothy Skipper

## References

For a detailed list of references, see the bibliography in the [paper](./Superfluids.pdf).

