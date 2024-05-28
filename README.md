# Femwell

![logo](https://raw.githubusercontent.com/HelgeGehring/femwell/main/logo_inline.svg)

[![Docs](https://github.com/helgegehring/femwell/actions/workflows/docs.yml/badge.svg)](https://HelgeGehring.github.io/femwell/)
[![Build](https://github.com/helgegehring/femwell/actions/workflows/build.yml/badge.svg)](https://github.com/HelgeGehring/femwell/actions/workflows/build.yml)
[![PiPy](https://img.shields.io/pypi/v/femwell)](https://pypi.org/project/femwell/)
[![Downloads](https://static.pepy.tech/badge/femwell/month)](https://pepy.tech/project/femwell)

## Welcome to FEMWELL!

FEMWELL is a physics simulation tool that utilises the Finite Element Method (FEM). With FEMWELL, you can simulate integrated circuits, electronic and photonic systems, and so much more. 
The project is created to provide an Open-Source FEM solver. You are welcome to contribute to FEMWELL or just use it. Any feedback and input are valuable and will make FEMWELL better!

## What is a FEM simulation?

FEM is a method to solve various types of differential equations that appear in physics, maths or engineering. FEM methods are used when describing how fields behave in an inhomogeneous setting, for example are electromagnetic fields in structured matter described by Maxwell’s equations. FEM uses a mesh to discretize the space and solve Maxwell’s equations via these “finite elements”.

A FEM simulation typically involves the following steps:

1.	Discretization (or Meshing) of the structure
2.	Calculation of the problem in each element
3.	Assembly accounting for boundary conditions between the elements
4.	Solution via iterative solvers or direct solution methods
   
FEM can be applied to various problems, from Maxwell’s equations to fluid simulation, heat transport or structural analysis.

## What is special about FEMWELL?

First and foremost: FEMWELL is open source! You can just use FEMWELL, you can contribute to FEMWELL and you can modify FEMWELL to fit your specific problem. 

At the moment we focus on photonic and electronic problems, meaning we concentrate on solving Maxwell’s equation. This is useful to understand the physics in modern devices used in classical or quantum computing technologies. 

We are actively working on extending FEMWELL to address other questions. You can find a list of examples below. 

## How can I use FEMWELL?

The simplest thing it to try out the examples in the browser! Hover the rocket at the top on the example pages and click live code. (Might take some time to load).
For more involved calculations, we recommend installing FEMWELL following the instructions.
If you can to improve FEMWELL, please get in touch with us. We are looking forward to your contributions! 

### Please note:
The documentation is lagging behind the state of code, so there's several features for which there are only examples in the code.

## Features

- Photonic eigenmode solver
- Periodic photonic eigenmode solver
- Electric eigenmode solver
- Thermal mode solver (static and transient)
- Coulomb solver

## Possible Simulations

### Photonic problems
  
- Eigenmodes of waveguides and determining their effective refractive index
- Coupling between neighboring waveguides
- Eigenmodes of bent waveguides
- Propagation loss of circular bends and mode mismatch loss with straight waveguides
- Calculation of the group velocity and its dispersion
- Calculation of overlap-integrals and confinement-factors
- Bragg grating cells
- Grating coupler cells
- Overlap integrals between waveguide modes
- Overlap integral between a waveguide mode and a fiber mode
- Coupled mode theory - coupling between adjacent waveguides
- Heat based photonic phase shifters
- Pockels based photonic phase shifters
- PN junction depletion modulator (analytical)

### Heat transport 
- Static thermal profiles
- Transient thermal behavior

### Electronics problems

- Coplanar waveguide RF design
- Eigenmode of a coaxial cable and its specific impedance
- Eigenmodes of electric transmission lines
  and determining their propagation constant (in work)
- Static electric fields

Something missing? Feel free to open an [issue](https://github.com/HelgeGehring/femwell/issues) :)

## Contributors

- Helge Gehring (Google, WWU Münster)
- Simon Bilodeau (Google, Princeton University)
- Joaquin Matres (Google)
- Marc de Cea Falco (Google, Massachusetts Institute of Technology)
- Lodovico Rossi (Princeton University)
- Doris Reiter (TU Dortmund University)
- Yannick Augenstein (Google, Karlsruhe Institute of Technology)
- Niko Savola (Google, Aalto University)
- Rouven Glauert (Idalab)
- Markus DeMartini (Google)
- Lucas Grosjean (Google, Femto-ST Institute)
- Eliza Leung (University of Adelaide)
- Duarte Silva (Eindhoven University of Technology)

Happy about every form of contribution -
pull requests, feature requests, issues, questions, ... :)
