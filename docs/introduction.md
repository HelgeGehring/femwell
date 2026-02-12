# Introduction

Femwell is an open-source mode solver based on the *Finite Element Method* (FEM). Its main focus is finding the eigenmodes of optical waveguides by solving the eigenvalue problem associated to Maxwell equations on a 2D cross section of the waveguide exploiting its translational symmetry.

Femwell is mainly programmed in Python and compared to other similar project based on C++ finite element libraries, its accent is on user-friendliness and transparency to the user. So differently from commercial mode solvers the code is open sourced and available for inspection and customization, but it is also simple enough that you can understand what the code is doing even if you're not an expert in the finite element method.

## Installation

Femwell can be easily installed through pip with the command

```pip install femwell```

In this case the mode solver will use [scipy](https://docs.scipy.org/doc/scipy/index.html) eigensolver to estimate the modes of the waveguide. It is also possible to use the [`slepc`](https://slepc.upv.es/slepc4py-current/docs/) eigensolver. In that case you should also run the following conda command

```conda install slepc4py=*=complex*```

## Quickstart

Once you have installed femwell, the fastest way to work with it is explore the example library, choose the one which most fit your task and adapt the geometry of the waveguide to the simulation you want to run. Some useful page to start with are reported below

- [Basic functionalities: *optical waveguides mode*](sec-waveguide-modes)
- [Electromagnetic theory of optical waveguides]()
- [Introduction to the Finite Element Method]()