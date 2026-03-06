# Quick Introduction to the Finite Element Method
Femwell is a mode solver based on the *Finite Element Method* (FEM). This method is employed to numerically approximate the solution of complex Partial Differential Equaions (PDE), including Maxwell equations and those we can derive from them. Specifically, femwell computes solutions for the Helmoltz equation.

$$ \nabla^2 \mathbf{E}(x, y, z) + k^2 n^2(x, y, z) \mathbf{E}(x, y, z) = 0,$$

where $n$ is the refractive index distribution of the optical waveguide and $\mathbf{E}$ is the electric field.

Treating in the detail the solution of Eq. (1) require advanced mathematical tools from real and functional analysis, which often make the FEM intimidating compared to other numerical methods like finite-difference based methods. Moreover, compared to other PDE commonly treated in FEM books and tutorials, Eq (1) is a complex-valued vectorial equation, whose solution a specific mathematical formulation. However, using a FEM library should not necessarily require a master in mathematics.

In this tutorial we will not then enter into the technical details of the finite element method, but we will just give a high-level introduction. However, differently from other tutorials, we will focus on some specific aspect of the solution of Eq. (1), which are not common outside the field of computational electromagnetics. Those users and readers that would like to get a more formal treatment of the subject can look at \cite{FEM REFERENCES}.

We will take as example a scalar Poisson equation

```{math}
:label: eq-poisson
\nabla^2 u  = \rho 
```

where $u$ is a scalar function and we will highlight how to extend to the vector case.

## Weak formulation of Partial Differential Equations
Instead of trying to solve the Helmoltz equation directly in the form expressed by Eq. [](eq-poisson), we seek for numerical solution of the so called *variational (or weak) formulation* of Eq. [](eq-poisson).

This is obtained by multiplying the PDE we want to solve by a *test function* $v$ and then integrating it over the domain of interest $\Omega$. In the case of the poisson equation we get

```{math}
:label: eq-poisson-int
\int_\Omega v \, \nabla^2 u \, d^3 \mathbf{x} + \int_\Omega \rho \, v \, d^3 \mathbf{x} = 0. 
```

From Eq. [](eq-poisson-int) we can take a step further and use the [Green identities](https://mathworld.wolfram.com/GreensIdentities.html) obtaining

```{math}
:label: eq-poisson-int
- \int_\Omega \nabla v \cdot \nabla u \, d^3 \mathbf{x} + \int_{\partial\Omega} v (\nabla u) \cdot \hat{\mathbf{n}} \, dS  + \int_\Omega \rho \, v \, d^3 \mathbf{x} = 0. 
```

By picking test functions that vanish on the domain boundary $\partial\Omega$, the surface integral in the last equation disappears leaving us with

```{math}
:label: eq-poisson-weak
\int_\Omega \nabla v \cdot \nabla u \, d^3 \mathbf{x} = \int_\Omega \rho \, v \, d^3 \mathbf{x}. 
```

which can also be written in operator form

```{math}
:label: eq-poisson-weak-operator
a(v,   u) = b(v). 
```

where we have introduced the operators

```{math}
:label: eq-poisson-operator-def
a(v,   u) &= \int_\Omega \nabla v \cdot \nabla u \, d^3 \mathbf{x}, \\ 
f(u) &= \int_\Omega \rho \, v \, d^3 \mathbf{x}. 
```

Eq. [](eq-poisson-weak) and [](eq-poisson-weak-operator) are then the starting point for numerical approximations based on Galerking methods among which we find theFEM.

## Galerkin's methods
The finite element method belong to a broader class of numerical method for the approximate solution of PDEs in weak form known as *Galerkin's Methods*.

The Galerkin method then consists in assuming that the solution $u$ and the test functions $v$ are linear combination of known functions, called *basis function*, that we indicate with $\phi_i$. In formula

```{math}
:label: eq-basis-exp
u &= \sum_{i = 0}^N u_i \phi_i \\
v &= \sum_{i = 0}^N v_i \phi_i
```

In this way the final solution is known when the the vector of coefficients $\mathbf{u} = (u_i)$ is known. By substituting Eq. [](eq-basis-exp) into [](eq-poisson-weak) we get

```{math}
:label: eq-poisson-exp
\sum_{i, j} v_i a(\phi_i, \phi_j) u_j = \sum_i f(\phi_i) v_i.
```

Eq. [](eq-poisson-exp) is satisfied when

```{math}
:label: eq-galerkin
\sum_{j} a(\phi_i, \phi_j) u_j = f(\phi_i).
```

which is a linear system of equations whose expression in matrix form is

```{math}
:label: eq-galerkin-matrix
A \mathbf{u} = \mathbf{b}.
```

where we introduced $A = a(\phi_i, \phi_j)$ which is known as *stifness matrix* and $\mathbf{f} = (f(\phi_i))$.

In this way we have recasted a complex differential problem to the solution of a linear system of equations which can be approached with the main well-developed linear algebra libraries.

However, the Galerkin method remain abstract, in the sense that it does not explain how to choose the basis functions $\phi_i$, even if the accuracy of the final solution ultimately depends upon this choice. Instead the Galerkin methods leaves this task to specific numerical algorithms among which we find the Finite Element Method.

## The Finite Element Method

```{figure} element-examples.svg
:name: fig-element-examples
:alt: example of 1D and 2D elements
:align: center

example of elements obtained from the discretization of a 1D domain (a) and 2D domain (b).
```

The finite element method, as we mentioned before, is a Galerkin method which use polynomial basis functions defined on specific sub-domains that we call *elements* and give name to the method. These sub-domain are defined based on a discretization of the original domain in a procedure called *meshing*. In {numref}`fig-element-examples`

- geometry definitio
- [meshing](chapter-mesh)
- selection of the element family (basis function)
- assembly
- solution of the associated linear system
- post-processing