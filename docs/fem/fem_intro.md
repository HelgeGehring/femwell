# Quick Introduction to the Finite Element Method
Femwell is a mode solver based on the *Finite Element Method* (FEM). This method is employed to numerically approximate the solution of complex Partial Differential Equaions (PDE), including Maxwell equations and those we can derive from them.

A detailed and formal discussion of the FEM is outside the scope of femwell and is even not necessary to effectively use a mode solver. For this reason here we will just highligth the main concepts and idea behind the method focusing on those who are important for every-day use of a mode solver. (Here suggest some literature).

## Galerkin's methods
The finite element method belong to a broader class of numerical method for the approximate solution of PDEs known as *Galerkin's Methods*. All Galerkin methods start from a weak formulation of Partial Differential Equations. ...

The Galerkin method then consists in assuming that the solution $u$ of the PDE is a linear combination of known functions, called *basis function*, that we indicate with $\phi_i$. In formula

$$ u = \sum_{i = 0}^N u_i \phi_i$$

In this way the final solution is known when the the vector of coefficients $(u_i)$ is known. Since the PDE we deal with are linear, it is possible to show that under this assumption, solving the weak formulation of the PDE is equivalent to solving the linear system
$$A\mathbf{u} = \mathbf{f},$$
where $A$ is known as *stifness matrix* and has elements
$$a_{ij} = a(\phi_i, \phi_j),$$
while $f_i = F(\phi_i)$. 

In this way we have recasted a complex differential problem to the solution of a linear system of equations which can be approached with the main well-developed linear algebra libraries.

However, the Galerkin method remain abstract, in the sense that it does not explain how to choose the basis functions $\phi_i$, even if the accuracy of the final solution ultimately depends upon this choice. Instead the Galerkin methods leaves this task to specific numerical algorithms among which we find the Finite Element Method.

## The Finite Element Method
The finite element method, as we mentioned before, is a Galerkin method, which depend on

- geometry definition
- [meshing](chapter-mesh)
- selection of the element family (basis function)
- assembly
- solution of the associated linear system
- post-processing