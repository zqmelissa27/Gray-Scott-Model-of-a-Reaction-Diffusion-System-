# Gray-Scott-Model-of-a-Reaction-Diffusion-System
Numerical simulation of the Gray–Scott reaction–diffusion model using finite differences, focusing on pattern formation arising from nonlinear chemical kinetics and diffusion processes in 2D. 

The Gray–Scott model describes the evolution of two interacting chemical species u(x,y,t) and v(x,y,t):
$$
\begin{aligned}
\frac{\partial u}{\partial t} &= D_u \nabla^2 u - uv^2 + F(1-u), \\
\frac{\partial v}{\partial t} &= D_v \nabla^2 v + uv^2 - (F + k)v
\end{aligned}
$$

## Numerical Method
The system is solved on a two-dimensional Cartesian grid using the FDM
Second-order central differences are used for the spatial discretization of the Laplacian operator,
and time integration is performed using an explicit time-stepping scheme.
## Initial and Boundary Conditions
The concentrations are initialized with a small perturbation around a homogeneous steady state.
Boundary conditions are taken to be periodic (or specify Neumann / Dirichlet if different),
allowing the emergence of self-organized spatial patterns.


