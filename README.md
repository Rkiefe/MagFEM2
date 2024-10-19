# MagFEM2
This is a simple and efficient implementation of the finite element method in 2D for magnetostatics. I decided to create this for anyone who looks at FEMCE and wants to make something similar. 

This version was implemented for 2D and considering a linear magnetization-magnetic field response. The user inputs a source magnetic field ("Hext"), and this implementation will calculate the interaction of this source field with a paramagnet. The input magnetic field Hext does not need to be uniforme, it can for example be the stray field of a magnet.

## How it works:

$-\nabla \cdot \vec{H} = \nabla \cdot \vec{M}$ 

Using

$\vec{H} = \vec{H}_{ext} - \nabla u$ and $\vec{M} = \chi \vec{H}$ gets us $\nabla^2 u = \nabla \cdot \vec{M}$

Now by using the basic principle behind FEM: $u=\sum_i u_i \phi_i$:

$\int_\Omega (1+\chi)\nabla u \cdot \nabla \phi dV = \int_\Omega \chi \vec{H}_{ext} \cdot \nabla \phi dV$

$\Omega$ is the magnetic domain. Unfortunately, the surrouding region (air) also has to be solved in because of the open boundary problem. You can implement BEM to deal with this issue and only solve for the internal (magnetic) region.
