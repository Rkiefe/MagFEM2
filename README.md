# MagFEM2
This is a simple and efficient implementation of the finite element method in 2D for magnetostatics. I decided to create this for anyone who looks at FEMCE and wants to make somthing similar. This code uses a different approach and is only for 2D and for the linear case. Can't give all my secrets away.

## How it works:

$-\nabla \cdot \vec{H} = \nabla \cdot \vec{M}$

$\nabla^2 u = \nabla \cdot \vec{M}$

Now by using the basic principle behind FEM:

$int_\Omega (1+\chi)\nabla u \cdot \nabla \phi dV = \int_\Omega \chi \vec{H_{ext}}\phi dV$
