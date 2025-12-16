# Solving the incompressible 2d Navier Stokes Equations using the projection method algorithm. 

***For more details on this project and the solution methodology please read the pdf of the full project report.***

This project presents code to solve the 2D transient incompressible Navier Stokes Equations, using the Projection Method algorithm and Finite Difference techniques in MATLAB. The set of equations are solved in a square domain of side length 𝒍 = 𝟏, where the boundary conditions at each of the edges of the domain are given. The equations are solved until steady state is reached, and we are given initial conditions for the velocity components in the x and y directions. CDS is used to discretize the convection term, an Euler explicit time scheme, and all of the variables for pressure and velocity are stored in a collocated grid.

The setup presented in the problem is commonly referred to as the cavity lid driven flow. Usually, this lid is in the top of the square domain; however, for this project, it is located on the West boundary of the square, with a y velocity component of 1.

## Methods Used
This project uses the Projection Method Algorithm on a collocated grid. For each of the terms in the equation: 
- Convection term is discretized using a Central Difference Scheme (CDS)
- Diffusion term is discretized using a Central Difference Scheme (CDS)
- Inner derivatives are discretized using Backward Differences (BDS)
- Outer derivatives are discretized using Forward Differences (FDS)
- All systems were solved using Gauss Seidel

## Results
Velocity and Pressure fields converge with grid refinement. The solver produces physically consistent flow behavior and steady-state solutions for a flow with given Re = 100. 

