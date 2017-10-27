# Hierarchical matrix toolbox
The <code>hm-toolbox</code> is a toolbox implementing the arithmetic of H-matrices in MATLAB. 

Only the simplest case of the so-called HODLR matrices is handled, where the partitioning
is recursively done in 2 x 2 blocks, and the off-diagonal blocks are of low-rank. An implementation
of the HSS arithmetic is also available through the @hss class, although it's still not
as complete as the @hm one.

Routines to compute matrix functions and to solve matrix equations are included. 

Some features depend on external packages, namely chebfun for the construction of 
HODLR / HSS matrices which sample regular functions on a grid, and rktoolbox to
solve Lyapunov / Sylvester equations with D&C methods. 
