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

# Installation instructions

To install the toolbox download the latest revision from Github by running
```
 git clone https://github.com/numpi/hm-toolbox.git
```

or downloading the ZIP file from the webpage [github.com/numpi/hm-toolbox](https://github.com/numpi/hm-toolbox). 
Rename the folder to <code>hm-toolbox</code> if needed. Then, add it to your MATLAB path by running
```Matlab
 >> addpath /path/to/hm-toolbox
```

You are now ready to create new @hm and @hss objects. Check some examples in the
<code>examples/</code> folder. 
