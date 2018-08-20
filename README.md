# Hierarchical matrix toolbox

[![Build status](https://api.travis-ci.org/numpi/hm-toolbox.svg?branch=master)](https://travis-ci.org/numpi/hm-toolbox)

The <code>hm-toolbox</code> is a toolbox implementing the arithmetic of H-matrices in MATLAB. 

Only the simplest case of the so-called HODLR matrices is handled, where the partitioning
is recursively done in 2 x 2 blocks, and the off-diagonal blocks are of low-rank. An implementation
of the HSS arithmetic is also available through the @hss class. 

Routines to compute matrix functions [1] and to solve matrix equations are included [1,2]. 

Some features depend on external packages, namely chebfun for the construction of 
HODLR / HSS matrices which sample (piecewise) regular functions on a grid, and rktoolbox to
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

# References

[1] Massei, S., Palitta, D., & Robol, L. (2017). Solving rank structured Sylvester and Lyapunov equations. arXiv preprint arXiv:1711.05493.

[2] Kressner, D., Massei, S., & Robol, L. (2017). Low-rank updates and a divide-and-conquer method for linear matrix equations. arXiv preprint arXiv:1712.04349.
