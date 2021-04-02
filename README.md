# Hierarchical matrix toolbox

[![tests](https://github.com/numpi/hm-toolbox/actions/workflows/automatic-testing.yml/badge.svg)](https://github.com/numpi/hm-toolbox/actions/workflows/automatic-testing.yml)

The <code>hm-toolbox</code> is a toolbox implementing the arithmetic of HODLR and HSS matrices in MATLAB. 

The HODLR case is handled in the @hodlr class, and correspond to H-matrices with partitioning
recursively done in 2 x 2 blocks, where the off-diagonal blocks are of low-rank. The HSS 
arithmetic uses the same partitioning (with nested bases), and is available through the
@hss class. 

Routines to compute matrix functions [1] and to solve matrix equations are included [1,2]. 

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

You are now ready to create new @hodlr and @hss objects. Check some examples in the
<code>examples/</code> folder. 

## Octave compatibility
The toolbox is compatible with GNU Octave >= 4.4. The easiest way to get a
recent version of Octave on Linux is to enable Flatpak following the instructions
[here](https://flatpak.org/setup/), and then run the command 
```
 flatpak install flathub org.octave.Octave
```
You can then run Octave from the menu or by typing <code>flatpak run org.octave.Octave</code> in a terminal. 
Some functions require dependencies found in additional packages, namely <code>octave-control</code>
(for the solvers of Lyapunov and Sylvester equations) and <code>octave-statistics</code> (for
the ACA code that is used inside the HODLR handle constructor). 

# References

1. Massei, S., Palitta, D., & Robol, L. (2018). Solving Rank-Structured Sylvester and Lyapunov Equations. SIAM Journal on Matrix Analysis and Applications, 39(4), 1564-1590.
2. Kressner, D., Massei, S., & Robol, L. (2019). Low-rank updates and a divide-and-conquer method for linear matrix equations. SIAM Journal on Scientific Computing 41 (2), A848-A876.
3. Massei, S., Robol, L., & Kressner, D. (2020). hm-toolbox: MATLAB Software for HODLR and HSS Matrices. SIAM Journal on Scientific Computing, 42(2), C43-C68.
