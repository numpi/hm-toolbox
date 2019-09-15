%% Computing matrix exponential of banded matrices

%% Approximation with banded matrices
% It is known that the matrix exponential of banded matrix enjoy a decay
% property in the entries, as we move far from the diagonal. That is, 
%
% $$ |[e^{A}]_{ij}| \leq C \rho^{|j - i|},  \qquad i,j = 1, \ldots, n $$
%
% where the constant $C$ and the rate $\rho$ depend on how well the
% exponential is well approximated by polynomials on the spectrum of A (if
% A is normal), or on its field of values (for the general case). Indeed,
% this result holds for a generic matrix function $f(A)$. 
%
% The decay property might deteroriate if $A$ has a large norm, for
% instance when it is obtained by discretizing unbounded differential
% operators. 

%% The matrix under consideration
% As an example, we consider the discretization of the 1D
% Laplacian
%
% $$ \displaystyle
%  A := -\frac{1}{h^2} \left[ \begin{array}{cccc}
%    2  & -1 &  & \\
%    -1 &  2 & \ddots & \\
%    & \ddots & \ddots & -1 \\
%    & & -1 & 2 \\
%  \end{array} \right], \qquad 
%  h := \frac{1}{n - 1}
% $$
%
% This matrix has spectral norm that grows as $n^2$, and we assume that 
% we would like to compute $e^{A}$. 

%% Using HODLR and HSS arithmetic
% Even though, as we will see, the approximate sparse structure is quickly
% lost for large norms $\| A \|$, the matrix $e^A$ has a
% hierarchical rank structure that can be efficiently exploited using HSS
% or HODLR arithmetic. 
%
% HM-TOOLBOX provides a routine EXPM both for the HSS and HODLR classes,
% and this implements a Padè approximant coupled with a scaling and
% squaring strategy, similarly to the MATLAB function EXPM for dense
% matrices. 
%
% We now compute $e^A$ for different choices of $n$. This can be achieved
% in the toolbox as follows. First, let us choose a reasonably large value
% of $n$, but not too large so we can compare with the dense EXPM. 

n = 2048; h = 1 / (n - 1);

%%
% Then, we can define the matrix $A$ as a banded matrix, and using the
% banded constructore in HM-TOOLBOX to obtained an HSS and HODLR
% representation. 

sA = -h^2 \ spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n);
hodlrA = hodlr('banded', sA);
hssA = hss('banded', sA);

%% 
% The toolbox automatically determined the bandwidth of the matrix by
% calling bandwidth, but we could have provided it in case we know it. See 
% ''help hm'' and help hss for further details. 
%
% As usual, the ranks in the offdiagonal blocks can be inspected with a
% call to the SPY command. 

spy(hodlrA);


%% 
% We can now re-compute the matrix exponential by calling EXPM

hodlr_expA = expm(hodlrA);
hss_expA   = expm(hssA);


%% 
% Let us now check the accuracy by comparing with the dense computation. 

eA = expm(full(sA));
res = norm(eA - hodlr_expA)

%% 
% The result is accurate almost up to the default truncation tolerance 
% ($10^{-12}$), but the norm of $e^A$ is rather small, since $A$ is 
% negative definite, and it $e^{\tau A}$ goes to zero for large enough $\tau$.
%
% For a more interesting case, we may consider the exponential of $\tau A$,
% for some $\tau < 1$. For intstance, we may choose $\tau = 1/100$. 

tau = 1 / 100;

%% 
% We can now re-compute the matrix exponential by calling EXPM

hodlr_expA = expm(tau * hodlrA);
hss_expA   = expm(tau * hssA);

%%
% And we may inspect the off-diagonal ranks of the result. This time, we
% try with the HSS one. 

spy(hss_expA);

%% 
% Note that the ranks have increased, but are still moderate, and indeed
% the HSS and HODLR storage are quite effective for this problem. 
%
% Again, we check the residual by comparing with the dense solver. 

eA = expm(tau * full(sA));
res = norm(eA - hss_expA)

%% Sparse vs HODLR / HSS storage
%
% Now that we have compute the matrix exponential, the toolbox allows to
% convert these matrices to sparse ones by thresholding the entries; we can
% try, for instance, to obtain an approximation to $e^A$ accurate to 8
% digits by running:

sp_expA = sparse(hss_expA, 1e-8);

%%
% We may check the bandwidth to see if the resulting matrix is sparse; the
% theory predicts the decay in the entries, but it gets slower as the norm
% increases. 

bandwidth(sp_expA)

%%
% We immediately see that the bandwidth is almost equal to $n$, and
% therefore the sparse storage does not give any advantage for this case. 
%
% Indeed, we can now compute the accuracy and the memory usage of this banded
% approximation, in Megabytes, and we compare it with the memory used by
% the HSS approximation (an analogous result would hold for the HODLR one). 

mem_sparse = whos('sp_expA'); mem_sparse = mem_sparse.bytes / 1024^2
mem_hss    = whos('hss_expA'); mem_hss = mem_hss.bytes / 1024^2
res_sparse = norm(full(sp_expA) - eA)

%%
% On this example, the HSS wins both for storage and accuracy, since the
% residual for the rank structured approximation was order of magnitude
% slower than this one. Such gain would be even stronger as we increase $n$. 
