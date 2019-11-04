%% Examples of usage for hodlr-toolbox

%% Introduction
% These pages provide a list of possible usages for the hm-toolbox. Some
% examples apply both to the HSS and HODLR storage, some others are only
% given for one of the two; note that the functionalities are largely
% equivalent between the formats, so the steps in one format can usually be
% repeated in the other with no changes.
%
% A notable exception is in the 'handle' constructor, which uses Adaptive
% Cross Approximation for the HODLR format, and randomized sampling in the
% HSS one. The latter also requires a fast matrix vector multiplication,
% which the former doesn't.
%
% This difference can often be overcome relying on the conversion between
% the formats provided by the HM2HSS and HSS2HM functions.
%

%% Linear systems
%
% * <hodlr_linear_system.html Solving a linear system> with rank structured
%      matrix, using the HODLR format.
% * <hss_toeplitz_solver.html Solving a Toeplitz linear system> shows how
%     to construct a superfast Toeplitz solver using HSS matrices.

%% Matrix functions
% * <hm_expm.html Computing the matrix exponential> showcases
%      the functions EXPM implemented in the toolbox; we see that in certain
%      cases hierarchical formats can be superior to sparse approximation
%      for the matrix exponential, despite the guaranteed exponential decay
%      of the entries far from the diagonal.

%% Matrix equations
% * <hodlr_lyapunov.html Solving linear matrix equations> arising from 2D
%      discretization of PDES, exploiting the HODLR format in the coefficients.
% * <hss_lyapunov.html Solving linear matrix equations>, concerning the same
%      PDEs of the previous item, but using the HSS format.