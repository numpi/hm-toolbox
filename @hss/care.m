function X = care(A, B, C, spA, tol, debug, nrmtype)
%CARE computes the stabilizing solution X of
%    the continuous-time algebraic Riccati equation
%
%   A'X + XA - XBB'X + C = 0
%
% X = CARE(A, B, C) solves the CARE equation using a divide and conquer
%     method with the Extended Krylov for treating equations with rhs of low-rank. 
%     If a sparse version spA of A is provided then sparsity is exploited for
%     building Krylov subspaces
%
if ~exist('debug', 'var') || isempty(debug)
	debug = false;
end
if ~exist('tol', 'var') || isempty(tol)
	tol = hssoption('threshold');
end
if ~exist('nrmtype', 'var') || isempty(nrmtype)
	nrmtype = 2;
end
if ~exist('spA', 'var') 
	spA = [];
end
	
X = hss_dac_care(A, B, C, spA, tol, debug, nrmtype);
end

