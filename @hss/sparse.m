function A = sparse(H, tol)
%SPARSE Convert an HSS matrix to the sparse format.
%
% A = SPARSE(H) is a sparse representation of the HSS matrix A.
%
% A = SPARSE(H, TOL) specified a tolerance for theresholding the entries.

if nargin < 2
    tol = 0;
end

A = sparse(hss2hodlr(H), tol);

end

