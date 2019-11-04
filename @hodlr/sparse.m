function A = sparse(H, tol)
%SPARSE Converts a HODLR matrix into a sparse matrix.
%
% S = sparse(H) converts a HODLR matrix (with hopefully many zero
%     entries) into a sparse matrix.
%
% S = sparse(H, tol) regards all entries of absolute value not larger
%     than the tolerance tol as zero.

if nargin < 2,
    tol = 0;
end
if nargin == 2 & tol<0,
    error('Tolerance tol must be nonnegative.');
end

if is_leafnode(H),
    if tol == 0,
        A = sparse(H.F);
    else
        [i, j, vals] = find( ( abs(H.F)>tol ).*H.F );
        [m, n] = size(H.F);
        A = sparse(i, j, vals, m, n);
    end
else
    A = [ sparse(H.A11, tol), lr2sparse(H.U12, H.V12, tol);
        lr2sparse(H.U21, H.V21, tol), sparse(H.A22, tol) ];
end

