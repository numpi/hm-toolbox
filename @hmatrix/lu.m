function [L, U] = lu(H)
%LU Compute the LU factorization of H

L = hmatrix;
U = hmatrix;

L.sz = [ H.sz(1), min(H.sz) ];
U.sz = [ min(H.sz), H.sz(2) ];

if is_leafnode(H)
    if H.admissible
        error('Admissible block on the diagonal during LU factorization');
    else
        [L.F, U.F] = lu(H.F);
    end
else
    [m1, n1] = size(H.A11);
    [m2, n2] = size(H.A22);
    
    [L.A11, U.A11] = lu(H.A11);
    U.A12 = L.A11 \ H.A12;
    L.A21 = H.A21 / U.A11;
    % [L.A22, U.A22] = lu(H.A22 - L.A21 * U.A12);
    [L.A22, U.A22] = lu(compress_hmatrix(hmatrix_minus(H.A22, L.A21 * U.A12), []));
    
    % Since these matrices are lower and upper triangular we can set some
    % blocks to be empty, and in particular admissible.
    L.A12 = hmatrix;
    L.A12.sz = [m1, n2];
    L.A12.admissible = true;
    L.A12.U = zeros(m1, 0); L.A12.V = zeros(n2, 0);
    U.A21 = hmatrix;
    U.A21.sz = [m2 n1];
    U.A21.admissible = true;
    U.A21.U = zeros(m2, 0); U.A21.V = zeros(n1, 0);
end

end

