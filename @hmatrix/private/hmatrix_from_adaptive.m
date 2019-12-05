function [H] = hmatrix_from_adaptive(Afun, m, n, maxrank)
%UNTITLED 

tol = hmatrixoption('threshold');
bs = hmatrixoption('block-size');

H = hmatrix;
H.sz = [m n];

if (m <= bs && n <= bs) || min([m n]) == 1    
    H.F = Afun((1:m).', 1:n);
    H.admissible = false;
else
    % Go down by recursion if necessary. Try first to get a low-rank
    % representation of the full matrix
    [U, V] = aca_or_fail(Afun, m, n, tol, maxrank);
    if ~isempty(U)
        H.admissible = true;
    end
    
    if H.admissible
        H.U = U; H.V = V;
    else
        H.admissible = false;
        
        m1 = ceil(m / 2);
        n1 = ceil(n / 2);
        
        H.A11 = hmatrix_from_adaptive(@(i,j) Afun(i, j), m1, n1, maxrank);
        H.A12 = hmatrix_from_adaptive(@(i,j) Afun(i, j + n1), m1, n - n1, maxrank);
        H.A21 = hmatrix_from_adaptive(@(i,j) Afun(i + m1, j), m - m1, n1, maxrank);
        H.A22 = hmatrix_from_adaptive(@(i,j) Afun(m1 + i, n1 + j), m - m1, n - n1, maxrank);
        
        % Check if the blocks here are all full --- if they are, then we
        % can merge them and avoid having a deep tree for no gain.
        if ~isempty(H.A11.F) && ~isempty(H.A12.F) && ...
                ~isempty(H.A21.F) && ~isempty(H.A22.F)
            %H.F = [ H.A11.F , H.A12.F ; H.A21.F , H.A22.F ];
            %H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
        end
    end            
end

end

