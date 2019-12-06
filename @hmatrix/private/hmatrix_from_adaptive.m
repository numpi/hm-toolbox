function [H] = hmatrix_from_adaptive(Afun, m, n, maxrank)
%UNTITLED 

tol = hmatrixoption('threshold');
bs = hmatrixoption('block-size');

compression = 'svd';

if isempty(maxrank)
    maxrank = min(m,n) / 10;
end

H = hmatrix;
H.sz = [m n];

if (m <= bs && n <= bs) || min([m n]) == 1
    switch compression
        case 'svd'
            [U, S, V] = tsvd(Afun(1:m, 1:n), tol);
            U = U * S;
            if size(U, 2) < min(maxrank, min(m, n) / 2)
                H.admissible = true;
                H.U = U;
                H.V = V;
            end
        case 'aca'
            [U, V] = aca_or_fail(Afun, m, n, tol, maxrank);
            if size(U, 1) ~= 0
                H.admissible = true;
                H.U = U; H.V = V;
            end
    end
    
    if ~H.admissible    
        H.F = Afun((1:m).', 1:n);
    end
else
    % Go down by recursion if necessary. Try first to get a low-rank
    % representation of the full matrix
    switch compression
        case 'aca'
            [U, V] = aca_or_fail(Afun, m, n, tol, maxrank);
            if size(U, 1) ~= 0
                H.admissible = true;
                H.U = U; H.V = V;
            end
        case 'svd'
            [U, S, V] = tsvd(Afun(1:m, 1:n), tol);
            U = U * S;
            if size(U, 2) < min(maxrank, min(m, n) / 2)
                H.admissible = true;
                H.U = U;
                H.V = V;
            end
    end
    
    if ~H.admissible
        
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

function [U,S,V] = tsvd(A,tol)
if min(size(A)) == 0
    U = zeros(size(A, 1), 0); S = []; V = zeros(size(A, 2), 0); return;
end
[U,S,V] = svd(A);

t = diag(S);
% t = cumsum(t(end:-1:1));

r = sum(t > tol);

% r = sum(cumsum(diag(S(end:-1:1,end:-1:1))) > tol);
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);

end
