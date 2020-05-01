function [H] = hmatrix_from_adaptive(Afun, m, n, maxrank, support, nrm)
%HMATRIX_FROM_ADAPTIVE Construct an hmatrix of a function sampling
%
% Inputs:
%
%  - AFUN: Handle function such that AFUN(I,J) is the submatrix of
%      samplings corresponding to indices in I and J. 
%  - M, N Dimensions of A.
%  - MAXRANK: Maximum allowed ranks for admissible blocks; if the rank is
%      higher, then the block is recursively split. This argument is
%      optional, if omitted or [] then MAXRANK = MIN(M,N) / 10.
%  - SUPPORT: Optional handle function SUPPORT(I1, J1, I2, J2) that
%      identifies if the submatrix AFUN(I1:I2, J1:J2) is either inside,
%      outside, or across the support. This function is expected to return
%      0, 1, or 2 in these cases, respectively. Default implementation
%      returns 1. 
%  

% The variable norm can be used to track the norm used to perform relative
% truncation, and is an estimate of the norm of the entire matrix. At the
% first call, is set to zero and then estimated by the first call to
% aca_or_fail.
if ~exist('nrm', 'var')
    nrm = 0;
end

tol = hmatrixoption('threshold');
bs = hmatrixoption('block-size');

if ~exist('support', 'var')
    support = @(i1, j1, i2, j2) 1;
end

compression = 'svd';

if isempty(maxrank) || maxrank > round(min(m,n) / 16)
    maxrank = round(min(m,n) / 16);
end

H = hmatrix;
H.sz = [m n];

if (m <= bs && n <= bs) || min([m n]) == 1
    switch compression
        case 'svd'
            [U, S, V] = tsvd(Afun(1:m, 1:n), tol, nrm);
            U = U * S;
            if size(U, 2) < maxrank
                H.admissible = true;
                H.U = U;
                H.V = V;
            end
        case 'aca'
            [U, V] = aca_or_fail(Afun, m, n, tol, maxrank, nrm);
            if size(U, 1) ~= 0
                H.admissible = true;
                H.U = U; H.V = V;
            end
    end
    
    if ~H.admissible    
        H.F = Afun((1:m).', 1:n);
    end
else
    % Depending on the support, we may choose to split the blocks anyway
    switch support(1, 1, m, n)
        case 0
            H.admissible = true;
            H.U = zeros(m, 0); H.V = zeros(n, 0);
        case 1
            compression_method = 'aca';
            
            switch compression_method
                case 'aca'            
                    [U, V, newnrm] = aca_or_fail(Afun, m, n, tol, maxrank, nrm);
                    if size(U, 1) ~= 0
                        H.admissible = true;
                        H.U = U; H.V = V;
                    else
                        H.admissible = false;
                        nrm = max(nrm, newnrm);
                    end
                case 'svd'
                    [U, S, V] = tsvd(Afun(1:m, 1:n), tol);
                    U = U * S;
                    if size(U, 2) < min(maxrank, min(m, n) / 16)
                        H.admissible = true;
                        H.U = U;
                        H.V = V;
                    else
                        H.admissible = false;
                    end
            end
        case 2
            H.admissible = false;
    end
            
    % Go down by recursion if necessary. Try first to get a low-rank
    % representation of the full matrix
    if ~H.admissible
        
        m1 = ceil(m / 2);
        n1 = ceil(n / 2);
        
        H.A11 = hmatrix_from_adaptive(@(i,j) Afun(i, j), m1, n1, ...
            maxrank, support, nrm);
        H.A12 = hmatrix_from_adaptive(@(i,j) Afun(i, j + n1), m1, ...
            n - n1, maxrank, ...
            @(i1, j1, i2, j2) support(i1, j1 + n1, i2, j2 + n1), nrm);
        H.A21 = hmatrix_from_adaptive(@(i,j) Afun(i + m1, j), ...
            m - m1, n1, maxrank, ...
            @(i1, j1, i2, j2) support(i1 + m1, j1, i2 + m1, j2), nrm);
        H.A22 = hmatrix_from_adaptive(@(i,j) Afun(m1 + i, n1 + j), ...
            m - m1, n - n1, maxrank, ...
            @(i1, j1, i2, j2) support(i1 + m1, j1 + n1, i2 + m1, j2 + n1), nrm);
        
        % Check if the blocks here are all full --- if they are, then we
        % can merge them and avoid having a deep tree for no gain.
        if ~isempty(H.A11.F) && ~isempty(H.A12.F) && ...
                ~isempty(H.A21.F) && ~isempty(H.A22.F)
            H.F = [ H.A11.F , H.A12.F ; H.A21.F , H.A22.F ];
            H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
        end
        
        if ~is_leafnode(H) && (H.A11.admissible && H.A12.admissible && H.A21.admissible && H.A22.admissible)
            % if we have 4 children of low-rank let's try to merge them
            % first we compute the low-rank representation of the 4 blocks together
            [U, V] = merge_low_rank(H.A11, H.A12, H.A21, H.A22, nrm);

            if size(U, 2) < min([size(H) / 2, maxrank]) % if the rank is not too high...
                mem = sum(H.sz) * size(U, 2);
                mem11 = sum(H.A11.sz) * size(H.A11.U, 2);
                mem12 = sum(H.A12.sz) * size(H.A12.U, 2);
                mem21 = sum(H.A21.sz) * size(H.A21.U, 2);
                mem22 = sum(H.A22.sz) * size(H.A22.U, 2);

                if mem <= 2 * (mem11 + mem12 + mem21 + mem22) % ...and we do not lose too much in memory, then we merge them
                    H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
                    H.U = U;
                    H.V = V;
                    H.admissible = true;
                end
            end
        end
    end            
end

end

function [U,S,V] = tsvd(A,tol,nrm)
if min(size(A)) == 0
    U = zeros(size(A, 1), 0); S = []; V = zeros(size(A, 2), 0); return;
end
[U,S,V] = svd(full(A));

t = diag(S);
% t = cumsum(t(end:-1:1));

if nrm == 0
    nrm = t(1);
end

r = sum(t > tol * nrm);

% r = sum(cumsum(diag(S(end:-1:1,end:-1:1))) > tol);
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);

end

function [U, V] = merge_low_rank(A11, A12, A21, A22, nrm)
[U, V] = compress_factors(blkdiag([ A11.U, A12.U ], [ A21.U, A22.U ]), ...
    [ blkdiag(A11.V, A12.V), blkdiag(A21.V, A22.V) ], nrm);
end
