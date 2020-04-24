function H = hmatrix_from_aca(H, Afun, m, n, progress_fcn)

if isempty(H)
    H = hmatrix_build_default_tree(m, n, hmatrixoption('block-size'));
end

H = hmatrix_from_aca_rec(H, Afun, progress_fcn);

end

function H = hmatrix_from_aca_rec(H, Afun, progress_fcn)
[m, n] = size(H);
if is_leafnode(H)
    progress_fcn(1, 1, 1, 'in progress');
    if H.admissible
        if min(m, n) <= 512
            [U, S, V] = tsvd(Afun(1:m, 1:n), hmatrixoption('threshold'));
            U = U * S;
            H.admissible = true;
            H.U = U;
            H.V = V;
        else
            [H.U, H.V] = aca_or_fail(Afun, m, n, hmatrixoption('threshold'), []);
        end
        if size(H.U, 1) == 0 
            error('Empty representation')
        end
    else
        H.F = Afun((1:m).', (1:n));
    end    
    progress_fcn(1, 1, 1, 'done');
else
    [m1, n1] = size(H.A11);
    
    H.A11 = hmatrix_from_aca_rec(H.A11, Afun, @(l,i,j,s) progress_fcn(l+1,i,j,s));
    H.A12 = hmatrix_from_aca_rec(H.A12, @(i, j) Afun(i, j + n1), @(l,i,j,s) progress_fcn(l+1,i,j+1,s));
    H.A21 = hmatrix_from_aca_rec(H.A21, @(i, j) Afun(i + m1, j), @(l,i,j,s) progress_fcn(l+1,i+1,j,s));
    H.A22 = hmatrix_from_aca_rec(H.A22, @(i, j) Afun(i + m1, j + n1), @(l,i,j,s) progress_fcn(l+1,i+1,j+1,s));
end
end

function [U,S,V] = tsvd(A,tol)
if min(size(A)) == 0
    U = zeros(size(A, 1), 0); S = []; V = zeros(size(A, 2), 0); return;
end
[U,S,V] = svd(full(A));

t = diag(S);
% t = cumsum(t(end:-1:1));

r = sum(t > tol);

% r = sum(cumsum(diag(S(end:-1:1,end:-1:1))) > tol);
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);

end