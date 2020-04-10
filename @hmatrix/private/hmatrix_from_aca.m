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
        [H.U, H.V] = aca_or_fail(Afun, m, n, hmatrixoption('threshold'), []);
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
