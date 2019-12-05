function varargout = subsref(H, S)
% Overloaded subsref function for HMATRIX matrices

switch S(1).type
    case '()'
        [m,n] = size(H);
        subs = S.subs;
        mind = subs{1}; nind = subs{2};
        if ischar(mind), mind = 1:m; end
        if ischar(nind), nind = 1:n; end
        if ~issorted(mind) || ~issorted(nind)
            error('Indices need to be sorted integers.');
        end
        if min(mind)<1 || max(mind)>m || min(nind)<1 || max(nind)>n
            error('Indices not in range.');
        end
        varargout = {purge_tree(hmatrix_sub(H,mind,nind))};
    otherwise
        varargout = {builtin('subsref',H,S)};
end
end

function H = hmatrix_sub(H, mind, nind)

H.sz = [length(mind), length(nind)];

if is_leafnode(H)
    if H.admissible
        H.U = H.U(mind, :);
        H.V = H.V(nind, :);
    else
        H.F = H.F(mind,nind);
    end
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    mind1 = my_intersect(m1, mind);
    nind1 = my_intersect(n1, nind);
    
    mind = mind - m1; nind = nind - n1;
    mind2 = my_intersect(H.A22.sz(1), mind);
    nind2 = my_intersect(H.A22.sz(2), nind);
    
    
    H.A11 = hmatrix_sub(H.A11, mind1, nind1);
    H.A22 = hmatrix_sub(H.A22, mind2, nind2);    
    H.A12 = hmatrix_sub(H.A12, mind1, nind2);
    H.A21 = hmatrix_sub(H.A21, mind2, nind1);
end

end

function s = my_intersect(m, mind)
    s = mind(mind <= m);
    s = s(s > 0);
end

