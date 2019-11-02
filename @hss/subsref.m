function varargout = subsref(H, S)
% Overloaded subsref function for HSS matrices

switch S(1).type
    case '()'
        [m,n] = size(H);
        subs = S.subs;
        mind = subs{1}; nind = subs{2};
        if ischar(mind), mind = 1:m; end
        if ischar(nind), nind = 1:n; end
        if ~issorted(mind) | ~issorted(nind),
            error('Indices need to be sorted integers.');
        end
        if min(mind)<1 | max(mind)>m | min(nind)<1 | max(nind)>n,
            error('Indices not in range.');
        end
        %       !!! purge_tree for hss ??? varargout = {purge_tree(hss_sub(H,mind,nind))};
        varargout = {hss_sub(H,mind,nind)};
    otherwise
        varargout = {builtin('subsref',H,S)};
end
end

function H = hss_sub(H, mind, nind);

[m,n] = size(H);

if H.leafnode
    H.D = H.D(mind,nind);
    H.U = H.U(mind,:);
    H.V = H.V(nind,:);
else
    [m1,n1] = size(H.A11);
    mind1 = intersect( 1:m1, mind );
    nind1 = intersect( 1:n1, nind );
    H.A11 = hss_sub(H.A11, mind1, nind1);
    H.ml = length(mind1); H.nl = length(nind1);
    
    [m2,n2] = size(H.A22);
    mind = mind - m1; nind = nind - n1;
    mind2 = intersect( 1:m2, mind );
    nind2 = intersect( 1:n2, nind );
    H.A22 = hss_sub(H.A22, mind2, nind2);
    H.mr = length(mind2); H.nr = length(nind2);
end

end
