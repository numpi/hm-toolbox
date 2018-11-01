function H = hss_build_hss_tree(m, n, block_size)
    H = build_hss_tree_rec(m, n, block_size);
    H.topnode = 1;
end

function H = build_hss_tree_rec(m, n, block_size)
H = hss();

%if m ~= n
%    error('Rectangular HSS matrices are not supported');
%end

H.topnode  = 0;

if max(m, n) > block_size && min(m, n) > min(8, block_size)
    [m1, m2] = split_indices(m);
    [n1, n2] = split_indices(n);
    
    H.ml = m1; H.mr = m2;
    H.nl = n1; H.nr = n2;
    
    H.A11 = build_hss_tree_rec(m1, n1, block_size);
    H.A22 = build_hss_tree_rec(m2, n2, block_size);
    
    H.leafnode = 0;
else
    H.leafnode = 1;
    H.D = zeros(m, n);
    H.U = zeros(m, 0);
    H.V = zeros(n, 0);
end
end

function [n1, n2] = split_indices(n)
n1 = ceil(n / 2);
n2 = n - n1;
end
