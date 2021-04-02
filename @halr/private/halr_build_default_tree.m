function H = halr_build_default_tree(m, n, block_size) % constructs a default geometrical clustering
H = halr;
H.sz = [m, n];
if max(m, n) > block_size && min(m, n) > 1
    [m1, m2] = split_indices(m);
    [n1, n2] = split_indices(n);
    
    H.A11 = halr_build_default_tree(m1, n1, block_size);
    H.A22 = halr_build_default_tree(m2, n2, block_size);
    H.A12 = offdiagonal_block12(m1, n2, block_size);
    H.A21 = offdiagonal_block21(m2, n1, block_size);
else
    H.F = zeros(m, n);
end
end


function [n1, n2] = split_indices(n)
n1 = ceil(n / 2);
n2 = n - n1;
end

function H = offdiagonal_block12(m, n, block_size)
H = halr;
H.sz = [m, n];
if max(m, n) > block_size && min(m, n) > 1
    [m1, m2] = split_indices(m);
    [n1, n2] = split_indices(n);
    H.A11 = lowrank_block(m1, n1);
    H.A12 = lowrank_block(m1, n2);
    H.A22 = lowrank_block(m2, n2);
    H.A21 = offdiagonal_block12(m2, n1, block_size);
else
    H.F = zeros(m, n);
end
end

function H = offdiagonal_block21(m, n, block_size)
H = halr;
H.sz = [m, n];
if max(m, n) > block_size && min(m, n) > 1
    [m1, m2] = split_indices(m);
    [n1, n2] = split_indices(n);
    H.A11 = lowrank_block(m1, n1);
    H.A21 = lowrank_block(m2, n1);
    H.A22 = lowrank_block(m2, n2);
    H.A12 = offdiagonal_block21(m1, n2, block_size);
else
    H.F = zeros(m, n);
end
end

function H = lowrank_block(m, n)
H = halr;
H.sz = [m, n];
H.U = zeros(m, 0);
H.V = zeros(n, 0);
H.admissible = true;
end


