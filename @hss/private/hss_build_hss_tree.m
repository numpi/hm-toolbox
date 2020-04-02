function H = hss_build_hss_tree(m, n, block_size, partitionm, partitionn)
if isempty(partitionm)
    H = build_hss_tree_rec(m, n, block_size);
else
    check_partitioning(partitionm, m);
    
    if ~exist('partitionn', 'var')
        partitionn = partitionm;
    else
        check_partitioning(partitionn, n);
    end
    
    
    if length(partitionn) ~= length(partitionm)
        error('Row and column partitioning must have the same number of elements');
    end
    
    H = build_hss_tree_partition_rec(m, n, partitionm, partitionn);
end

H.topnode = 1;
end

function check_partitioning(partition, n)
for j = 2 : length(partition)
    if partition(j) < partition(j-1)
        error('The cluster vector must be non decreasing');
    end
end

if partition(end) ~= n
    error('The last element of the cluster must be the dimension');
end
end


function H = build_hss_tree_rec(m, n, block_size)

H = hss();
H.topnode  = 0;

if max(m, n) > block_size && min(m, n) > 1
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

function H = build_hss_tree_partition_rec(m, n, partitionm, partitionn)

H = hss();
H.topnode = 0;

lp = length(partitionm);

if lp > 1
    if mod(lp, 2) ~= 0
        error('The partitions vector must have length equal to a power of 2');
    end
    
    H.ml = partitionm(lp/2);
    H.nl = partitionn(lp/2);
    
    H.mr = m - H.ml;
    H.nr = n - H.nl;
    
    if (H.ml > 0 && H.nl > 0) && (H.mr > 0 && H.nr > 0)
        H.leafnode = 0;
        H.A11 = build_hss_tree_partition_rec(H.ml, H.nl, partitionm(1:lp/2), partitionn(1:lp/2));
        H.A22 = build_hss_tree_partition_rec(H.mr, H.nr, partitionm(lp/2+1:end) - H.ml, partitionn(lp/2+1:end) - H.nl);
    elseif (H.ml > 0 && H.nl > 0) && (H.mr == 0 || H.nr == 0)
        H = build_hss_tree_partition_rec(H.ml, H.nl, partitionm(1:lp/2), partitionn(1:lp/2));
    elseif (H.ml == 0 || H.nl == 0) && (H.mr > 0 && H.nr > 0)        
        H = build_hss_tree_partition_rec(H.mr, H.nr, partitionm(lp/2+1:end) - H.ml, partitionn(lp/2+1:end) - H.nl);
    else
        H.D = zeros(m, n);
        H.U = zeros(m, 0);
        H.V = zeros(n, 0);
        H.leafnode = 1;
    end
else
    H.D = zeros(m, n);
    H.leafnode = 1;
end

end

function [n1, n2] = split_indices(n)
n1 = ceil(n / 2);
n2 = n - n1;
end
