function H = hm_build_hm_tree(m, n, block_size, partitionm, partitionn)
    if isempty(partitionm)
        H = build_hm_tree_rec(m, n, block_size);
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
        
        H = build_hm_tree_partition_rec(m, n, partitionm, partitionn);
    end
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

function H = build_hm_tree_rec(m, n, block_size)

H = hm();
    
H.sz = [m n];

if max(m, n) > block_size && min(m, n) > min(8, block_size)
    [m1, m2] = split_indices(m);
    [n1, n2] = split_indices(n);
    
    H.A11 = build_hm_tree_rec(m1, n1, block_size);
    H.A22 = build_hm_tree_rec(m2, n2, block_size);
    
    H.U12 = zeros(m1, 0);
    H.V12 = zeros(n2, 0);
    H.U21 = zeros(m2, 0);
    H.V21 = zeros(n1, 0);
else
    H.F = zeros(m, n);
end

end

function H = build_hm_tree_partition_rec(m, n, partitionm, partitionn)

H = hm();
H.sz = [m n];

lp = length(partitionm);

if lp > 1
    if mod(lp, 2) ~= 0
        error('The partitions vector must have length equal to a power of 2');
    end
    
    m1 = partitionm(lp/2);
    n1 = partitionn(lp/2);
    
    m2 = m - m1;
    n2 = n - n1;
    
    if m1 == 0 || m2 == 0 || n1 == 0 || n2 == 0
        H.F = zeros(m, n);
    else
        H.A11 = build_hm_tree_partition_rec(m1, n1, partitionm(1:lp/2), partitionn(1:lp/2));
        H.A22 = build_hm_tree_partition_rec(m2, n2, partitionm(lp/2+1:end) - m1, partitionn(lp/2+1:end) - n1);
        
        H.U12 = zeros(m1, 0); H.V12 = zeros(n2, 0);
        H.U21 = zeros(m2, 0); H.V21 = zeros(n1, 0);
    end
else
   H.F = zeros(m, n); 
end
    
end

function [n1, n2] = split_indices(n)
n1 = ceil(n / 2);
n2 = n - n1;
end