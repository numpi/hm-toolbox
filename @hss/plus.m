function C = plus(A,B)

if ~check_cluster_equality(A, B)    
    error('Cluster or dimension mismatch in A and B');
end

C = hss_sum(A, B);
end
