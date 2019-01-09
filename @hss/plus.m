function H = plus(H1, H2)

if ~check_cluster_equality(H1, H2)    
    error('H1 + H2: Cluster or dimension mismatch in H1 and H2');
end

H = hss_sum(H1, H2);

end
