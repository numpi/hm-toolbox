function eq = check_cluster_equality(H1, H2)
%CHECK_CLUSTER_EQUALITY Check if row and column cluster of H match
%
% R = CHECK_CLUSTER_EQUALITY(H) returns true if the row and column cluster
%    of H match. 
%
% R = CHECK_CLUSTER_EQUALITY(H1, H2) checks if H1 and H2 have the same  row
%     and column clusters, respectively. 
%
% R = CHECK_CLUSTER_EQUALITY(C1, C2) checks if the vectors C1 and C2
%     containing the integer of the integer partitions of the clusters 
%     are equal. 

if isa(H1, 'hodlr') || isa(H1, 'hss') 
    if ~exist('H2', 'var')
        [r, c] = cluster(H1);
        eq = check_cluster_equality(r, c);
    else
        [r1, c1] = cluster(H1);
        [r2, c2] = cluster(H2);

        eq = check_cluster_equality(r1, r2) && ...
             check_cluster_equality(c1, c2);
    end
else
    % In this case H1 and H2 are the actual clusters, that is vectors of 
    % the integer partitions.
    eq = ( length(H1) == length(H2) ) && ...
         ( prod(H1 == H2) == 1 );
end

end

