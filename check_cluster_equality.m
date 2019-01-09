function eq = check_cluster_equality(H1, H2)
%CHECK_CLUSTER_EQUALITY Check if row and column cluster of H match
%
% R = CHECK_CLUSTER_EQUALITY(H) returns true if the row and column cluster
%    of H match. 
%
% R = CHECK_CLUSTER_EQUALITY(H1, H2) checks if H1 and H2 have the same  row
%     and column clusters, respectively. 

if ~exist('H2', 'var')
    [r, c] = cluster(H1);

    if length(r) ~= length(c)
        eq = false;
        return;
    end

    eq = prod(r == c) == 1;
else
    [r1, c1] = cluster(H1);
    [r2, c2] = cluster(H2);

    if length(r1) ~= length(r2) || length(c1) ~= length(c2)
        eq = false;
        return;
    end

    eq = (( prod(r1 == r2)  ) * ( prod(c1 == c2)  )) == 1;
end

end

