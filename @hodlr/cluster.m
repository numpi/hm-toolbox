function [p, q] = cluster(H)
%CLUSTER Obtain the partitioning of the indices at the lowest level.
%
% [R, C] = CLUSTER(H) obtains a description of the clustering of the matrix
%     H as an array of indices [P(1), ..., P(L)] such that the partitioning
%     of the matrix at the lowest level is:
%
%             (1, P(1))    (P(1), P(2))   ...   (P(L-1), P(L))
%
%     The vector R encodes the partitioning of the rows, while C encodes
%     the one of the columns.

if is_leafnode(H)
    p = size(H, 1);
    q = size(H, 2);
else
    [p1, q1] = cluster(H.A11);
    [p2, q2] = cluster(H.A22);
    
    % Add elements in the case the tree is slightly unbalanced
    m = max(length(p1), length(p2));
    
    p1 = [ p1, ones(1, m - length(p1)) * p1(end) ];
    p2 = [ p2, ones(1, m - length(p2)) * p2(end) ];
    q1 = [ q1, ones(1, m - length(q1)) * q1(end) ];
    q2 = [ q2, ones(1, m - length(q2)) * q2(end) ];
    
    p = [ p1, p2 + p1(end) ];
    q = [ q1, q2 + q1(end) ];
end

end

