function [rc, cc] = cluster(T)
%CLUSTER 

rc = arrayfun(@(i) size(T.U{i}, 1), 1 : length(T.U));
cc = arrayfun(@(i) size(T.V{i}, 1), 1 : length(T.V));

end

