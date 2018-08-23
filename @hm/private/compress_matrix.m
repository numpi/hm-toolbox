function [U,V] = compress_matrix(A)
%COMPRESS_MATRIX Compress a dense matrix into low rank format.

maxk = 64;

if max(size(A)) < 2048
    [U,S,V] = svd(A);
else
    [U,S,V] = svds(A, maxk);
end

threshold = hmoption('threshold');

k = sum(abs(diag(S)) > S(1,1) * threshold);

U = U(:,1:k) * S(1:k,1:k);
V = V(:,1:k);


end

