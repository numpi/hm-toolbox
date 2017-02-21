function [U, V] = compress_matrix(Uold)
%COMPRESS_MATRIX Get a low-rank representation for Uold.  

threshold = eps;

[U,S,V] = svd(Uold);

rk = sum(diag(S) > threshold);

U = U(:,1:rk) * sqrt(S(1:rk,1:rk));
V = V(:,1:rk) * sqrt(S(1:rk,1:rk));

end

