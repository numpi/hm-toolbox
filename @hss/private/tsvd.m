function [U,S,V] = tsvd(A,tol)
if min(size(A)) == 0
	U = []; S = []; V = []; return;
end
[U,S,V] = svd(A);
r = sum(diag(S) > tol * S(1,1));
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);
