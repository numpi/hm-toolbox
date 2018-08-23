function [U,S,V] = tsvd(A,tol)
if min(size(A)) == 0
    U = []; S = []; V = []; return;
end
[U,S,V] = svd(A);

t = diag(S);
% t = cumsum(t(end:-1:1));

r = sum(t > tol);

% r = sum(cumsum(diag(S(end:-1:1,end:-1:1))) > tol);
U = U(:,1:r);
V = V(:,1:r);
S = S(1:r,1:r);
