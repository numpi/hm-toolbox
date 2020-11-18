function [X, Y] = sqrtm(A)
%SQRTM Compute the matrix square root of X.
%
% Y = SQRT(X) computes a matrix Y such that Y^2 = X. The computed square
% root is the principal branch of the square root.

converged = false;

threshold = hssoption('threshold') * norm(A);

X = A;
Y = eye(size(A, 1), 'like', A);

Xnew = .5 * (X + Y);
Y = .5 * (Y + inv(X));
X = Xnew;

it = 1;
max_its = 100;

while ~converged && it < max_its
    it = it + 1;
    
    Xnew = .5 * (X + inv(Y));
    Ynew = .5 * (Y + inv(X));
    
    if norm(Xnew - X) < threshold
        converged = true;
    end
    
    X = Xnew;
    Y = Ynew;
end

if it == max_its
    warning('Maximum number of iterations exceeded');
end

end

