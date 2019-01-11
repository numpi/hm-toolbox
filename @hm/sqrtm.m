function X = sqrtm(H)
%SQRTM Compute the matrix square root of H.
%
% X = SQRT(H) computes a matrix Y such that X^2 = H. The computed square
% root is the principal branch of the square root.

converged = false;

threshold = hmoption('threshold') * norm(H);
        
X = H;
Y = eye(size(H), 'like', H);

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

