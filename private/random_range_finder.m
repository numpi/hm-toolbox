function [Q, X, Omega] = random_range_finder(Afun, n, tol)
%RANDOM_SVD Compute an approximate SVD with random sampling.


if ~exist('tol', 'var')
    tol = eps;
end

if ~isa(Afun, 'function_handle')
    Afun = @(v) Afun * v;
end

% Oversampling parameter
p = 10;

% Estimated rank
k = 5;

converged = false;
Omega = randn(n, p);
X = Afun(Omega);

while ~ converged
    Omega_temp = randn(n, k);
    Omega = [Omega, Omega_temp];
    X = [X, Afun(Omega_temp) ];
    [Q, R, ~] = qr(X, 0);
    
    if abs(R(size(R,2),end)) / abs(R(1,1)) < tol
        converged = true;
        
        rk = sum(diag(abs(R)) > tol * abs(R(1,1)));
        
        Q = Q(:,1:rk);
    end
end


end

