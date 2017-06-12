function X = sign_lyap(A, G)
%SIGN_LYAP Solve the Lyapunov equation with the sign function. 
%
% X = SIGN_LYAP(A, G) solves the Lyapunov equation AX + XA' + G = 0 
% where the the matrix A is negative definite and G is an H-matrix.
%
% This function using the sign function. 

threshold = 1e-6;

% Compute the sign function of the matrix
%   
%  [ A'  G ]
%  [ 0  -A ]
% 
converged = false;
it = 1;
max_it = inf;

n = size(A, 1);

while ~converged && it < max_it
    Gold = G;
    Aold = A;
    
    As = inv(A);
	
	mu = sqrt( norm(As) / norm(A) );
	
    G = .5 * (As' * G * As + G) / mu;
    A = .5 * (mu * A + As / mu);
    
    corr = ( norm(A - Aold) + norm(G - Gold) ) / ( norm(A) + norm(G) );
    converged = corr < threshold;  
	
    % Enable for sign method debugging
    % fprintf ('SIGN :: Iteration %d, Relative correction: %e\n', it, corr);  
    it = it + 1;
end

X = .5 * G;

end

