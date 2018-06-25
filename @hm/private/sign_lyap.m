function X = sign_lyap(A, G)
%SIGN_LYAP Solve the Lyapunov equation with the sign function. 
%
% X = SIGN_LYAP(A, G) solves the Lyapunov equation AX + XA' + G = 0 
% where the the matrix A is negative definite and G is an H-matrix.
%
% This function using the sign function. 

threshold = 1e-8;

% We need to change sign to the right hand side
G = -G;

% Compute the sign function of the matrix
%   
%  [ A'  G ]
%  [ 0  -A ]
% 
converged = false;
it = 1;
max_it = inf;

n = size(A, 1);

% The next iteration at which the residue will be checked.
next_check = 3;

while ~converged && it < max_it
    Gold = G;
    Aold = A;
    
    As = inv(A);
if it == 1
    mu = sqrt( norm(As) / norm(A) );
else
    mu=1;
end
    G = .5 * (As * G * As.' / mu + G * mu);
    A = .5 * (mu * A + As / mu);        

    
    % We assume that we need at least two steps for convergence. 
    if it == next_check || 1
        corr = ( norm(A - Aold, 'fro') + norm(G - Gold, 'fro') ) ...
            / ( norm(A, 'fro') + norm(G, 'fro') );
        converged = corr < threshold;  
        
        % Newton is quadratically convergent, so if the above gives a good
        % estimate of the error we have a good guess of the next time we'll
        % need to check the norm. 
        next_check = it + max(1, floor(log2(log2(corr / threshold))));
	end
	
    % Enable for sign method debugging
    % fprintf ('SIGN :: Iteration %d, Relative correction: %e\n', it, corr);  
    % fprintf('SIGN :: Iteration %d, hmrank %d\n', it, hmrank(A));
    it = it + 1;
end

X = .5 * G;

end

