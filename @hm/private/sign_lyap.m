function X = sign_lyap(varargin)
%SIGN_LYAP Solve the Lyapunov equation with the sign function. 
%
% X = SIGN_LYAP(A, G) solves the Lyapunov equation AX + XA' + G = 0 
% where the the matrix A is negative definite and G is an H-matrix.
%
% This function using the sign function. 

if length(varargin) == 2
    A = varargin{1};
    C = varargin{2};
    is_lyapunov = true;
elseif length(varargin) == 3
    A = varargin{1};
    B = varargin{2};
    C = varargin{3};
    is_lyapunov = false;
else
    error('SIGN_LYAP :: this function requires 2 or 3 arguments');
end

threshold = 1e-8;

% Compute the sign function of the matrix
%   
%  [ A'  G ]
%  [ 0   A ]
% 
converged = false;
it = 1;
max_it = inf;

% The next iteration at which the residue will be checked.
next_check = 3;

while ~converged && it < max_it
    Cold = C;
    Aold = A;
    
    if ~is_lyapunov; Bold = B; end
    
    As = inv(A);
    
    if ~is_lyapunov
        Bs = inv(B);
    end
    
    if it == 1
        if is_lyapunov
            mu = sqrt( norm(As) / norm(A) );
        else
            mu = sqrt(max([ norm(As), norm(Bs) ]) / ...
                      min([ norm(A),  norm(B) ]));
        end
    else
        mu = 1;
    end
    
    if is_lyapunov
        C = .5 * (As * C * As.' / mu + C * mu);
    else
        C = .5 * (As * C * Bs / mu + mu * C);
        B = .5 * (mu * B + Bs / mu);
    end
    
    A = .5 * (mu * A + As / mu);
    
    % We assume that we need at least two steps for convergence. 
    if it == next_check || 1
        if is_lyapunov
            corr = ( norm(A - Aold, 'fro') + norm(C - Cold, 'fro') ) ...
                / ( norm(A, 'fro') + norm(C, 'fro') );
        else
            corr = ( norm(A - Aold, 'fro') + norm(C - Cold, 'fro') + norm(B - Bold, 'fro') ) ...
                / ( norm(A, 'fro') + norm(C, 'fro') + norm(B, 'fro') );
        end
        
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

% The sign of the solution should be chosen accordingly depending if A and
% B converged to the identity or -identity (posdef or negdef coeffs). 
v = rand(size(A, 2), 1);
s = - sign(v' * (A * v));

X = (.5 * s) * C;

end

