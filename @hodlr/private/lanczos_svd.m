function [U, S, V] = lanczos_svd(Afun, m, n, tol)
%LANCZOS_SVD Two-sided Lanczos for SVD.

V = randn(n, 1);
V(:,1) = V(:,1) / norm(V);

U = zeros(m, 0);

res = inf;

if ~exist('tol', 'var')
    tol = eps;
end

alfa = [];
beta = [];

it = 0;

nrm = 0.0;
nrm_converged = false;

while res > tol
    it = it + 1;
    
    if it > min(m,n) / 2
        if m < n
            A = Afun(eye(m), 'transp')';
        else
            A = Afun(eye(n), 'notransp');
        end
        
        [U, S, V] = svd(A, 'econ');
        rk = sum(diag(S) > tol * S(1,1));

        U = U(:,1:rk);
        V = V(:,1:rk);
        S = S(1:rk, 1:rk);
        return;
    end
    
    w = Afun(V(:,end), 'notransp');
    
    % Reorthogonalize w against U
    w = w - U * (U' * w);
    w = w - U * (U' * w);
    
    alfa = [ alfa, norm(w) ];

    if alfa(end) == 0
        break;
    end
    
    U = [ U , w / alfa(end) ];
    
    % Compute beta
    w = Afun(U(:,end), 'transp');
    w = w - V * (V' * w);
    w = w - V * (V' * w);
    
    beta = [ beta, norm(w) ];
    
    if beta(end) == 0
        break;
    end
    
    V = [ V, w / beta(end) ];
    
    % Estimate the norm: if we have a good estimate, evaluate the
    % possibility of stopping the iteration.
    if ~nrm_converged
        nrm_est = norm(diag(alfa) + diag(beta(1:end-1), 1));
        if (nrm_est - nrm) < tol * nrm
            nrm_converged = true;
        end
        nrm = nrm_est;
    else
        res = max(abs([ alfa(end), beta(end) ])) / nrm;
    end
    
    % fprintf('Iteration %d, res = %e, nrm_conv = %d, nrm = %e\n', it, res, nrm_converged, nrm);
end

if length(alfa) == length(beta)
    [Ul, S, Vl] = svd(diag(alfa) + diag(beta(1:end-1), 1));
else
    [Ul, S, Vl] = svd([ diag(alfa(1:end-1)), zeros(length(beta),1)  ] + ...
                      [ zeros(length(beta), 1), diag(beta) ], 'econ');
end

if isempty(S)
    U = zeros(m, 0);
    S = [];
    V = zeros(n, 0);
    return;
end

% Possibly perform recompression
rk = sum(diag(S) > tol * S(1,1));

U = U * Ul(:,1:rk);
V = V(:,1:size(Vl, 1)) * Vl(:,1:rk);
S = S(1:rk, 1:rk);

end

