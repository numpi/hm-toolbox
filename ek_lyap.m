function [Xu, VA, D, it] = ek_lyap(A, u, k, tol, debug, varargin)
%EK_LYAP Approximate the solution of a AX + XA' + u * u' = 0
%
% XU = EK_LYAP(A, U, K) approximates the solution of the Lyapunov equation
%     in the factored form XU * XU' in the case of X positive definite.
%
% [XU, VA] = EK_LYAP(A, U, K, TOL, DEBUG) also returns the basis VA, and
%     the optional parameters TOL and DEBUG control the stopping criterion
%     and the debugging during the iteration.

% [XU, VA, D, it] = EK_LYAP(A, U, K, TOL, DEBUG, 'kernel', K) solves A X + X A' + u * K * u' = 0,
%     where K is a symmetric matrix. It also returns the K D such that
%     X = XU * D * XU' and the number of iterations needed for convergence.

if ~exist('debug', 'var')
    debug = false;
end
p = inputParser;
addParameter(p, 'kernel', 1);
parse(p, varargin{:});
K = p.Results.kernel;

if ~issymmetric(K)
    error('EK_LYAP:: kernel of the RHS is not symmetric')
end

% Check if the rktoolbox is in the path
if ~exist('rat_krylov', 'file')
    error('rktoolbox not found. Did you forget to add it to the path?');
end

if ~isstruct(A)
    if issparse(A)
        AA = ek_struct(A, issymmetric(A));
        nrmA = normest(A, 1e-2);
    else
        AA = ek_struct(A, false);
        nrmA = norm(A);
    end
    
    AA.nrm = nrmA;
else
    AA = A;
    nrmA = AA.nrm;
end

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, n) r < tol_eps * n;
end

% Start with the initial basis
[VA, KA, HA] = rat_krylov(AA, u, inf);

% Dimension of the space
sa = size(u, 2);

bsa = sa;

Cprojected = ( VA(:,1:size(u,2))' * u ) * K * ( u' * VA(:,1:size(u,2)) );

if ~exist('tol', 'var')
    tol = 1e-8;
end

it = 1;

while sa - 2*bsa < k
    if ( size(VA, 2) + 2 * bsa >= size(VA, 1) )
        
        na = size(VA, 1);
        
        A = AA.multiply(1.0, 0.0, eye(na));
        
        Y = lyap(A, u * K * u');
        
        [QQ, DD] = eig(.5 * (Y + Y'));
        
        rk = sum(arrayfun(@(s) tol(s, max(abs(diag(DD))) / nrmA), ...
            abs(diag(DD))) == false);
        [~,ii] = sort(diag(abs(DD))); ii = ii(end:-1:end-rk+1);
        
        Xu = VA(:,1:size(QQ,1)) * QQ(:,ii) * sqrt(abs(DD(ii,ii)));
        D = diag(sign(diag(DD(ii, ii))));
        
        if debug
            fprintf('%d Dense solver: rank %d, size %d\n', it, rk, max(na,nb));
        end
        
        return;
    end
    [VA, KA, HA] = rat_krylov(AA, VA, KA, HA, [0 inf]);
    
    sa = size(VA, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA(1:end-bsa,:) / KA(1:end-bsa,:);
    Cs = zeros(size(As, 1), size(As, 2));
    
    % FIXME: The above steps can be carried out much more efficiently, just
    % computing the new columns and rows of As and Bs.
    
    Cs(1:size(u,2), 1:size(u,2)) = Cprojected;
    
    [Y, res] = lyap_galerkin(As, Cs, 2*bsa);
    
    %Cs = VA(:,1:size(Y, 1))' * u * K * u' * VA(:,1:size(Y,1));
    % keyboard
    %YY = lyap(As(1:end-2*bsa,1:end-2*bsa), -Cs);
    %keyboard;
    
    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e\n', it, res / norm(Y));
    end
    
    if tol(res, norm(Y)) % res < norm(Y) * tol
        break
    end
    it = it + 1;
end

[QQ, DD] = eig(.5 * (Y+ Y'));

rk = sum(arrayfun(@(s) tol(s, max(abs(diag(DD))) / nrmA), ...
    abs(diag(DD))) == false);
[~,ii] = sort(diag(abs(DD))); ii = ii(end:-1:end-rk+1);

Xu = VA(:,1:size(QQ,1)) * QQ(:,ii) * sqrt(abs(DD(ii,ii)));
D = diag(sign(diag(DD(ii, ii))));

end
