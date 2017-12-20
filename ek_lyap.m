function [Xu, VA, AA] = ek_lyap(A, u, k, tol, debug)
%EK_LYAP Approximate the solution of a AX + XA' = u * u'
%
% XU = EK_LYAP(A, U, K) approximates the solution of the Lyapunov equation
%     in the factored form XU * XU'. 
%
% [XU, VA] = EK_LYAP(A, U, K, TOL, DEBUG) also returns the basis VA, and
%     the optional parameters TOL and DEBUG control the stopping criterion
%     and the debugging during the iteration. 

if ~exist('debug', 'var')
    debug = false;
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

Cprojected = ( VA(:,1:size(u,2))' * u ) * ( u' * VA(:,1:size(u,2)) );
	
if ~exist('tol', 'var')
    tol = 1e-8;
end

it = 1;

while sa - 2*bsa < k
    [VA, KA, HA] = rat_krylov(AA, VA, KA, HA, [0 inf]);
    
    sa = size(VA, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA(1:end-bsa,:) / KA(1:end-bsa,:);
    Cs = zeros(size(As, 1), size(As, 2));
    
    % FIXME: The above steps can be carried out much more efficiently, just
    % computing the new columns and rows of As and Bs. 
    
    Cs(1:size(u,2), 1:size(u,2)) = Cprojected;    
    
    [Y, res] = lyap_galerkin(As, Cs, 2*bsa);

    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e\n', it, res / norm(Y));
    end

    if tol(res, norm(Y)) % res < norm(Y) * tol		
        break
    end        
it = it + 1;   
end

[QQ, DD] = eig(Y);
nn = max(abs(diag(DD)));
ii = find(arrayfun(@(s) tol(s, nn / nrmA), diag(DD)) == false );
%ii = find(arrayfun(@(i) tol(
%diag(DD) > max(diag(DD)) * tol / nrmA);

Xu = VA(:,1:size(QQ,1)) * QQ(:,ii) * sqrt(DD(ii,ii));

end
