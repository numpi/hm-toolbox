function [Xu, Xv, VA, VB] = SylvKrylov2(A, B, u, v, k, tol, debug)
%SYLVKRYLOV2 Approximate the solution of a Lyapunov equation. 

% Check if the Sylvester equation is really a Lyapunov equation. 
% FIXME: Notice the sign change, due to the different handling of the
% constant term in LyapKrylov and SylvKrylov.

if ~exist('debug', 'var')
    debug = false;
end

nrmA = normest(A,1e-2);
nrmB = normest(B,1e-2);

if normest(A - B') < eps * max(nrmA, nrmB) && false
    % Determine if A is posdef or negdef
    xx = randn(size(A,1), 1);
    s = sign(xx' * A * xx);
    
	if norm(u + v, 1) < eps * length(u)
		disp('lyap -1')
		[Xu, VA] = LyapKrylov2(A, u, k, tol, nrmA);
		Xv = -s * Xu;
		VB = VA;
		return;
	elseif norm(u - v, 1) < eps * length(u)
        disp('lyap 1')
		[Xu, VA] = LyapKrylov2(A, u, k, tol, nrmA);
		Xv = s * Xu;
		VB = VA;
		return;
	end
end

%disp('sylv')

if ~isstruct(A) 
	[LA, UA, pA, qA] = lu(A);
	AA = struct(... 
	'solve', @(nu, mu ,x) (nu * A - mu * eye(size(A), 'like', A)) \ x, ... 
    'solve_new', @(nu, mu, x) sparse_solve(nu, mu, LA, UA, pA, qA, x, A), ... 
    'multiply', @(rho, eta, x) rho * A * x - eta * x, ...
    'isreal', true);
else
	AA = A;
end

if ~isstruct(B) 
	[LB, UB, pB, qB] = lu(B');
	BB = struct(... 
	'solve', @(nu, mu ,x) (nu * B' - mu * eye(size(B), 'like', B)) \ x, ... 
	'solve_new', @(nu, mu, x) sparse_solve(nu, mu, LB, UB, pB, qB, x, B'), ... 
	'multiply', @(rho, eta, x) rho * B' * x - eta * x, ...
	'isreal', true);
else
	BB = B';
end

% Start with the initial basis
[VA, KA, HA] = rat_krylov(AA, u, inf);
[VB, KB, HB] = rat_krylov(BB, v, inf);

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

Cprojected = ( VA(:,1:size(u,2))' * u ) * ( v' * VB(:,1:size(v,2)) );
	
if ~exist('tol', 'var')
    tol = 1e-8;
end

it=1;
while max(sa-2*bsa, sb-2*bsb) < k
    [VA, KA, HA] = rat_krylov(AA, VA, KA, HA, [0 inf]);
    [VB, KB, HB] = rat_krylov(BB, VB, KB, HB, [0 inf]);
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA(1:end-bsa,:) / KA(1:end-bsa,:);
    Bs = HB(1:end-bsa,:) / KB(1:end-bsb,:);
    Cs = zeros(size(As, 1), size(Bs, 1));
    
    % FIXME: The above steps can be carried out much more efficiently, just
    % computing the new columns and rows of As and Bs. 
    
    Cs(1:size(u,2), 1:size(v,2)) = Cprojected;    
    
    [Y, res] = lyap_galerkin(As, Bs, Cs, 2*bsa, 2*bsb);

    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e\n', it, res / norm(Y));
    end

    if res < norm(Y) * tol
        break
    end        
 it=it+1;   
end

% it
[UU,SS,VV] = svd(Y);

rk = sum(diag(SS) > SS(1,1) * tol / max(nrmA, nrmB));

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

