function [Xu, Xv, VA, VB] = ek_sylv(A, B, u, v, k, tol, debug)
%EK_SYLV Approximate the solution of a Sylvester equation AX + XB' = U*V'.
%
% [XU,XV] = EK_SYLV(A, B, U, V, K) approximates the solution of the 
%     Sylvester equation in the factored form XU * XV'. 
%
% [XU, VA] = EK_SYLV(A, B, U, V, K, TOL, DEBUG) also returns the bases VA 
%     and VB, and the optional parameters TOL and DEBUG control the 
%     stopping criterion and the debugging during the iteration. 
%
% The tolerance TOL can also be specified as a function TOL(R, N) that
% takes as input the residual norm and the norm of the solution (R and N,
% respectively), and returns true if the solution can be accepted. 

if ~exist('debug', 'var')
    debug = false;
end

if ~exist('tol', 'var')
    tol = 1e-8;
end

% Check if the rktoolbox is in the path
if ~exist('rat_krylov', 'file')
    error('rktoolbox not found. Did you forget to add it to the path?');
end

if ~isstruct(A)
	AA = ek_struct(A, issymmetric(A));
else
	AA = A;
end

if ~isstruct(B) 
	BB = ek_struct(B', issymmetric(B));
else
	BB = B';
end

nrmA = AA.nrm;
nrmB = BB.nrm;

% Start with the initial basis
[VA, KA, HA, param_A] = rat_krylov(AA, u, inf);
[VB, KB, HB, param_B] = rat_krylov(BB, v, inf);

param_A.deflation_tol = tol;
param_B.deflation_tol = tol;

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

Cprojected = ( VA(:,1:size(u,2))' * u ) * ( v' * VB(:,1:size(v,2)) );

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, n) r < tol_eps * n;
end

it=1;
while max(sa-2*bsa, sb-2*bsb) < k
    [VA, KA, HA, param_A] = rat_krylov(AA, VA, KA, HA, [0 inf], param_A);
    [VB, KB, HB, param_B] = rat_krylov(BB, VB, KB, HB, [0 inf], param_B);
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA(1:end-bsa,:) / KA(1:end-bsa,:);
    Bs = HB(1:end-bsa,:) / KB(1:end-bsb,:);
    Cs = zeros(size(As, 1), size(Bs, 1));
    
    Cs(1:size(u,2), 1:size(v,2)) = Cprojected;    
    
    [Y, res] = lyap_galerkin(As, Bs, Cs, 2*bsa, 2*bsb);

    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e\n', it, res / norm(Y));
    end

    if tol(res, norm(Y)) % res < norm(Y) * tol
        break
    end        
 it=it+1;
end
 it
% fprintf('lyap its = %d, nrmA = %e\n', it, nrmA)
[UU,SS,VV] = svd(Y);

% rk = sum(diag(SS) > SS(1,1) * tol / max(nrmA, nrmB));
rk = sum(arrayfun(@(s) tol(s, SS(1,1) / max(nrmA, nrmB)), diag(SS)) == false);

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

