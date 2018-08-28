function [Xu, Xv, As, Bs] = rk_sylv(poles, A, B, u, v, k, tol, debug, nrm_type)
%RK_SYLV Approximate the solution of a Sylvester equation AX + XB' = U*V'.
%
% [XU,XV] = RK_SYLV(POLES, A, B, U, V, K) approximates the solution of the
%     Sylvester equation in the factored form XU * XV'. The variable POLES
%     is a 2 x N matrix containing the poles to use in the rational Krylov
%     method. The poles for A are on the first row, the ones for B on the
%     second one.
%
% [XU, VA] = RK_SYLV(POLES, A, B, U, V, K, TOL, DEBUG) also returns the
%     bases VA and VB, and the optional parameters TOL and DEBUG control
%     the stopping criterion and the debugging during the iteration.
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

if ~exist('nrm_type', 'var')
    nrm_type = 2;
end

% Check if the rktoolbox is in the path
if ~exist('rat_krylov', 'file')
    error('rktoolbox not found. Did you forget to add it to the path?');
end

if ~isstruct(A)
        AA = rk_struct(A);
else
        AA = A;
end

if ~isstruct(B)
    BB = rk_struct(B');
else
    BB = B';
end
nrmA = AA.nrm;
nrmB = BB.nrm;

if size(poles, 1) == 1
	poles = [poles; poles];
end

% Start with the initial basis
% [VA, KA, HA, param_A] = rat_krylov(AA, u, inf);
% [VB, KB, HB, param_B] = rat_krylov(BB, v, inf);

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

% Cprojected = ( VA(:,1:size(u,2))' * u ) * ( v' * VB(:,1:size(v,2)) );

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, nrm) r < tol_eps * nrm;
end

it=1;

% Counter for the vector of poles
counter = 1;

while max(sa-bsa, sb-bsb) < k
    % fprintf('Using poles: (%e, %e)\n', poles(1,counter), poles(2,counter));
    next_inf = counter + find(min(poles(:, counter:end), [], 1) == inf) - 1;
    if isempty(next_inf)
	    next_inf = size(poles, 2);
    else
        next_inf = next_inf(1);
    end
    if ~exist('VA', 'var')
        [VA, KA, HA, param_A] = rat_krylov(AA, u, poles(1, counter:next_inf));
        [VB, KB, HB, param_B] = rat_krylov(BB, v, poles(2, counter:next_inf));
    else
        [VA, KA, HA, param_A] = rat_krylov(AA, VA, KA, HA, poles(1, counter:next_inf), param_A);
        [VB, KB, HB, param_B] = rat_krylov(BB, VB, KB, HB, poles(2, counter:next_inf), param_B);
    end
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    if poles(1, next_inf) == inf && poles(2, next_inf) == inf
        
        % Compute the solution and residual of the projected Lyapunov equation
        As = HA / KA(1:end-bsa,:);
        Bs = HB / KB(1:end-bsb,:);
        Cs = zeros(size(As, 1), size(Bs, 1));
        
        if ~exist('Cprojected', 'var')
            Cprojected = (VA(:,1:bsa)' * u) * (VB(:,1:bsb)'*v)';
        end
        
        Cs(1:size(u,2), 1:size(v,2)) = Cprojected;
        
        [Y, res] = lyap_galerkin(As, Bs, Cs, bsa, bsb);
        
        % You might want to enable this for debugging purposes
        if debug
            fprintf('%d Residue: %e\n', it, res / norm(Y));
        end
        
        if tol(res, norm(Y, nrm_type)) % res < norm(Y) * tol
            break
        end
    end
    
    it = it + 1;
    
    % Switch to the next poles for the next round
    counter = mod(next_inf, length(poles)) + 1;
end

[UU,SS,VV] = svd(Y);

rk = sum(arrayfun(@(s) tol(s, SS(1,1) / max(nrmA, nrmB)), diag(SS)) == false); 

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

