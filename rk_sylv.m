function [Xu, Xv, VA, VB] = rk_sylv(poles, A, B, u, v, k, tol, debug, nrmtype)
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

if ~exist('nrmtype', 'var')
    nrmtype = 2;
end
if size(u, 2) ~= size(v, 2)
	error('RK_SYLV:: different number of columns in the factorization of the rhs');
end
if size(u, 2) == 0
    Xu = u;
    Xv = v;
    VA = u;
    VB = v;
    return;
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

if isfloat(poles) && size(poles, 1) == 1
    poles = [poles; poles];
end

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);

bsa = sa;
bsb = sb;

% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, nrm) r < tol_eps * nrm;
end

it=1;

% Counter for the vector of poles
counter = 1;

residuals = [];

% rk_krylov = @rat_krylov;


[VA, KA, HA, param_A] = rk_krylov(AA, u, inf);
[VB, KB, HB, param_B] = rk_krylov(BB, v, inf);
Cprojected = (VA(:,1:bsa)' * u) * (VB(:,1:bsb)'*v)';

while max(sa-bsa, sb-bsb) < k
    
    poleA = poles(1, counter);
    poleB = poles(2, counter);    
    
    [VA, KA, HA, param_A] = rk_krylov(AA, VA, KA, HA, poleA, param_A);
    [VB, KB, HB, param_B] = rk_krylov(BB, VB, KB, HB, poleB, param_B);
    
    sa = size(VA, 2);
    sb = size(VB, 2);
    
    % We only compute the solution and residual of the projected Lyapunov 
    % equation in the first three iterations, and then use the fact that
    % the convergence is expected to be linear to estimate the number of
    % iterations requires to converge. This saves some Lyapunov dense
    % solutions, which are relatively expensive.
    if  it <= 3 || ~exist('tol_eps', 'var')
        check_residual = true;
    else
        if it == 4
            pp = polyfit(1:3, log(residuals(1:3)), 1);
            r1 = residuals(1) / norm(Y, nrmtype);
            needed_iterations = ceil(...
                (log(tol_eps) - log(r1)) / pp(1));
            needed_iterations = min(20, needed_iterations);
        end
        check_residual = it >= needed_iterations;
    end
    
    check_residual = true;
    
    if check_residual
        if poleA ~= inf
            [~, KA2, HA2] = rk_krylov(AA, VA, KA, HA, inf, param_A);
        else
            KA2 = KA; HA2 = HA;
        end
        
        if poleB ~= inf        
            [~, KB2, HB2] = rk_krylov(BB, VB, KB, HB, inf, param_B);
        else
            KB2 = KB; HB2 = HB;
        end
        
        % Compute the solution and residual of the projected Lyapunov equation
        %fprintf('cond KA = %e, cond KB = %e\n', cond(KA2(1:end-bsa,:)), cond(KB2(1:end-bsb,:)))
        As = HA2 / KA2(1:end-bsa,:);
        Bs = HB2 / KB2(1:end-bsb,:);
        Cs = zeros(size(As, 1), size(Bs, 1));
        Cs(1:size(u,2), 1:size(v,2)) = Cprojected;
        
        [Y, res] = lyap_galerkin(As, Bs, Cs, bsa, bsb, nrmtype);
        temp=KA2(1:end-bsa, :)\Y;
        keyboard
        fprintf('res = %e\n', norm(temp(end-bsa+1:end, :)))
        residuals(it) = res;
        
        % You might want to enable this for debugging purposes
        if debug
            fprintf('%d Residue: %e\n', size(Y,1), res / norm(Y, nrmtype));
        end
        
        if tol(res, norm(Y, nrmtype)) % res < norm(Y) * tol
            break
        end
    end
    
    it = it + 1;
    
    % Switch to the next poles for the next round
    if isfloat(poles)
        counter = mod(counter, size(poles, 2)) + 1;
    else
        counter = counter + 1;
    end
end

[UU,SS,VV] = svd(Y);

switch nrmtype
    case 2
        s = diag(SS);
        rk = sum( arrayfun(@(ss) tol(ss, s(1) / (nrmA + nrmB)), s) == 0);
    case 'fro'
        d = sort(diag(SS));
        s = sqrt(cumsum(d.^2));
        rk = sum( arrayfun(@(ss) tol(ss, s(end) / (nrmA + nrmB)), s) == 0 );
end

Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(SS(1:rk,1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(SS(1:rk,1:rk));

end

