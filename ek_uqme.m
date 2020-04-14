function [Xu, Xv, VA, VB] = ek_uqme(A, B, X0, u, v, k, tol, debug, nrmtype, sA, sB)
%EK_UQME Approximate the solution the quadratic equation A X^2 + (A X0 + B) X + A X X0 + u v' = 0
%       by means of a rational krylov projection method.
%
% [XU,XV] = EK_UQME(A, B, X0, U, V) approximates the solution of the
%     quadratic equation in the factored form XU * XV'.
%
% [Xu, Xv, VA, VB] = EK_UQME(A, B, X0, U, V, K, TOL, DEBUG, NRMTYPE) also returns the bases VA
%     and VB, and the optional parameters K, TOL, DEBUG and NRMTYPE control the maximum dimension
%     of the Krylov subspace, the threshold for the stopping criterion, the debugging during the iteration
%	  and the norm to which the solution has to be accurate
%           
%

warning('off','MATLAB:nearlySingularMatrix');
if ~exist('debug', 'var')
    debug = false;
end

if ~exist('k', 'var')
    k = 1000;
end

if ~exist('tol', 'var')
    tol = 1e-8;
end
% tol can be function tol(r, n) that is given the residual and the norm, or
% a scalar. In the latter case, we turn it into a function
if isfloat(tol)
    tol_eps = tol;
    tol = @(r, n) r < tol_eps * n;
end

if ~exist('nrmtype', 'var')
    nrmtype = 2;
end

% Poles
tA = .9999; tB = 1.0001;
zA = [tA, -tA];
zB = [tB, -tB];

%-----------Create the structure for (X0 + A^-1 B) and X0'---------------------------
if ~isstruct(A) 
	if  isa(X0, 'hodlr') || isa(X0, 'hss')
		%AA = ek_uqme_struct(A, zA, A * X0 + B);
		AA = ek_uqme_struct(X0 + A\B, zB); % Both structures work; this is one is faster, the former might be more accurate if A is nearly singular
	else
		AA = X0 + A \ B; % Unstructured case
	end
else
    	AA = A;
		nrmA = norm(AA);
end
if ~exist('nrmA', 'var')
	nrmA = AA.nrm;
end

if ~isstruct(X0)
    if(isa(X0,'hodlr') || isa(X0,'hss')) 
		BB = ek_uqme_struct(X0', zB);
		%BB = X0'; 
		%nrmB = norm(BB);
    else
		BB = X0'; % Unstructured case
		nrmB = norm(BB);
    end
else
    BB = X0; % Note that in this case the struct which refers to X0' has to be passed
end
if ~exist('nrmB', 'var')
	nrmB = BB.nrm;
end

%---------------------------------------------------------------------------------------

% multiply the rhs by A^-1
u = A\u;

% Dimension of the space
sa = size(u, 2);
sb = size(v, 2);
n = size(A, 1);
bsa = sa;
bsb = sb;

it = 1;


while max(sa-2*bsa, sb-2*bsb) < k
    %----------------When the Krylov space is too big -----> dense solver-----------
    if exist('VA', 'var') && exist('HA2', 'var') && ( ( size(HA2, 1) + 2 * bsa >= size(VA, 1) ) || ...
            ( size(HB2, 1) + 2 * bsb >= size(VB, 1) ) ) 
        warning('EK_UQME:: Rational Krylov space is equal to the whole space: using dense solver\n');    
		n = size(A,1);
		S = [eye(n), eye(n); zeros(n), full(-X0 - A\B)]; S = S\[full(X0), zeros(n); -u*v', eye(n)];	
        E = S(1:n, 1:n); G = -S(1:n, n+1:2*n);
        P = -S(n+1:2*n, 1:n); F = S(n+1:2*n, n+1:2*n);
        Y = sda(E, F, G, P);
        if debug
            fprintf('Dense solver: It = %d, size = %d, residual = %e\n', it, n, norm(A * Y^2 + (A * X0 + B) * Y + A * Y * X0 + A * u * v'));
        end
        [UU, DD, VV] = svd(Y);        
        switch nrmtype
    		case 2
        		s = diag(DD);
        		rk = sum( arrayfun(@(ss) tol(ss, max(s) / (nrmA + nrmB)), s) == 0);
    		case 'fro'
        		d = sort(diag(DD));
        		s = sqrt(cumsum(d.^2));
        		rk = sum( arrayfun(@(ss) tol(ss, s(end) / (nrmA + nrmB)), s) == 0 );
		end

        
        Xu = UU(:,1:rk) * sqrt(DD(1:rk, 1:rk));
        Xv = VV(:,1:rk) * sqrt(DD(1:rk, 1:rk));

        return;
    end
	%----------------------------------------------------------------------------------
	%---------------Generation/augmentation of Krylov subspaces------------------------
    if ~exist('VA', 'var')	
        [VA, KA, HA, param_A] = rk_krylov(AA, u, zA);
        [VB, KB, HB, param_B] = rk_krylov(BB, v, zB);	
		[~, KA2, HA2, ~] = rk_krylov(AA, VA, KA, HA, inf, param_A);
       	[~, KB2, HB2, ~] = rk_krylov(BB, VB, KB, HB, inf, param_B);	

    else
   		[VA, KA, HA, param_A] = rk_krylov(AA, VA, KA, HA, zA, param_A);
        [VB, KB, HB, param_B] = rk_krylov(BB, VB, KB, HB, zB, param_B);
		[~, KA2, HA2, ~] = rk_krylov(AA, VA, KA, HA, inf, param_A);
        [~, KB2, HB2, ~] = rk_krylov(BB, VB, KB, HB, inf, param_B);
    end    
    sa = size(VA, 2);
    sb = size(VB, 2);
	%----------------------------------------------------------------------------------
    %------Compute the solution and residual of the projected Lyapunov equation--------
    As = HA2(1:end-bsa, :)/KA2(1:end-bsa, :); 
    Bs = HB2(1:end-bsb, :)/KB2(1:end-bsb, :); 
    Cs = zeros(size(As, 1), size(Bs, 1));
    Ar = HA2/KA2(1:end-bsa, :);
    Br = HB2/KB2(1:end-bsb, :);
    if ~exist('Cprojected', 'var')
        Cprojected = ( VA(:,1:size(u,2))' * u ) * ( v' * VB(:,1:size(v,2)) );
    end
    if size(As,1) > size(u,2)
    	Cs(1:size(u,2), 1:size(v,2)) = Cprojected;  
    else
    	Cs = Cprojected(1:size(As, 1), 1:size(Bs, 1));
    end

	if ~exist('Ds', 'var')
		Ds = VB(:,1:size(As,1))'*VA(:,1:size(As,1));
	else
		Ds = [Ds, VB(:, 1:size(Ds,1))' * VA(:, size(Ds,2)+1:size(As,2)); VB(:, size(Ds,1)+1:size(As,1))' * VA(:, 1:size(As,2))];
	end
    % projected solution
    [Y, converged] = newton_nare(Bs', -Ds, Cs, As, zeros(size(As)), debug);
	%----------------------------------------------------------------------------------
%-----------DEBUG----------------------
% YY = VA(:, 1:size(Y,1)) * Y * VB(:,1:size(Y,2))'; % projected back solution
% fprintf('dim. space = %d res. corr = %e res. corr2 = %e\n',size(Y,1), norm((YY)^2+ (X0 + A\B) * YY + YY * X0 + u*v')/norm(YY), norm(A*(YY)^2+ (A*X0 + B) * YY + A * YY * X0 + A*u*v')/norm(YY))
%--------------------------------------

    if converged % if the newton's method did not converge then we enlarge the space moving to next iteration
    	%--------------------Compute and check the residual--------------------------------
    	res = max(norm(Ar(end-bsa+1:end, 1:size(Y,1)) * Y), ...
        	norm(Y * Br(end-bsb+1 : end, 1:size(Y,2))'));
    	if debug
        	fprintf('%d Relative Residual: %e, Absolute residual: %e\n', it, res/norm(Y), res);
    	end
    	if  tol(res, norm(Y, nrmtype))
        	break
    	end
	end
	%----------------------------------------------------------------------------------
    it = it + 1;
end

% Recompression
[UU, DD, VV] = svd(Y);
switch nrmtype
    	case 2
        	s = diag(DD);
        	rk = sum( arrayfun(@(ss) tol(ss, max(s) / (nrmA + nrmB)), s) == 0);
    	case 'fro'
        	d = sort(diag(DD));
        	s = sqrt(cumsum(d.^2));
        	rk = sum( arrayfun(@(ss) tol(ss, s(end) / (nrmA + nrmB)), s) == 0 );
end
Xu = VA(:,1:size(Y,1)) * UU(:,1:rk) * sqrt(DD(1:rk, 1:rk));
Xv = VB(:,1:size(Y,2)) * VV(:,1:rk) * sqrt(DD(1:rk, 1:rk));


end
%----------------------------------------------------------------------
% Auxiliary functions
%----------------------------------------------------------------------

function S = ek_uqme_struct(A, poles, AXB)
% Build a struct for building rational Krylov subspaces for the correction equation of UQMEs.
%
if exist('AXB', 'var')
	if isa(A, 'hodlr')
		[LA, UA] = lu(A);
	elseif isa(A, 'hss')
		ulvA = ulv(A);
	end
	if poles(1) ~= inf
		if isa(A, 'hodlr')
			[LA1, UA1]  = lu(AXB -  poles(1) * A);
		elseif isa(A, 'hss')
			ulvA1  = ulv(AXB -  poles(1) * A);
		end
	else
		LA1 = 0; UA1 = 0;
	end
	if poles(2) ~= inf
		if isa(AXB, 'hodlr')
			[LAm1,UAm1] = lu(AXB -  poles(2) * A);
		elseif isa(AXB, 'hss')
			ulvAm1  = ulv(AXB -  poles(2) * A);
		end
	else
		LA1 = 0; UA1 = 0;
	end
	if isa(AXB, 'hodlr')
		S = struct(...
        		 'solve', @(nu, mu, x) aux(A, x, mu, nu, LA1, UA1, LAm1, UAm1, LA, UA, poles), ...
        		 'multiply', @(rho, eta, x) rho * AXB * x - eta * A * x, ...
        		 'isreal', isreal(AXB), ...
        		 'nrm', normest(AXB, 1e-2));
	elseif isa(AXB, 'hss')
		S = struct(...
        		 'solve', @(nu, mu, x) aux(A, x, mu, nu, ulvA1, [], ulvAm1, [], ulvA, [], poles), ...
        		 'multiply', @(rho, eta, x) rho * AXB * x - eta * A * x, ...
        		 'isreal', isreal(AXB), ...
        		 'nrm', normest(AXB, 1e-2));
	end
else
	if poles(1) ~= inf
		if isa(A, 'hodlr')
			[LA1, UA1]  = lu(A - poles(1) * eye(size(A), 'like', A));
		elseif isa(A, 'hss')
			ULV1 = ulv(A - poles(1) * eye(size(A), 'like', A));
		else
			error('Unsupported structure')
		end
	else
		LA1 = 0; UA1 = 0; poles(1) = poles(2) + 3;
	end
	if poles(2) ~= inf
		if isa(A, 'hodlr')
			[LAm1, UAm1] = lu(A - poles(2) * eye(size(A), 'like', A));
		elseif isa(A, 'hss')
			ULVm1 = ulv(A - poles(2) * eye(size(A), 'like', A));
		else
			error('Unsupported structure')
		end
	else
		LAm1 = 0; UAm1 = 0; poles(2) = poles(1) + 3;
	end
	if isa(A, 'hodlr')
		S = struct(...
            		'solve', @(nu, mu, x) aux2(A, x, mu, nu, LA1, UA1, LAm1, UAm1, poles), ...
            		'multiply', @(rho, eta, x) rho * A * x - eta * x, ...
            		'isreal', isreal(A), ...
            		'nrm', normest(A, 1e-2));
	elseif isa(A, 'hss')
		    S = struct(...
        		... % 'solve', @(nu, mu, x) (nu * A - mu * eye(size(A), 'like', A)) \ x, ...
        		'solve', @(nu, mu, x) aux2(A, x, mu, nu, ULV1, [], ULVm1, [], poles),...
        		'multiply', @(rho, eta, x) rho * (A * x) - eta * x, ...
        		'isreal', isreal(A), ...
        		'nrm', normest(A, 1e-2));
	end
end
end
%----------------------------------------------------------------------
%----------------------------------------------------------------------

function y = aux(A, x, mu, nu, LA1, UA1, LAm1, UAm1, LA, UA, poles)  
if nu == 0
	if isa(A, 'hodlr')
		y = -lu_solve(1.0, 0.0, LA, UA, x) / mu;
		return
	elseif isa(A, 'hss')
		y = -ulv_solve(LA, x) / mu;
		return
	end
end
if mu == poles(1)
	if isa(A, 'hodlr')
		y = lu_solve(1.0, 0, LA1, UA1, x);
	elseif isa(A, 'hss')
		y = ulv_solve(LA1, x);
	end
elseif mu == poles(2)
	if isa(A, 'hodlr')
		y = lu_solve(1.0, 0, LAm1, UAm1, x);
	elseif isa(A, 'hss')
		y = ulv_solve(LAm1, x);		
	end
end
end
%----------------------------------------------------------------------
%----------------------------------------------------------------------

function y = aux2(A, x, mu, nu, LA1, UA1, LAm1, UAm1, poles)  
if nu == 0
	y = -x/mu;
	return
end
if nu == 1/poles(1)
	if isa(A, 'hodlr')
		y = lu_solve(1.0, 0, LA1, UA1, x) * poles(1);
	elseif isa(A, 'hss')
		y = ulv_solve(LA1, x) * poles(1);
	end
elseif nu == 1/poles(2)
	if isa(A, 'hodlr')
		y = lu_solve(1.0, 0, LAm1, UAm1, x) * poles(2);
	elseif isa(A, 'hss')
		y = ulv_solve(LAm1, x) * poles(2);
	end
else
	y = (nu * A - mu * eye(size(A), 'like', A))\x;
end
end
        
