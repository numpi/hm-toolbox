function [Xu, VA, D, it] = ek_lyap(A, u, maxdim, tol, debug, varargin)
%EK_LYAP Approximate the solution of a AX + XA' + u * u' = 0
%
% [XU, VA, D, it] = EK_LYAP(A, U, MAXDIM, TOL, DEBUG, 'kernel', K) solves
%
%        A X + X A' + u * K * u' = 0,
%
%     where K is a symmetric matrix. It also returns the basis VA, the D such that
%     X = XU * D * XU' and the number of iterations needed for convergence. maxdim
%     is the maximum dimension allowed for the projection space.
%     The optional parameters TOL and DEBUG control the stopping criterion
%     and the debugging during the iteration.

if ~exist('debug', 'var')
    debug = false;
end

p = inputParser;

addOptional(p, 'nrmtype', 2, ...
    @(x) (isnumeric(x) && x == 2) || strcmp(x, 'fro'));

addParameter(p, 'kernel', 1);

parse(p, varargin{:});

K = p.Results.kernel;
nrmtype = p.Results.nrmtype;

if ~issymmetric(K)
    error('EK_LYAP:: kernel of the RHS is not symmetric')
end

if ~isstruct(A)
    if issparse(A)
        AA = ek_struct(A, issymmetric(A));
        nrmA = normest(A, 1e-2);
    else
        AA = ek_struct(A, false);
        nrmA = norm(A, nrmtype);
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

% Dimension of the space
sa = size(u, 2);

bsa = sa;

if ~exist('tol', 'var')
    tol = 1e-8;
end

it = 1;

while sa - 2*bsa < maxdim
    if exist('VA', 'var') && ( size(VA, 2) + 2 * bsa >= size(VA, 1) )
        
        na = size(VA, 1);
        
        A = AA.multiply(1.0, 0.0, eye(na));
        
        Y = lyap(A, u * K * u');
        
        [QQ, DD] = eig(.5 * (Y + Y'));
        
        switch nrmtype
    		case 2
        		s = abs(diag(DD));
        		rk = sum( arrayfun(@(ss) tol(ss, max(s) / (2 *nrmA)), s) == 0);
    		case 'fro'
        		d = sort(abs(diag(DD)));
        		s = sqrt(cumsum(d.^2));
        		rk = sum( arrayfun(@(ss) tol(ss, s(end) / (2 * nrmA)), s) == 0 );
	end
        [~,ii] = sort(diag(abs(DD))); ii = ii(end:-1:end-rk+1);
        
        Xu = VA(:,1:size(QQ,1)) * QQ(:,ii) * sqrt(abs(DD(ii,ii)));
        D = diag(sign(diag(DD(ii, ii))));
        
        if debug
            fprintf('%d Dense solver: rank %d, size %d\n', it, rk, max(na,nb));
        end
        
        return;
    end
    
    if ~exist('VA', 'var')
        [VA, KA, HA, params] = ek_krylov(AA, u);
    else
        [VA, KA, HA, params] = ek_krylov(VA, KA, HA, params);
    end
    
    sa = size(VA, 2);
    
    % Compute the solution and residual of the projected Lyapunov equation
    As = HA / KA(1:end-bsa,:);
    Cs = zeros(size(As, 1), size(As, 1));
    
    if ~exist('Cprojected', 'var')
        Cprojected = ( VA(:,1:size(u,2))' * u ) * K * ( u' * VA(:,1:size(u,2)) );
    end
    
    Cs(1:size(u,2), 1:size(u,2)) = Cprojected;
    
    [Y, res] = lyap_galerkin(As, Cs, bsa, nrmtype);
    
    %Cs = VA(:,1:size(Y, 1))' * u * K * u' * VA(:,1:size(Y,1));
    % keyboard
    %YY = lyap(As(1:end-2*bsa,1:end-2*bsa), -Cs);
    %keyboard;
    
    % You might want to enable this for debugging purposes
    if debug
        fprintf('%d Residue: %e\n', it, res / norm(Y, nrmtype));
    end
    
    if tol(res, norm(Y, nrmtype)) % res < norm(Y) * tol
        break
    end
    it = it + 1;
end

[QQ, DD] = eig(.5 * (Y+ Y'));

switch nrmtype
    	case 2
        	s = abs(diag(DD));
        	rk = sum( arrayfun(@(ss) tol(ss, max(s) / (2 *nrmA)), s) == 0);
    	case 'fro'
        	d = sort(abs(diag(DD)));
        	s = sqrt(cumsum(d.^2));
        	rk = sum( arrayfun(@(ss) tol(ss, s(end) / (2 * nrmA)), s) == 0 );
end

[~,ii] = sort(diag(abs(DD))); ii = ii(end:-1:end-rk+1);

Xu = VA(:,1:size(QQ,1)) * QQ(:,ii) * sqrt(abs(DD(ii,ii)));
D = diag(sign(diag(DD(ii, ii))));

end
