function [Y, res] = lyap_galerkin(varargin)
%LYAP_GALERKIN Solve the reduced Lyapunov equation and check the residual.
%
% Y = LYAP_GALERKIN(HA, HB, C, bsa, bsb) solves the (projected) Lyapunov
%        equation HA1 * Y + Y * HB1' = C1, where HA1 is obtained shrinking
%        HA by BSA columns and rows, HB1 is obtained shrinking HB1 by BSB
%        columns and rows, and C1 is obtained cutting C accordingly.
%
% [Y, RES] = LYAP_GALERKIN(HA, HB, C, bsa, bsb) does the same computation
%        but also computes the residual of the unprojected equation,
%        assuming that HA, HB, and C have been projected using a
%        Krylov-type method and that the action of A on the basis excluding
%        the last BSA (resp. BSB) columns is contained in the full basis.

if length(varargin) >= 3 && length(varargin) <= 4
    HA = varargin{1};
    C = varargin{2};
    bsa = varargin{3};
    if nargin == 4
	nrmtype = varargin{4};
    else
	nrmtype = 2;
    end
    is_lyapunov = true;
else
    HA = varargin{1};
    HB = varargin{2};
    C = varargin{3};
    bsa = varargin{4};
    bsb = varargin{5};
    if nargin == 6
	nrmtype = varargin{6};
    else
	nrmtype = 2;
    end
    is_lyapunov = false;
end

% Consider the projected matrices at the previous step, which is needed to
% check the Galerkin condition
HA1 = HA(1 : end - bsa, :);

if ~is_lyapunov
    HB1 = HB(1 : end - bsb, :);
end

% Compute the solution of the Lyapunov equation (word of warning: please
% check the sign of C in the implementation of SylvKrylov).
if ~is_lyapunov
    Y = lyap(HA1, HB1', C(1:end-bsa,1:end-bsb));
else
    Y = lyap(HA1, C(1:end-bsa,1:end-bsa));
end

% Check the residual
if ~is_lyapunov
    if nrmtype == 2
    	res = max(norm(HA(end-bsa+1:end, :) * Y), ...
        norm(Y * HB(end-bsb+1 : end, :)'));
    else
	res = norm( [HA(end-bsa+1:end, :) * Y; ...
        HB(end-bsb+1 : end, :) * Y'], 'fro');
    end
else
    if nrmtype == 2
    	res = norm(HA(end-bsa+1:end, :) * Y);
    else
	res = sqrt(2) * norm(HA(end-bsa+1:end, :) * Y, 'fro');
    end
end

end

