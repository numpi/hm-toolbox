function [Y, res] = care_galerkin(varargin)
%LYAP_GALERKIN Solve the reduced CARE and check the residual.
%
% Y = CARE_GALERKIN(HA, HB, C, bsa, nrmtype) solves the (projected) CARE
%        HA1' * Y + Y * HA1 - Y * HB * HB' * Y + C1 = 0, where HA1 is obtained shrinking
%        HA by BSA columns and rows, 
%        columns and rows, and C1 is obtained cutting C accordingly.
%
% [Y, RES] = CARE_GALERKIN(HA, HB, C, bsa, nrmtype) does the same computation
%        but also computes the residual of the unprojected equation,
%        assuming that HA, HB, and C have been projected using a
%        Krylov-type method and that the action of A on the basis excluding
%        the last BSA columns is contained in the full basis.

if length(varargin) >= 4 && length(varargin) <= 5
    HA = varargin{1};
    HB = varargin{2};
    C = varargin{3};
    bsa = varargin{4};
    if nargin == 5
	nrmtype = varargin{5};
    else
	nrmtype = 2;
    end
    is_lyapunov = true;
else
    error('CARE_GALERKIN:: too many or too few input arguments')
end

% Consider the projected matrices at the previous step, which is needed to
% check the Galerkin condition
HA1 = HA(1 : end - bsa, :);

% Compute the solution of the CARE equation (word of warning: please
% check the sign of C in the implementation of SylvKrylov).
rhs = C(1:end-bsa,1:end-bsa); 
if ~exist('icare', 'file')
	Y = care(HA1', HB, rhs);
else
	%Y = small_care_solve(HA1', HB, rhs, hodlroption('threshold'), 50);
	[Y, converged] = newton_care(HA1', HB * HB', rhs, zeros(size(HA1)), hodlroption('threshold'), 50);
	%Y = icare(HA1', HB, rhs);
	if false && ~converged
		res = inf;
		return
	end
end


if isempty(Y)
	error('CARE_GALERKIN:: There is no finite stabilizing solution for the projected CARE')
end

% Compute the residual
if nrmtype == 2
    res = norm(HA(end-bsa+1:end, :) * Y);
else
    res = sqrt(2) * norm(HA(end-bsa+1:end, :) * Y, 'fro');
end

end

