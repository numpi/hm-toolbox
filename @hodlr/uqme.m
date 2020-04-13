function [G,F] = uqme(A2, A1, A0, tol, debug, nrmtype)
% [G,F]=hodlr_cr(A0,A1,A2) solves the UQMEs A0 + A1 X + A2 X^2 = 0 and 
% A2 + A1 Y + A0 Y^2 = 0 by means of cyclic reduction
% A0, A1, A2: matrix coefficients of the UQMEs
% G: solution of the UQME A0 + A1 X + A2 X^2 = 0
% F: solution of the UQME A2 + A1 Y + A0 Y^2 = 0
if ~exist('debug', 'var') || isempty(debug)
	debug = false;
end
if ~exist('tol', 'var') || isempty(tol)
	tol = hodlroption('threshold');
end
if ~exist('nrmtype', 'var') || isempty(nrmtype)
	nrmtype = 2;
end

if nargout == 2	
	[G, F] = hodlr_cr(A2, A1, A0, tol, debug, nrmtype);
elseif nargout == 1
	%G = hodlr_cr(A2, A1, A0, tol, debug, nrmtype);
	G = hodlr_dac_uqme(A2, A1, A0, tol, debug, nrmtype);
else
	error('UQME:: Unvalid number of output arguments')
end

