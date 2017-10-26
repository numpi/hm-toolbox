function [Y, res] = lyap_galerkin(HA, HB, C, bsa, bsb)
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

% Consider the projected matrices at the previous step, which is needed to
% check the Galerking condition
HA1 = HA(1 : end - bsa, 1 : end - bsa);
HB1 = HB(1 : end - bsb, 1 : end - bsb);

% Compute the solution of the Lyapunov equation (word of warning: please
% check the sign of C in the implementation of SylvKrylov).
if norm(HA1 - HB1', 1) > 0
	Y = lyap(HA1, HB1', -C(1:end-bsa,1:end-bsb));
else
	Y = lyap(HA1, -C(1:end-bsa,1:end-bsb));
end

% Check the residual
res = max(norm(HA(end-bsa+1:end, 1 : end - bsa) * Y), ...
          norm(Y * HB(end-bsb+1 : end, 1 : end - bsb)'));

end

