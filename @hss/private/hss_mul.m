function b=hss_mul(hss,x)
%b=hss_mul(hss,x)
% a fast Multiplication of a HSS matrix with a vector x from right or with
% multiple vectors combined to the matrix x
% According to:
% S. Chandrasekaran, P. Dewilde, M. Gu, W. Lyons, and T. Pals, A fast
% solver for HSS representations via sparse matrices , SIAM J. Matrix Anal.
% Appl. 29 (2006/07), no. 1, 67--81 (electronic).
% Chapter 3: 'Fast multiplication'
%
% inputs:
%   hss -  a HSS matrix
%   x - a single vector or several vectors combined to a matrix
% outputs:
%   b - result of the multiplication b = HSS*x
%
% Example
%   m = 10;
%   n = 11;
%   x = rand(n,1);
%   A = rand(m,n);
%   hss = hss_buildA(A);
%   b = hss_mul(hss,x);
%   max(abs(A*x-b))

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% hss_mul: Multiplication of a HSS matrix with a vector from right
% Copyright (C) 2009  Stefan Pauli (stefan.pauli@alumni.ethz.ch)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if (hss.leafnode==0)            %parent node
    vb = upSweep(hss, x);
    f=0;
    b = downSweep(vb, hss, f, x);
else
    b= hss.D*x;
end


end

function [vb, g] = upSweep(hss, x)
% this is the up-sweep recursion from the fast multiplication. eg. equation
% (3) and (4). The resulting G are stored in a Tree.
% comment: G_0_1 is not used, and therefore not calculated
%
% input:
% vb - the tree in which the G are stored
% hss - the tree in which the HSS matrix is stored
% x - the vector for the multiplication (HSS * x)
vb = struct();

if (hss.leafnode==0)   % parent node
    
    [vb.vbl, vb.gl] = upSweep(hss.A11, x(1:hss.nl,:));
    [vb.vbr, vb.gr] = upSweep(hss.A22, x(hss.nl+1:end,:));
    if hss.topnode==0    % not top node
        g = hss.Wl'*vb.gl + hss.Wr'*vb.gr;      % (4)
    end
else        %leav node
    g = hss.V'*x;                           % (3)
end

end

function [b] = downSweep(vb, hss, f, x)
% this is the down-sweep recursion of the fast multiplication. eg. equation
% (5), (6), and (7).
%
% input:
% vb - the tree in which the G are stored
% hss - the tree in which the HSS matrix is stored
% f - the f of eq. (5), (6), and (7).
% x - the vector for the multiplicaton (HSS * x)
%
% output
% b - the result of the multiplicaton (b = HSS * x)

if (hss.leafnode==0)   % parent node
    if hss.topnode==1 % top node
        fl = hss.B12*vb.gr;
        fr = hss.B21*vb.gl;
    else % not top node
        fl = hss.B12*vb.gr+hss.Rl*f;
        fr = hss.B21*vb.gl+hss.Rr*f;
    end
    
    bl = downSweep(vb.vbl, hss.A11, fl, x(1:hss.nl,:));
    br = downSweep(vb.vbr, hss.A22, fr, x(hss.nl+1:end,:));
    b = [bl; br];
    
else        %leaf node
    b = hss.D*x + hss.U*f;
end
end
