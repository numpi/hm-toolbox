function [hssC] = HSS_HSS_Multiply(hssA,hssB)
%function [hssC] = HSS_HSS_Multiply(hssA,hssB)
%Takes in an hss representation of a matrix A (hssA) and multiplies this with
% the hss representation of another matrix B (hssB).  The output is the hss
% representation of a matrix C (hssC), such that A*B=C.
%
%INPUT:             hssA   (struct) a matrix in HSS form as given by the
%                           function Gen_HSS_Matrix().
%                   hssB   (struct) a matrix in HSS form as given by the
%                           function Gen_HSS_Matrix().
%
%OUTPUT:            hssC    (struct) a matrix in HSS form; a result of the
%                          multiplication of hss1*hss2.
%
% Author: Kristen Lessel, klessel@engineering.ucsb.edu
%
% Copyright (C) 2013 Kristen Lessel, (klessel@engineering.ucsb.edu)
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
tol = hssoption('threshold');
ftree.f = zeros(size(hssA.Rr,2),size(hssB.Wl,2)); %size(hss.Rr,2) = size(hss.Rl,2)

%downsweep to calculate all of the g's
[gtree, ~] = get_gtree(hssA,hssB);

%Main call - contains an upsweep to calculate f's
[hssC, ~] = HSS_HSS_Multiply_iter(hssA, hssB, ftree, gtree);
hssC.topnode = 1;
hssC = hss_compress(hssC, tol);
end

function [gtree,g] = get_gtree(hssA,hssB)

if hssA.leafnode == 1 %if leaf node
    g = hssA.V'*hssB.U; %You must ensure the 2nd dimension of U from hssA and V from hssB match.
    gtree= struct;
else %if not leaf node
    
    [gtreeL, gl] = get_gtree(hssA.A11,hssB.A11);
    [gtreeR, gr] = get_gtree(hssA.A22,hssB.A22);
    
    %Compute new g
    if hssA.topnode == 0
        g = hssA.Wl'*gl*hssB.Rl +hssA.Wr'*gr*hssB.Rr;
    else
        g = [];
    end
    
    gtreeL.g = gl;
    gtreeR.g = gr;
    
    %store g's in tree
    gtree.gtreeL = gtreeL;
    gtree.gtreeR = gtreeR;
end

end

%gtree is hss structure with g's stored at each  node as well
function [hssC, ftree] = HSS_HSS_Multiply_iter(hssA,hssB,ftree,gtree)

hssC = hss();

hssC.ml = hssA.ml;
hssC.mr = hssA.mr;
hssC.nl = hssB.nl;
hssC.nr = hssB.nr;

hssC.topnode = 0;

%if leaf node
if hssA.leafnode == 1
    %Calculate U, V and D at the leaf nodes
    hssC.U = [hssA.U  hssA.D*hssB.U];
    hssC.V = [hssB.D'*hssA.V hssB.V];
    hssC.D = hssA.D*hssB.D + hssA.U*ftree.f*hssB.V';
    hssC.leafnode = 1;
else
    %compute f
    if hssA.topnode == 0
        ftreeL.f = hssA.B12*gtree.gtreeR.g*hssB.B21 + hssA.Rl*ftree.f*hssB.Wl';
        ftreeR.f = hssA.B21*gtree.gtreeL.g*hssB.B12 + hssA.Rr*ftree.f*hssB.Wr';
    else
        ftreeL.f = hssA.B12*gtree.gtreeR.g*hssB.B21;
        ftreeR.f = hssA.B21*gtree.gtreeL.g*hssB.B12;
        
    end
    
    %left call
    [hssL, ftreeL] = HSS_HSS_Multiply_iter(hssA.A11,hssB.A11,ftreeL,gtree.gtreeL);
    
    %right call
    [hssR, ftreeR] = HSS_HSS_Multiply_iter(hssA.A22,hssB.A22,ftreeR,gtree.gtreeR);
    
    %store f's in tree
    ftree.ftreeL = ftreeL;
    ftree.ftreeR = ftreeR;
    
    %store recursive HSS structure
    hssC.A11 = hssL;
    hssC.A22 = hssR;
    
    %Store new B's
    if hssA.topnode == 1
        hssC.B21 = blkdiag(hssA.B21, hssB.B21);
        hssC.B12 = blkdiag(hssA.B12, hssB.B12);
    else
        hssC.B21 = [hssA.B21 hssA.Rr*ftree.f*hssB.Wl'; zeros(size(hssB.B21,1),size(hssA.B21,2)) hssB.B21];
        hssC.B12 = [hssA.B12 hssA.Rl*ftree.f*hssB.Wr'; zeros(size(hssB.B12,1), size(hssA.B12,2)) hssB.B12];
    end
    
    %Store new W's
    if hssA.topnode == 0
        hssC.Wl = [hssA.Wl zeros(size(hssA.Wl,1),size(hssB.Wl,2)); hssB.B21'*gtree.gtreeR.g'*hssA.Wr hssB.Wl];
        hssC.Wr = [hssA.Wr zeros(size(hssA.Wr,1),size(hssB.Wr,2)); hssB.B12'*gtree.gtreeL.g'*hssA.Wl hssB.Wr];
        
        %Store new R's
        hssC.Rl = [hssA.Rl hssA.B12*gtree.gtreeR.g*hssB.Rr; zeros(size(hssB.Rl,1),size(hssA.Rl,2)) hssB.Rl];
        hssC.Rr = [hssA.Rr hssA.B21*gtree.gtreeL.g*hssB.Rl; zeros(size(hssB.Rr,1),size(hssA.Rr,2)) hssB.Rr];
    end
    
    hssC.leafnode = 0;
end

end
