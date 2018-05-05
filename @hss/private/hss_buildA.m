function [obj] = hss_buildA(A,tol,kMax,p)
%hss = hss_buildA(A)
% Builds a HSS representation of any matrix A. This method can transform
% every matrix A in to HSS, but this method is very slow.
%
% inputs:
%       1. A - The matrix to be converted to HSS (m x n)
% optional inputs:
%       2. tol - the conversion includes SVD. All singular values below the
%       relative tolerance 'tol' are neglected.
%          (default = 0);
%       3. kMax - the depth of the HSS representation will not exceed kMax.
%          (default = Inf);
%       4. p - probability to accept a next level in the HSS matrix, if
%       all other parameter accept a new level. 0 < p < 1 can lead to very
%       unbalanced trees, which can be suitable for testing.
%          (default = 1);
%
% output:
%       1. hss -  the HSS representation of the matrix
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
% hss_buildA: Build a HSS representation out of any matrix A
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

if nargin<4
    p=1;
    if nargin<3
        kMax = Inf;
        if nargin<2
            tol=0;
        end
    end
end
% call a iterative function to calculate the HSS representation of A
obj = hss_buildA_iter(A,tol,kMax,p);
end

function [obj] = hss_buildA_iter(A,tol,kMax,p,k,m,n,U,V)
% Builds a HSS representation of any matrix A. This iterative method can transform
% every matrix A in to HSS, but this method is very slow.
%
% inputs:
%       1. A - The matrix to be converted to HSS (m x n)
% optional inputs:
%       2. tol - the conversion includes SVD. All singular values below the
%       relative tolerance 'tol' are neglected.
%          (default = 0);
%       3. kMax - the depth of the HSS representation will not exceed kMax.
%          (default = Inf);
%       4. p - probability to accept a next level in the HSS matrix, if
%       all other parameter accept a new level. 0 < p < 1 can lead to very
%       unbalanced trees, which can be suitable for testing.
%          (default = 1);
% internal inputs (do not set them!):
%       5. k - represents the depth of the actual node in the tree
%          (default = 1)
%       6. m,n - size of the actual sub-matrix
%          (default: see in the code)
%       7. U,V - U and V of the HSS representation
%          (default = [])
%
%
% output:
%       1. hss -  the HSS representation of the matrix
%

obj = hss();

if nargin<9
    k=1;
    U=[];
    V=[];
    [m,n] = size(A);
    m = [1:m];
    n = [1:n];
else
    k=k+1;
end

if (length(m)>2 && length(n)>2 && p>rand(1) && k<kMax) % parent node
    % hss = struct('ml',floor((m(end)-m(1)+1)/2),'mr',ceil((m(end)-m(1)+1)/2),'nl',floor((n(end)-n(1)+1)/2),'nr',ceil((n(end)-n(1)+1)/2));
    obj.ml = floor((m(end)-m(1)+1)/2);
    obj.mr = ceil((m(end)-m(1)+1)/2);
    obj.nl = floor((n(end)-n(1)+1)/2);
    obj.nr = ceil((n(end)-n(1)+1)/2);
    obj.leafnode=0;
    obj.topnode=0;

    % indices of the current sub-block      
    %      nl nr
    %   ml
    %   mr
    ml = m(1):m(1)+obj.ml-1;
    mr = m(1)+obj.ml:m(end);
    nl = n(1):n(1)+obj.nl-1;
    nr = n(1)+obj.nl:n(end);

    % indices with out the indices of current sub-block
    butml=allExept(ml,size(A,1));
    butmr=allExept(mr,size(A,1));
    butnl=allExept(nl,size(A,2));
    butnr=allExept(nr,size(A,2));

    % Calculate the SVD of all the Hankel blocks

    [Ul S1 V1] = svdTol(A(ml,butnl),tol);

    [Ur S2 V1] = svdTol(A(mr,butnr),tol);

    [U12 S3 Vl] = svdTol(A(butml,nl),tol);

    [U21 S4 Vr] = svdTol(A(butmr,nr),tol);

    
    % out of the SVD, calculate the Bu and Bl
    Bl = Ur'*U12(m(1):m(1)+size(Ur,1)-1,:)*S3;
    Bu = Ul'*U21(m(1):m(1)+size(Ul,1)-1,:)*S4;
   
    obj.B21=Bl;
    obj.B12=Bu;
   
    %recursion part
    obj.A11 = hss_buildA_iter(A,tol, kMax, p,k,ml,nl,Ul,Vl);
    obj.A22 = hss_buildA_iter(A,tol, kMax, p,k,mr,nr,Ur,Vr);


    if k ~=1  % not top node
        obj.Rl = Ul\U(1:obj.ml,:);
        obj.Rr = Ur\U(obj.ml+1:end,:);

        obj.Wl = Vl\V(1:obj.nl,:);
        obj.Wr = Vr\V(obj.nl+1:end,:);
    else % top node
        obj.topnode=1;
    end


else    % leaf node
    obj = hss(); obj.leafnode = 1; % struct('leafnode', 1);
    obj.topnode=0;
    if k ~=1  % not top node
        obj.U = U;
        obj.V = V;
    else % leaf and top node -> U,V,B not neaded
        obj.topnode=1;
    end
    obj.D = A(m,n);
end

end

function butn = allExept(n,max)
butn=[1:n(1)-1 , n(end)+1 :max];
end

function [U S V] = svdTol(A,tol)
% SVD which neglect all the singular values smaller than the tolerance tol
if nargin<2
    tol=0;
end
    [U S V] = svd(A,'econ');
    i = min(size(S));
    tol=tol*S(1,1);
    i = i-length(find(diag(S(1:i,1:i))<=tol));
    U = U(:,1:i);
    S = S(1:i,1:i);
    V = V(:,1:i);
end
