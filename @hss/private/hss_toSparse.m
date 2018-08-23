function [S, indb, indx, n] = hss_toSparse (hss,b,run,I, J, S)
%[S, bHSS, indx] = hss_toSparse (hss,b)
% Build a sparse matrix in nested dissection ordering out of the HSS
% representation. This can be used to solve the equation HSS*x=b for x in
% linear time.
% Theoretical framework:
% S. Chandrasekaran, P. Dewilde, M. Gu, W. Lyons, and T. Pals, A fast
% solver for HSS representations via sparse matrices , SIAM J. Matrix Anal.
% Appl. 29 (2006/07), no. 1, 67--81 (electronic).
% Chapter 4: 'Sparse represenation'
%
% inputs:
%       1. hss -  the HSS representation of the matrix.
%       2. b - right hand side of the equation.
% optional inputs:
%       3. run - if true -> run the function
%                if false -> precompute the memory usage for I,J and S.
%                n = ( length of I, J and S).
%          (default = true)
%       4. I,J,S - preallocated Vector of the size (nx1).
%          (default: allocated during the process)
%
% output:
%       1. S - the sparse matrix in nested dissection ordering of the HSS
%       matrix.
%       2. indb, indx - The matrix S is bigger than the original matrix A.
%       indb and indx are used to place the entries of b to the right
%       places in the new bigger b and to extract x out of the bigger
%       solution vector.
%       3. n - minimal length of I,J,S.
%
% Example
%   m = 10;
%   n = 11;
%   b = rand(m,1);
%   A = rand(m,n);
%   hss = hss_buildA(A);
%   [S, indb, indx] = hss_toSparse(hss,b);
%   bHSS =  sparse(size(indb,1),size(b,2));
%   bHSS(indb,:)=b;
%   x=S\bHSS;
%   x = x(indx,:);
%   max(abs(A*x-b))
%

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% hss_toSparse: Build a sparse matrix the HSS representation.
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

% Default case for run: run = true
if nargin < 3
    run = true;
end

if nargin < 6
    if run==true
        % Compute the needed size of I, J, and S and allocate them
        [S, indb, indx, n] = hss_toSparse(hss,b,false);
        I = zeros(n,1);
        J = zeros(n,1);
        S = zeros(n,1);
    else
        I=[];
        J=[];
        S=[];
    end
end


[I, J, S, n, indb, indx] = hss_toSparse_iter (hss,true(size(b,1),1),run,I, J, S);
if run==true
    % make a standard matlab sparse matrix
    S = sparse(I(1:n), J(1:n), S(1:n));
end
end


function [I, J, S, n, indb, indx, x0, y0, xB0, yB0] = hss_toSparse_iter (hss,b,run,I, J, S,n,indb,x0,y0)
% Build a sparse matrix in nested dissection ordering out of the HSS
% representation. This can be used to solve the equation HSS*x=b for x in
% linear time.
%
% inputs:
%       1. hss -  the HSS representation of the matrix.
%       2. b - vector from one to the length(rigth hand side of the equation).
%       3. run - if true -> run the function
%                if false -> precompute the memory usage for I,J and S.
%                n = ( length of I, J and S).
%          (default = true)
%       4. I,J,S - preallocated Vector of the size (nx1).
%          (default: allocated during the process)
%       5. n - 1 plus # elements stored in S.
% optional inputs:
%       6. indb, indx - The matrix S is bigger than the original matrix A.
%       indb and indx are used to place the entries of b to the right
%       places in the new bigger b and to extract x out of the bigger
%       solution vector.
%       7. x0,y0 - starting from this coordinate the node can writhe his
%       matrices D,U V' and two I in tho the sparse matrix S.
% output:
%       1. I,J,S -  the sparse matrix in a form where every nonzero has
%       its value S and the coordinats I,J.
%       2. # elements stored in S plus 1.
%       3. indb, indx - The matris S is bigger than the original matrix A.
%       indb and indx are used to place the entries of b to the right
%       places in the new bigger b and to extract x out of the bigger
%       result vector.
%       4. x0,y0 - starting from this coordinate the node can writhe his
%       matrices D,U V' and two I in tho the sparse matrix S.
%       5. xB0,yB0 - Index where the parent node can write W, B in tho the
%       sparse matrix S.

% set some values, just in the top node of the tree.
if nargin < 10
    indb = sparse(logical([]));
    x0=1;
    y0=1;
    n=0;
end

if hss.topnode==true
    n=0;
end


if hss.leafnode==false
    [I, J, S, n, indb, indxl, x0, y0, xBl0, yBu0] = hss_toSparse_iter(hss.A11,b(1:hss.ml),run,I, J, S, n,indb,x0,y0);
    [I, J, S, n, indb, indxr, x0, y0, xBu0, yBl0] = hss_toSparse_iter(hss.A22,b(hss.ml+1:end),run,I, J, S, n,indb,x0,y0);
    
    indx = [indxl; indxr];
    
    y1 = size(hss.B12,1);
    y2 = size(hss.B21,1);
    x1 = size(hss.B21,2);
    x2 = size(hss.B12,2);
    
    %Bu
    %oldS(yBu0:yBu0+y1-1 , xBu0:xBu0+x2-1) = hss.B12;
    [I, J, S, n] = addSparce(I, J, S, n, yBu0:yBu0+y1-1, xBu0:xBu0+x2-1, hss.B12,run);
    %Bl
    %oldS(yBl0:yBl0+y2-1 , xBl0:xBl0+x1-1) = hss.B21;
    [I, J, S, n] = addSparce(I, J, S, n, yBl0:yBl0+y2-1 , xBl0:xBl0+x1-1, hss.B21,run);
    
    
    if hss.topnode==false
        y3 = size(hss.Wl',1);
        x4 = size(hss.Rl,2);
        x3 = y3;
        y4 = x4;
        
        indx = [indx; sparse(false(x3+x4,1))];
        
        %Rl
        %oldS(yBu0:yBu0+y1-1 , x0+x3:x0+x3+x4-1) = hss.Rl;
        [I, J, S, n] = addSparce(I, J, S, n, yBu0:yBu0+y1-1 , x0+x3:x0+x3+x4-1, hss.Rl,run);
        %Rr
        %oldS(yBl0:yBl0+y2-1 , x0+x3:x0+x3+x4-1) = hss.Rr;
        [I, J, S, n] = addSparce(I, J, S, n, yBl0:yBl0+y2-1 , x0+x3:x0+x3+x4-1, hss.Rr,run);
        %Wl
        %oldS(y0:y0+y3-1 , xBl0:xBl0+x1-1) = hss.Wl';
        [I, J, S, n] = addSparce(I, J, S, n, y0:y0+y3-1 , xBl0:xBl0+x1-1, hss.Wl',run);
        %Wr
        %oldS(y0:y0+y3-1 , xBu0:xBu0+x2-1) = hss.Wr';
        [I, J, S, n] = addSparce(I, J, S, n, y0:y0+y3-1 , xBu0:xBu0+x2-1, hss.Wr',run);
        %-I
        %oldS(y0:y0+y3-1 , x0:x0+x3-1) = -spdiags(ones(y3,1), 0, y3,y3);
        [I, J, S, n] = addSparceDiag(I, J, S, n, y0:y0+y3-1 , x0:x0+x3-1, -1,run);
        %-I
        %oldS(y0+y3:y0+y3+y4-1 , x0+x3:x0+x3+x4-1) = -spdiags(ones(y4,1), 0, y4,y4);
        [I, J, S, n] = addSparceDiag(I, J, S, n, y0+y3:y0+y3+y4-1 , x0+x3:x0+x3+x4-1, -1,run);
        
        indb(y0:y0+y3+y4-1) = sparse(false(y3+y4,size(indb,2)));
        
        yB0 = y0+y3;
        xB0 = x0;
        y0 = y0+y3+y4;
        x0 = x0+x3+x4;
        
    else
        yB0=[];
        xB0=[];
    end
    
    
else %leaf node
    
    y1 = size(hss.D,1);
    y2 = size(hss.V',1);
    y3 = size(hss.U,2);
    x1 = size(hss.D,2);
    x2 = y2;
    x3 = y3;
    
    yB0 = y0+y1+y2;
    xB0 = x0+x1;
    
    %D
    %oldS(y0:y0+y1-1, x0:x0+x1-1)=hss.D;
    [I, J, S, n] = addSparce(I, J, S, n, y0:y0+y1-1, x0:x0+x1-1, hss.D,run);
    %V'
    %oldS(y0+y1:y0+y1+y2-1 , x0:x0+x1-1)=hss.V';
    [I, J, S, n] = addSparce(I, J, S, n, y0+y1:y0+y1+y2-1 , x0:x0+x1-1, hss.V',run);
    %-I
    %oldS(y0+y1:y0+y1+y2-1 , x0+x1:x0+x1+x2-1) = -spdiags(ones(y2,1), 0, y2,y2);
    [I, J, S, n] = addSparceDiag(I, J, S, n, y0+y1:y0+y1+y2-1 , x0+x1:x0+x1+x2-1, -1,run);
    %U
    %oldS(y0:y0+y1-1 , x0+x1+y2:x0+x1+y2+x3-1) = hss.U;
    [I, J, S, n] = addSparce(I, J, S, n, y0:y0+y1-1 , x0+x1+y2:x0+x1+y2+x3-1, hss.U,run);
    %-I
    %oldS(y0+y1+y2:y0+y1+y2+y3-1 , x0+x1+y2:x0+x1+y2+x3-1) = -spdiags(ones(x3,1), 0, x3,x3);
    [I, J, S, n] = addSparceDiag(I, J, S, n, y0+y1+y2:y0+y1+y2+y3-1 , x0+x1+y2:x0+x1+y2+x3-1, -1,run);
    
    indb(y0:y0+y1-1,:) = b;
    indb(y0+y1:y0+y1+y2+y3-1,:) = sparse(false(x2+x3,size(indb,2)));
    
    indx = sparse([true(x1,1); false(x2+x3,1)]);
    
    y0 = y0+y1+y2+y3;
    x0 = x0+x1+x2+x3;
end
end

function [I, J, S, n] = addSparce(I, J, S, n, i,j,s,run)
% Add the matrix s in the position i,j to the sparse Matrix. The sparse
% Matrix has already n nonzero elements.
% start with: [I, J, S, n] = addSparce([], [], [], 0, i,j,s)
%
% inputs:
%       1. I,J,S - preallocated Vector of the size (nx1).
%          (default: allocated during the process)
%       2. n - 1 plus # elements stored in S.
%       3. i,j,s - store the matrix s in S(i,j), where i and j can be a
%       vector.
%       4. run - if true -> run the function
%                if false -> precompute the memory usage for I,J and S.
%                n = ( length of I, J and S).
% output:
%       1. I,J,S -  the sparse matrix in a form where every nonzero has
%       its value S and the coordinates I,J.
%       2. n - 1 plus # elements stored in S.

if run==true
    i=i(:);
    j=j(:);
    num = length(i)*length(j);
    if num+n >length(I)
        warning('slow down by dynamically allocated memory for I,J and S');
    end
    
    I(n+1:n+num) =repmat(i,length(j),1);
    temp = repmat(j,1,length(i))';
    J(n+1:n+num) = temp(:);
    S(n+1:n+num) = s(:);
    n=n+num;
else
    num = length(i)*length(j);
    n = n+num;
end

end

function [I, J, S, n] = addSparceDiag(I, J, S, n, i,j,s,run)
% Add s * identity matrix in the position i,j to the sparse Matrix. the sparse
% Matrix has already n nonzero elements.
%
% inputs:
%       1. I,J,S - preallocated Vector of the size (nx1).
%          (default: allocated during the process)
%       2. n - 1 plus # elements stored in S.
%       3. i,j,s - store the matrix s in S(i,j), where i and j can be a
%       vector.
%       4. run - if true -> run the function
%                if false -> precompute the memory usage for I,J and S.
%                n = ( length of I, J and S).
% output:
%       1. I,J,S -  the sparse matrix in a form where every nonzero has
%       its value S and the coordinates I,J.
%       2. n - 1 plus # elements stored in S.
if run==true
    i=i(:);
    j=j(:);
    num = length(j);
    if num+n >length(I)
        warning('slow down by dynamically allocated memory for I,J and S');
    end
    I(n+1:n+num) = i;
    J(n+1:n+num) = j;
    S(n+1:n+num) = ones(num,1)*s;
    n=n+num;
    
else
    num = length(i)*length(j);
    n = n+num;
end

end
