function [L, U, indb, indx, time ,n]= hss_lu(hss,b,run,I, J, S)
%[L, U, bHSS, indx]= hss_lu(hss,b)
% returns the LU decomposition of the linear system of equation Ax = b if A
% is in HSS format.
%
% inputs:
%         1. hss -  the HSS representation of the matrix A.
%         2. b - the right hand side of Ax = b
% optional inputs:
%         6. run - if true -> run the function
%                  if false -> precompute the memory usage for I,J and S.
%                  n = ( length of I, J and S).
%            (default = true)
%         7. I,J,S - preallocated Vector of the size (nx1).
%            (default: allocated during the process)
% outputs:
%         1. L,U,indb,indx - the LU decomposition of the HSS Matrix. (see help
%         lu) in order to solve the system in linear time, new variables
%         are introduced, and the system gets bigger. to get x make the
%         following:
%                   bHSS =  indb;
%                   bHSS(indb)=b;
%                   x =  U\(L\bHSS);
%                   x = x(indx,:);
%         2. time - struct of different execution times measured.
%         3. n - minimal length of I,J,S.
%
% Example
%   m = 10;
%   n = 11;
%   b = rand(m,1);
%   A = rand(m,n);
%   hss = hss_buildA(A);
%   [L, U, indb, indx]= hss_lu(hss,b);
%   bHSS = zeros(size(indb));
%   bHSS(indb) = b;
%   x = U\(L\bHSS);
%   x = x(indx);
%   max(abs(A*x-b))
%

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% hss_lu: LU decomposition of a HSS matrix
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
    run=true;
end

if run == true
    if nargin < 6
        % Compute the needed size of I, J, and S and allocate them
        [I, J, S, n] = hss_toSparse(hss,b,false);
        I = zeros(n,1);
        J = zeros(n,1);
        S = zeros(n,1);
    end
    
    % prepare the struct for the time measurement
    time=struct();
    
    % Compose the matrix S in the format where every non zero element has
    % two indexes (I,J) and a value (S)
    tic % time.composeS = cputime;
    [S, indb, indx, n] = hss_toSparse(hss,b,run,I, J, S);
    time.composeS = toc; %cputime-time.composeS;
    
    % Do the LU factorization and apply it to b
    tic %time.LU = cputime;
    [L,U] = lu(S);
    time.LU = toc; %cputime-time.LU;
    
else %run == false
    tic
    [I, J, S, n] = hss_toSparse(hss,b,run);
    time.composeS = toc;
    L=[];
    U=[];
    indb=[];
    indx=[];
end
end
