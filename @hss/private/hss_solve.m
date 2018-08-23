function [x, time ,n]= hss_solve(hss,b,run,I, J, S)
%x= hss_solve(hss,b)
% Solves the linear system of equation Ax = b if A is in HSS format.
% Theoretical framework:
% S. Chandrasekaran, P. Dewilde, M. Gu, W. Lyons, and T. Pals, A fast
% solver for HSS representations via sparse matrices , SIAM J. Matrix Anal.
% Appl. 29 (2006/07), no. 1, 67--81 (electronic).
% Chapter 4: 'Sparse represenation'
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
%         1. x - the solution of Ax = b.
%         2. time - struct of different execution times measured.
%         3. n - minimal length of I,J,S.
%
% Example
%   m = 10;
%   n = 11;
%   b = rand(m,1);
%   A = rand(m,n);
%   hss = hss_buildA(A);
%   x = hss_solve(hss,b);
%   max(abs(A*x-b))

% Authors:  Stefan Pauli, stefan.pauli@alumni.ethz.ch
%           Karthik Jayaraman Raghuram, jrk@ece.ucsb.edu
% v1.0 Created 21-Okt-09
%
% hss_solve: Solves the linear system of equation Ax = b if A is in HSS format.
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
    
    
    % Compose the sparse matrix S out of the HSS representation
    % tic % time.composeS = cputime;
    [S, indb, indx, n] = hss_toSparse(hss,b,run,I, J, S);
    % time.composeS = toc; %cputime-time.composeS;
    
    % use \ to calculate x
    % tic %time.backslash = cputime;
    bHSS =  sparse(size(indb,1),size(b,2));
    bHSS(indb,:)=b;
    x=S\bHSS;
    % time.backslash = toc; %cputime-time.backslash;
    
    
    % Select the x
    x = x(indx,:);
    
else %run == false
    % tic
    [I, J, S, n] = hss_toSparse(hss,b,run);
    % time.composeS = toc;
    x=[];
end
end
