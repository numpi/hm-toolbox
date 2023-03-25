function x = hss_chol_solve(A, b)
% HSS_CHOL_SOLVE computes the solution of the linear system A x = b
%              where A is a positive semidefinite HSS matrix and b a block colum.
%	       The procedure is based on the implicit generalized Cholesky factorization
%	       of A and is taken from the paper:
%
%	       Xia, J., Chandrasekaran, S., Gu, M., & Li, X. S. (2010). Fast algorithms for
%              hierarchically semiseparable matrices. Numerical Linear Algebra with Applications, 17(6), 953-976.
%
if A.leafnode == 1
    x = full(A) \ b;
    return
end

if ~issymmetric(A)
    A = hss_symmetrize(A);
end

x = hss_chol_solve_rec(A,b,{},{},{}, [1:size(A,1)]);
end


function [x, D, U, b, ind, cind,  ort, right_ind, triang] = hss_chol_solve_rec(A, b, ort, right_ind, triang, cur_ind)
%------------------------------------------------------------------------------------------------------------------------------------------------
% Recursively compute the solution of the linear system
%
% INPUT  : A coefficient matrix (HSS format), b right hand side, ort, triang and right_ind see description of outputs, cur_ind indices of the current submatrix
%
% OUTPUT : if the function is called on the root then it returns only the solution x
%	   otherwise it returns:
%
%	   - x         contains the computed entries of the solution and 0 elsewhere (the set of variables is the one corresponding to the node)
%	   - D         remaining part of diagonal block
%	   - U         remaining rows of the column generator
%	   - b         part of the right hand side corresponding to the equation that have not be solved yet
%	   - ind       indices of the component of x, corresponding to computed entries up to now
%	   - cind      complementary of the set ind
%	   - ort       array of cells containing the orthogonal transformations on the right
%	   - triang    array of cells containing the triangular transformations on the right
%	   - right_ind   array of cells containing the indices where the right transformations have to be applied
%
%------------------------------------------------------------------------------------------------------------------------------------------------

if A.leafnode == 1  % LEAF
    k = size(A.U, 2);
    [row, col] = size(A.D);
    if k == row
        x = zeros(col, size(b,2)); D = A.D; U = A.U; ind = zeros(1,0); cind = [1:col];
        return
    end
    ind = [1:row-k];  % Indices of computed entries in the solution, at the end
    cind = [row - k + 1: row]; % Complementary of ind
    Q = zeros(row);
    % Compute the Cholesky factor of diagonal block
    LD = chol(A.D,'lower');
    % Compress off-diagonal block row
    [Q(:, end:-1:1), A.U(end:-1:1, :)] = qr(LD\A.U);
    U = A.U(end-k+1:end,:);
    b = Q'*(LD\  b);
    % Retrieve part of the variables
    x = b(1:end-k, :);
    x = [x; zeros(col - size(x,1), size(x,2))];
    % Update transformations on the right
    ort = [ort, Q];
    triang = [triang, LD];
    right_ind = [right_ind, cur_ind];
    % Reduce the diagonal block
    D = eye(k);
    b = b(end-k+1:end, :);
    
else    % NOT A LEAF
    [xl,  Dl, Ul, bl, indl, cindl, ort, right_ind, triang] = hss_chol_solve_rec(A.A11, b(1:A.ml, :), ort, right_ind, triang, cur_ind(1:A.nl));                    % recursive call on left  child
    [xr,  Dr, Ur, br, indr, cindr, ort, right_ind, triang] = hss_chol_solve_rec(A.A22, b(A.ml + 1:A.ml + A.mr, :), ort, right_ind, triang, cur_ind(A.nl+1:end));  % recursive call on right child
    
    x = [xl; xr]; % Store the computed variables
    
    % Merge nodes
    D = [Dl, Ul * A.B12 * Ur'; Ur * A.B21 * Ul', Dr]; % possibile inghippo 2.0
    b = [bl; br];
    
    if A.topnode == 1 % If we are in the root we solve the remaining dense sytem and we apply right transformations
        U = [];
        x([cindl, A.nl + cindr],:) = D \ b;
        % Apply all the transformations on the right
        for j = length(ort):-1:1
            x(right_ind{j}, :) = ort{j} * x(right_ind{j}, :);
            x(right_ind{j}, :) = triang{j}' \ x(right_ind{j}, :);
            
        end
        return
    end
    
    % Otherwise compress the off-diagonal block row in the current level
    U = [Ul * A.Rl; Ur * A.Rr];
    k = size(U,2);
    if k == size(U,1) % full rank, no compression needed
        ind = [indl, A.nl + indr];
        cind = [cindl, A.nl + cindr];
        return
    end
    Q = zeros(size(U,1));
    LD = chol(D,'lower');
    [Q(:, end:-1:1), U(end:-1:1, :)] = qr(LD\U);
    U = U(end-k+1:end,:);
    b = Q'*(LD \ b);
    
    % Retrieve part of the variables
    z = b(1:end-k, :);
    % Update transformations on the right
    ort = [ort, Q];
    triang = [triang, LD];
    right_ind = [right_ind, [cur_ind(cindl), cur_ind(A.nl + cindr)]];
    
    nz = size(z,1);
    n_cindl = length(cindl);
    
    if nz < n_cindl
        x(cindl(1:nz), :) = z;
        %indl = union(indl, cindl(1:nz)); % updates indices of computed entries in the current level
        indl = [indl, cindl(1:nz)];
        cindl = cindl(nz+1:end);
    else
        x(cindl,:) = z(1:n_cindl ,:);
        x(A.nl + cindr(1: nz-n_cindl), :) = z(n_cindl + 1:end, :);
        indl = [1:size(xl,1)];   % updates indices of computed entries in the current level
        indr = [indr, cindr(1: nz-n_cindl)];
        cindl = zeros(1,0);
        cindr = cindr(nz-n_cindl+1:end);
    end
    ind = [indl, A.nl + indr];  % updates indices of computed entries up to the current level
    cind = [cindl, A.nl + cindr]; % Complementary of ind;
    % Reduce the rows in the system
    D = eye(k);
    b = b(end-k+1:end, :);
end
end
