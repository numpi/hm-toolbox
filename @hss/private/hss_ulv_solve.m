function x = hss_ulv_solve(A, b)
% HSS_ULV_SOLV computes the solution of the linear system A x = b
%              where A is a HSS matrix and b a block colum. The
%	       procedure is based on the implicit ULV factorization
%	       of A and is taken from the paper:
%
%	       Chandrasekaran, Shiv, Ming Gu, and Timothy Pals.
%	       "A fast ULV decomposition solver for hierarchically semiseparable representations."
%	       SIAM Journal on Matrix Analysis and Applications 28.3 	(2006): 603-622.
%
if A.leafnode == 1
    x = full(A) \ b;
    return
end

x = hss_ulv_solve_rec(A,b,{},{}, [1:size(A,1)]);
end

%TODO Check whether it is the most efficient implementation, for example it seems that at the moment we re-construct the generators for the columns (the Vs)
%     See if there is some redundant data (e.g.: set of indices) that can be avoided

function [x, D, U, V, b, ind, cind,  ort, ort_ind] = hss_ulv_solve_rec(A, b, ort, ort_ind, cur_ind)
%------------------------------------------------------------------------------------------------------------------------------------------------
% Recursively compute the solution of the linear system
%
% INPUT  : A coefficient matrix (HSS format), b right hand side, ort and ort_ind see description of outputs, cur_ind indices of the current submatrix
%
% OUTPUT : if the function is called on the root then it returns only the solution x
%	   otherwise it returns:
%
%	   - x         contains the computed entries of the solution and 0 elsewhere (the set of variables is the one corresponding to the node)
%	   - D         remaining compressed diagonal block
%	   - U         remaining rows of the column generator
%	   - V         updated column generator
%	   - b         part of the right hand side corresponding to the equation that have not be solved yet
%	   - ind       indices of the component of x, corresponding to computed entries up to now
%	   - cind      complementary of the set ind
%	   - ort       array of cells containing the orthogonal transformations on the right
%	   - ort_ind   array of cells containing the indices where the right orthogonal transformations have to be applied
%
%------------------------------------------------------------------------------------------------------------------------------------------------

if A.leafnode == 1  % LEAF
    k = size(A.U, 2);
    [row, col] = size(A.D);
    if k == row
        x = zeros(col, size(b,2)); D = A.D; U = A.U; V = A.V; ind = zeros(1,0); cind = [1:col];
        return
    end
    ind = [1:row-k];  % Indices of computed entries in the solution, at the end
    cind = [row - k + 1: row]; % Complementary of ind
    Q = zeros(row);
    % Compress off-diagonal block row
    [Q(:, end:-1:1), A.U(end:-1:1, :)] = qr(A.U);
    U = A.U(end-k+1:end,:);
    A.D = Q' * A.D;
    b = Q' * b;
    % Triangularize part of diagonal block
    [Q, L] = qr(A.D(1:end-k,:)');
    % Retrieve part of the variables
    x = (L(1:end-k,:)')\b(1:end-k, :);
    x = [x; zeros(col - size(x,1), size(x,2))];
    % Update transformations on the right
    V = Q' * A.V;
    ort = [ort, Q];
    ort_ind = [ort_ind, cur_ind];
    % Reduce the rows in the system
    D = A.D(end-k+1:end, :) * Q;
    b = b(end-k+1:end, :) - D(:, 1:end-k) * x(1:end-k,:);
    % Reduce the columns of the diagonal block
    D = D(:,end-k+1:end);
    
else    % NOT A LEAF
    [xl,  Dl, Ul, Vl, bl, indl, cindl, ort, ort_ind] = hss_ulv_solve_rec(A.A11, b(1:A.ml, :), ort, ort_ind, cur_ind(1:A.nl));                    % recursive call on left  child
    
    [xr,  Dr, Ur, Vr, br, indr, cindr, ort, ort_ind] = hss_ulv_solve_rec(A.A22, b(A.ml + 1:A.ml + A.mr, :), ort, ort_ind, cur_ind(A.nl+1:end));  % recursive call on right child
    
    x = [xl; xr]; % Store the computed variables
    
    % Reduce the columns in the system
    bl = bl - Ul * A.B12 * (Vr(indr,:)'* xr(indr, :));
    br = br - Ur * A.B21 * (Vl(indl,:)'* xl(indl, :));
    % Merge nodes
    D = [Dl, Ul * A.B12 * Vr(cindr, :)'; Ur * A.B21 * Vl(cindl, :)', Dr];
    b = [bl; br];
    
    if A.topnode == 1 % If we are in the root we solve the remaining dense sytem and we apply right transformations
        U = []; V = [];
        x([cindl, A.nl + cindr],:) = D \ b;
        % Apply all the transformations on the right
        for j = length(ort):-1:1
            x(ort_ind{j}, :) = ort{j} * x(ort_ind{j}, :);
        end
        return
    end
    
    % Otherwise compress the off-diagonal block row in the current level
    V = [Vl * A.Wl; Vr * A.Wr];
    U = [Ul * A.Rl; Ur * A.Rr];
    k = size(U,2);
    if k >= size(U,1) % full rank, no compression needed
        ind = [indl, A.nl + indr];
        cind = [cindl, A.nl + cindr];
        return
    end
    Q = zeros(size(U,1));
    [Q(:, end:-1:1), U(end:-1:1, :)] = qr(U);
    U = U(end-k+1:end,:);
    
    D = Q' * D;
    b = Q' * b;
    % Triangularize (1,1) submatrix of diagonal block
    [Q, L] = qr(D(1:end-k,:)');
    
    % Retrieve part of the variables
    z = (L(1:end-k,:)')\b(1:end-k, :);
    % Update transformations on the right
    V([cindl, A.nl + cindr],:) = Q' * V([cindl, A.nl + cindr],:);
    ort = [ort, Q];
    ort_ind = [ort_ind, [cur_ind(cindl), cur_ind(A.nl + cindr)]];
    
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
        %indr = union(indr, cindr(1: nz-n_cindl));
        indr = [indr, cindr(1: nz-n_cindl)];
        cindl = zeros(1,0);
        cindr = cindr(nz-n_cindl+1:end);
    end
    ind = [indl, A.nl + indr];  % updates indices of computed entries up to the current level
    cind = [cindl, A.nl + cindr]; % Complementary of ind;
    % Reduce the rows in the system
    D = D(end-k+1:end, :) * Q;
    b = b(end-k+1:end, :) - D(:, 1:end-k) * z;
    % Reduce the columns of the diagonal block
    D = D(:,end-k+1:end);
end
end
