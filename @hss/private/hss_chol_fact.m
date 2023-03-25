function F = hss_chol_fact(A)
% HSS_ULV_SOLV computes the  generalized Cholesky factorization
%	       of A  taken from the paper:
%
%	       Xia, J., Chandrasekaran, S., Gu, M., & Li, X. S. (2010). Fast algorithodlrs for
%              hierarchically semiseparable matrices. Numerical Linear Algebra with Applications, 17(6), 953-976.
%


F = struct();
F.U = {};
F.ind = {};
F.L = {};

if A.leafnode == 1
    F.L = {chol(full(A))};
    F.ind = {[1:size(A,1)]};
    return
end

F = hss_chol_fact_rec(A, F, [1:size(A,1)]);
end



function [F, D, U, ind, cind] = hss_chol_fact_rec(A, F, cur_ind)
%------------------------------------------------------------------------------------------------------------------------------------------------
%
% INPUT  : A coefficient matrix (HSS format), ort, triang and right_ind see description of outputs, cur_ind indices of the current submatrix
%
% OUTPUT : if the function is called on the root then it returns only the solution x
%	   otherwise it returns:
%
%	   - D         remaining part of diagonal block
%	   - U         remaining rows of the column generator
%	   - ind       indices of the component of x, corresponding to computed entries up to now
%	   - cind      complementary of the set ind
%	   - F         structure containing transformation and indices
%
%------------------------------------------------------------------------------------------------------------------------------------------------

if A.leafnode == 1  % LEAF
    k = size(A.U, 2);
    [row, col] = size(A.D);
    if k == row
        D = A.D; U = A.U; ind = zeros(1,0); cind = [1:col];
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
    %b = Q'*(LD\  b);
    % Retrieve part of the variables
    %x = b(1:end-k, :);
    %x = [x; zeros(col - size(x,1), size(x,2))];
    % Update transformations on the right
    F.U = [F.U, Q];
    F.L = [F.L, LD];
    F.ind = [F.ind, cur_ind];
    
    % Reduce the diagonal block
    D = eye(k);
    %b = b(end-k+1:end, :);
    
else    % NOT A LEAF
    [F, Dl, Ul, indl, cindl] = hss_chol_fact_rec(A.A11,  F, cur_ind(1:A.nl));                    % recursive call on left  child
    
    [F, Dr, Ur, indr, cindr] = hss_chol_fact_rec(A.A22,  F, cur_ind(A.nl+1:end));  % recursive call on right child
    
    %x = [xl; xr]; % Store the computed variables
    
    % Merge nodes
    D = [Dl, Ul * A.B12 * Ur'; Ur * A.B21 * Ul', Dr];
    %b = [bl; br];
    
    if A.topnode == 1 % If we are in the root we solve the remaining dense sytem and we apply right transformations
        U = [];
        F.L = [F.L, chol(D,'lower')];
        F.ind = [F.ind, [cur_ind(cindl), cur_ind(A.nl + cindr)]];
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
    %b = Q'*(LD \ b);
    
    % Retrieve part of the variables
    %z = b(1:end-k, :);
    % Update transformations on the right
    F.U = [F.U, Q];
    F.L = [F.L, LD];
    F.ind = [F.ind, [cur_ind(cindl), cur_ind(A.nl + cindr)]];
    
    nz = size(D,1) - k;
    n_cindl = length(cindl);
    
    if nz < n_cindl
        %x(cindl(1:nz), :) = z;
        %indl = union(indl, cindl(1:nz)); % updates indices of computed entries in the current level
        indl = [indl, cindl(1:nz)];
        cindl = cindl(nz+1:end);
    else
        %x(cindl,:) = z(1:n_cindl ,:);
        %x(A.nl + cindr(1: nz-n_cindl), :) = z(n_cindl + 1:end, :);
        indl = [1:A.ml];   % updates indices of computed entries in the current level
        indr = [indr, cindr(1: nz-n_cindl)];
        cindl = zeros(1,0);
        cindr = cindr(nz-n_cindl+1:end);
    end
    ind = [indl, A.nl + indr];  % updates indices of computed entries up to the current level
    cind = [cindl, A.nl + cindr]; % Complementary of ind;
    % Reduce the rows in the system
    D = eye(k);
    %b = b(end-k+1:end, :);
end
end
