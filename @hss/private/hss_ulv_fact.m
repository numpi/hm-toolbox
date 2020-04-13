function F = hss_ulv_fact(A)
% HSS_ULV_SOLV computes the  ULV factorization
%	       of A  taken from the paper:
%
%	       Chandrasekaran, Shiv, Ming Gu, and Timothy Pals.
%	       "A fast ULV decomposition solver for hierarchically semiseparable representations."
%	       SIAM Journal on Matrix Analysis and Applications 28.3 	(2006): 603-622.
%


F = struct();
F.A = A;
F.U = {};
F.V = {};
F.L = {};

if A.leafnode == 1
    return
end

F = hss_ulv_fact_rec(F, A, {}, {}, [1:size(A,1)]);
end



function [F, D, U, V, ind, cind,  ort, ort_ind] = hss_ulv_fact_rec(F, A, ort, ort_ind, cur_ind)
%------------------------------------------------------------------------------------------------------------------------------------------------
% Recursively compute the facators of ULV decomposition
%
% INPUT  : A coefficient matrix (HSS format), ort and ort_ind see description of outputs, cur_ind indices of the current submatrix
%
% OUTPUT : if the function is called on the root then it returns only the solution x
%	   otherwise it returns:
%
%          - F	       struct containing the sequences of trasformations
%	   - D         remaining compressed diagonal block
%	   - U         remaining rows of the column generator
%	   - V         updated column generator
%	   - ind       indices of the component of x, corresponding to computed entries up to now
%	   - cind      complementary of the set ind
%	   - ort       array of cells containing the orthogonal transformations on the right
%	   - ort_ind   array of cells containing the indices where the right orthogonal transformations have to be applied
%
%------------------------------------------------------------------------------------------------------------------------------------------------

if A.leafnode == 1  % LEAF
    k = size(A.U,2);
    [row, col] = size(A.D);
    if k == row
        D = A.D; U = A.U; V = A.V; ind = zeros(1,0); cind = [1:col];
        return
    end
    ind = [1:row-k];  % Indices of computed entries in the solution, at the end
    cind = [row - k + 1: row]; % Complementary of ind
    Q = zeros(row);
    % Compress off-diagonal block row
    [Q(:, end:-1:1), A.U(end:-1:1, :)] = qr(A.U);
    F.U = [F.U, Q];
    U = A.U(end-k+1:end,:);
    A.D = Q' * A.D;
    % Triangularize part of diagonal block
    [Q, L] = qr(A.D(1:end-k,:)');
    F.L = [F.L, {L(1:end-k,:)'} ];
    % Update transformations on the right
    V = Q' * A.V;
    ort = [ort, Q];
    ort_ind = [ort_ind, cur_ind];
    % Reduce the rows in the system
    D = A.D(end-k+1:end, :) * Q;
    F.L = [F.L, {D(:, 1:end-k)}];
    
    % Reduce the columns of the diagonal block
    D = D(:,end-k+1:end);
    
else    % NOT A LEAF
    [F,  Dl, Ul, Vl, indl, cindl, ort, ort_ind] = hss_ulv_fact_rec(F, A.A11, ort, ort_ind, cur_ind(1:A.nl));                    % recursive call on left  child
    
    [F,  Dr, Ur, Vr, indr, cindr, ort, ort_ind] = hss_ulv_fact_rec(F, A.A22, ort, ort_ind, cur_ind(A.nl+1:end));  % recursive call on right child
    
    % Reduce the columns in the system
    F.L = [F.L, {Ul} , {A.B12 * Vr(indr,:)'}];
    F.L = [F.L, {Ur} , {A.B21 * Vl(indl,:)'}];

    % Merge nodes
    D = [Dl, Ul * A.B12 * Vr(cindr, :)'; Ur * A.B21 * Vl(cindl, :)', Dr];
    
    if A.topnode == 1 % If we are in the root we solve the remaining dense sytem and we apply right transformations
        F.L = [F.L, {D}];
        U = []; V = [];
        % Apply all the transformations on the right
        for j = length(ort):-1:1
            F.V = [F.V, ort{j}, ort_ind{j}];
            % x(ort_ind{j}, :) = ort{j} * x(ort_ind{j}, :);
        end
        return
    end
    
    % Otherwise compress the off-diagonal block row in the current level
    V = [Vl * A.Wl; Vr * A.Wr];
    U = [Ul * A.Rl; Ur * A.Rr];
    k = size(U,2);
    if k == size(U,1) % full rank, no compression needed
        ind = [indl, A.nl + indr];
        cind = [cindl, A.nl + cindr];
        return
    end
    Q = zeros(size(U,1));
    [Q(:, end:-1:1), U(end:-1:1, :)] = qr(U);
    U = U(end-k+1:end,:);
    F.U = [F.U,Q];
    D = Q' * D;
    % Triangularize (1,1) submatrix of diagonal block
    [Q, L] = qr(D(1:end-k,:)');
    
    % Retrieve part of the variables
    F.L = [F.L, {L(1:end-k,:)'}];
    
    
    % Update transformations on the right
    V([cindl, A.nl + cindr],:) = Q' * V([cindl, A.nl + cindr],:);
    ort = [ort, Q];
    ort_ind = [ort_ind, [cur_ind(cindl), cur_ind(A.nl + cindr)]];
    
    nz = size(L,2);
    n_cindl = length(cindl);
    
    if nz < n_cindl
        % updates indices of computed entries in the current level
        indl = [indl, cindl(1:nz)];
        cindl = cindl(nz+1:end);
    else
        indl = [1:A.ml];   % updates indices of computed entries in the current level
        indr = [indr, cindr(1: nz-n_cindl)];
        cindl = zeros(1,0);
        cindr = cindr(nz-n_cindl+1:end);
    end
    ind = [indl, A.nl + indr];  % updates indices of computed entries up to the current level
    cind = [cindl, A.nl + cindr]; % Complementary of ind;
    % Reduce the rows in the system
    D = D(end-k+1:end, :) * Q;
    F.L = [F.L, {D(:, 1:end-k)}];
    
    % Reduce the columns of the diagonal block
    D = D(:,end-k+1:end);
end
end
