function x = hss_ulv_fact_solve(F,b)
% HSS_ULV_FACT_SOLV computes the solution of the linear system A X =B
%	       where F contains ULV decomposition of A
%
%


if F.A.leafnode == 1
    x = full(F.A) \ b;
    return
end

idx = [ 1 1 1 ];

x = hss_ulv_fact_solve_rec(F, F.A, [1:size(F.A,1)], b, idx);

end



function [x, F, ind, cind, b, idx] = hss_ulv_fact_solve_rec(F, A, cur_ind, b, idx)
%------------------------------------------------------------------------------------------------------------------------------------------------
% Recursively compute the solution of the linear system
%
% INPUT  : A coefficient matrix (HSS format), b right hand side, ort and ort_ind see description of outputs, cur_ind indices of the current submatrix
%
% OUTPUT : if the function is called on the root then it returns only the solution x
%	   otherwise it returns:
%
%	   - x         part of the computed solution
%          - F	       struct containing the sequences of trasformations
%	   - ind       indices of the component of x, corresponding to computed entries up to now
%	   - cind      complementary of the set ind
%	   - b         updated right hand side of the unsolved equations
%          - idx       array of indices for the stack of trasformations
%
%------------------------------------------------------------------------------------------------------------------------------------------------

if A.leafnode == 1  % LEAF
    k = size(A.U,2);
    [row, col] = size(A.D);
    if k == row
        x = zeros(col, size(b,2)); ind = zeros(1,0); cind = [1:col];
        return
    end
    ind = [1:row-k];  % Indices of computed entries in the solution, at the end
    cind = [row - k + 1: row]; % Complementary of ind
    % Compress off-diagonal block row
    b = F.U{idx(1)}' * b; idx(1) = idx(1) + 1;
    % Triangularize part of diagonal block
    x = F.L{idx(2)} \ b(1:end-k, :); idx(2) = idx(2) + 1;
    x = [x; zeros(col - size(x,1), size(x,2))];
    % Update transformations on the right
    b = b(end-k+1:end, :) - F.L{idx(2)} * x(1:end-k,:); idx(2) = idx(2) + 1;
    % Reduce the columns of the diagonal block
    
else    % NOT A LEAF
    [xl, F,   indl, cindl, bl, idx] = hss_ulv_fact_solve_rec(F, A.A11,  cur_ind(1:A.nl), b(1:A.ml, :), idx);                    % recursive call on left  child
    
    [xr, F,   indr, cindr, br, idx] = hss_ulv_fact_solve_rec(F, A.A22,  cur_ind(A.nl+1:end), b(A.ml + 1:A.ml + A.mr, :), idx);  % recursive call on right child
    
    x = [xl; xr]; % Store the computed variables
    
    % Reduce the columns in the system
    bl = bl - F.L{idx(2)} * (F.L{idx(2)+1} * xr(indr, :)); idx(2) = idx(2) + 2;
    br = br - F.L{idx(2)} * (F.L{idx(2)+1} * xl(indl, :)); idx(2) = idx(2) + 2;
    % Merge nodes
    b = [bl; br];
    if A.topnode == 1 % If we are in the root we solve the remaining dense sytem and we apply right transformations
        x([cindl, A.nl + cindr],:) = F.L{idx(2)} \ b; idx(2) = idx(2) + 1;
        % Apply all the transformations on the right
        for j = 1:2:length(F.V)
            x(F.V{idx(3)+1}, :) = F.V{idx(3)} * x(F.V{idx(3)+1}, :); idx(3) = idx(3) + 2;
        end
        return
    end
    
    % Otherwise compress the off-diagonal block row in the current level
    k = size(A.Rl,2);
    if k == length(cindl)+length(cindr) % full rank, no compression needed
        ind = [indl, A.nl + indr];
        cind = [cindl, A.nl + cindr];
        return
    end
    b = F.U{idx(1)}' * b; idx(1) = idx(1) + 1;
    % Retrieve part of the variables
    z = F.L{idx(2)} \ b(1:end-k, :); idx(2) = idx(2) + 1;
    
    % Update transformations on the right
    
    nz = size(z,1);
    n_cindl = length(cindl);
    
    if nz < n_cindl
        x(cindl(1:nz), :) = z;
        indl = [indl, cindl(1:nz)];
        cindl = cindl(nz+1:end);
    else
        x(cindl,:) = z(1:n_cindl ,:);
        x(A.nl + cindr(1: nz-n_cindl), :) = z(n_cindl + 1:end, :);
        indl = [1:A.ml];   % updates indices of computed entries in the current level
        indr = [indr, cindr(1: nz-n_cindl)];
        cindl = zeros(1,0);
        cindr = cindr(nz-n_cindl+1:end);
    end
    ind = [indl, A.nl + indr];  % updates indices of computed entries up to the current level
    cind = [cindl, A.nl + cindr]; % Complementary of ind;
    % Reduce the rows in the system
    b = b(end-k+1:end, :) - F.L{idx(2)} * z; idx(2) = idx(2) + 1;
end
end
