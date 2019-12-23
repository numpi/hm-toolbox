function H = hss_from_full(obj, A)
%HSS_FROM_FULL Build an HSS representation of a dense matrix.
%
% H = HSS_FROM_FULL(A) returns the HSS representation of A, where
%     truncation is done at the relative tolerance indicated by
%     hssoption('threshold'), using the compression method indicated by
%     hssoption('svd').
%
% This functions is equivalent to calling H = HSS(A) with no other options,
% which is the preferred syntax.
%
% The implemenetation and the algorithodlr is based on the one presented in
% the paper
%
%     Xia, Jianlin, et al. "Fast algorithodlrs for hierarchically
%     semiseparable matrices." Numerical Linear Algebra with Applications
%     17.6 (2010): 953-976.
%
% The complexity of the algorithodlr, assuming a low HSS rank, is quadratic in
% the dimension.

% When called with one variable only we are not invoked as a method, therefore,
% we shall set A equal to the first parameter.
if ~exist('A', 'var')
    A = obj;
end

m = size(A, 1);
n = size(A, 2);

if m ~= n
    % error('Rectangular HSS matrices are not (yet) supported');
end

tol = hssoption('threshold');

% Prepare the tree for the HSS structure -- leaving all the blocks empty
% H = hss_build_hss_tree(m, n, block_size);
H = obj;

% Create the stack
rs = createStack();
cs = createStack();

row_offsets = [];
col_offsets = [];

H = BuildHSS_iter(H, A, tol, rs, cs, ...
    row_offsets, col_offsets, 0, 0);

end

function [H, rs, cs, row_offsets, col_offsets, rs1, cs1] = BuildHSS_iter(...
    H, A, tol, rs, cs, row_offsets, col_offsets, mh, nh)

m = size(H, 1);
n = size(H, 2);

if H.leafnode
    % You might want to uncomment this for debugging purposes
    % fprintf('> Leafnode, mh = %d, nh = %d, m = %d, n = %d\n', mh, nh, m, n);
    
    % We need to extract all the blocks in the stack, and put them
    % together to assemble the matrix.
    
    % Let's do it for the HSS block row first.
    B = A(mh+1:mh+m, nh+n+1:end);
    
    for j = length(col_offsets) - 1 : -1 : 0
        Z = hss_elementAt(cs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+m,:) , B ];
    end
    
    [H.U, Z] = compress_hss_block(B, tol);
    
    % Push Z in the stack, and updates the elements in there as well
    counter = 0;
    for j = 0 : length(col_offsets) - 1
        W = hss_elementAt(cs, j);
        
        cs = hss_setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(col_offsets(j+1)+m+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
        % col_offsets(j+1) = size(Z, 2);
    end
    
    % Last but not least, push the new compressed stuff into the stack,
    % and add an offset entry for it
    rs = hss_push(rs, Z(counter+1:end, :));
    rs1 = size(Z, 2);
    
    % And now do the columns
    B = A(mh+m+1:end, nh+1:nh+n)';
    
    for j = length(row_offsets) - 1 : -1 : 0
        Z = hss_elementAt(rs, j);
        
        % FIXME: We should preallocate B
        B2 = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+n,:) , B ];
        B = B2;
    end
    
    [H.V, Z] = compress_hss_block(B, tol);
    
    % Push Z in the stack, and updates the elements in the as well
    counter = 0;
    for j = 0 : length(row_offsets) - 1
        W = hss_elementAt(rs, j);
        
        rs = hss_setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(row_offsets(j+1)+n+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
    end
    
    % Last but not least, push the new compressed stuff into the stack,
    % and add an offset entry for it
    cs = hss_push(cs, Z(counter+1:end, :));
    cs1 = size(Z, 2);
    
    col_offsets = [ col_offsets , 0 ];
    row_offsets = [ row_offsets , 0 ];
    
    H.D = A(mh+1:mh+m, nh+1:nh+n);
else
    % fprintf('> Non-leafnode, mh = %d, nh = %d, m = %d, n = %d\n', mh, nh, m, n);
    % Call the constructor recursively on the left and right childs
    [H.A11, rs, cs, row_offsets, col_offsets] = BuildHSS_iter(...
        H.A11, A, tol, rs, cs, ...
        row_offsets, col_offsets, mh, nh);
    
    [H.A22, rs, cs, row_offsets, col_offsets] = BuildHSS_iter(...
        H.A22, A, tol, rs, cs, row_offsets, col_offsets, ...
        mh + size(H.A11, 1), nh + size(H.A11, 2));
    
    % Extract Bl and Bu from the stacks, and merge the children
    [row1, rs] = hss_pop(rs); [row2, rs] = hss_pop(rs);
    [col1, cs] = hss_pop(cs); [col2, cs] = hss_pop(cs);
    
    rs1 = size(row2, 1) - size(row1, 1);
    cs1 = size(col2, 1) - size(col1, 1);
    
    [rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H);
    
    assert(rs1 == rrV)
    assert(cs1 == rrU)
    
    H.B12 = row2(1:rrV, 1:rlU)';
    H.B21 = col2(1:rrU, 1:rlV);
    
    if H.topnode
        return;
    end
    
    if length(col_offsets) > 2
        col_offsets = col_offsets(1:end-2) - col_offsets(end-2);
    else
        col_offsets = [];
    end
    
    if length(row_offsets) > 2
        row_offsets = row_offsets(1:end-2) - row_offsets(end-2);
    else
        row_offsets = [];
    end
    
    % Merge the rows and cols
    B = [ row2(rrV+1:end, :), row1 ]';
    
    for j = length(col_offsets) - 1 : -1 : 0
        Z = hss_elementAt(cs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+rU,:) , B ];
    end
    
    [U, Z] = compress_hss_block(B, tol);
    
    H.Rl = U(1:rlU, :);
    H.Rr = U(rlU+1:end, :);
    
    % Push Z in the stack, and updates the elements in there as well
    counter = 0;
    for j = 0 : length(col_offsets) - 1
        W = hss_elementAt(cs, j);
        
        cs = hss_setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(col_offsets(j+1)+rU+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
    end
    
    rs = hss_push(rs, Z(counter+1:end, :));
    
    B = [ col2(rrU+1:end, :), col1 ]';
    
    for j = length(row_offsets) - 1 : -1 : 0
        Z = hss_elementAt(rs, j);
        
        % FIXME: We should preallocate B
        B = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+rV,:) , B ];
    end
    
    [U, Z] = compress_hss_block(B, tol);
    
    H.Wl = U(1:rlV, :);
    H.Wr = U(rlV+1:end, :);
    
    % Push Z in the stack, and updates the elements in the as well
    counter = 0;
    for j = 0 : length(row_offsets) - 1
        W = hss_elementAt(rs, j);
        
        rs = hss_setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
            Z(counter+1:counter+size(W, 2), :)' ; ...
            W(row_offsets(j+1)+rV+1:end,:) ], j);
        
        counter = counter + size(W, 2);
        row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
    end
    
    cs = hss_push(cs, Z(counter+1:end, :));
    
    row_offsets = [ row_offsets, 0 ];
    col_offsets = [ col_offsets, 0 ];
end
end

function [rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H)
if H.A11.leafnode
    rlU = size(H.A11.U, 2);
    rlV = size(H.A11.V, 2);
else
    rlU = size(H.A11.Rl, 2);
    rlV = size(H.A11.Wl, 2);
end

if H.A22.leafnode
    rrU = size(H.A22.U, 2);
    rrV = size(H.A22.V, 2);
else
    rrU = size(H.A22.Rl, 2);
    rrV = size(H.A22.Wl, 2);
end

rU = rlU + rrU;
rV = rlV + rrV;
end



function [U, Z] = compress_hss_block(B, tol)

if isempty(B)
    U = zeros(size(B, 1), 0);
    Z = zeros(size(B, 2), 0);
    return
end

switch hssoption('compression')
    case 'qr'
        [Q, R, P] = prrqr(B, tol);        
        IP = zeros(1, length(P)); IP(P) = 1 : length(P);        
        U = Q;
        Z = R(:, IP)';
        
    case 'svd'
        [U, B, Z] = svd(B, 'econ');
        
        if tol == 0
            rk = size(B, 1);
        else
            switch hssoption('norm')
                case 2
                    rk = sum(diag(B) > tol * B(1,1));
                case 'fro'
                    t = diag(B);
                    tt = sqrt(cumsum(t(end:-1:1).^2));
                    rk = sum(tt > tol * tt(end));
            end
        end
        
        U = U(:, 1:rk);
        B = B(1:rk, 1:rk);
        
        Z = Z(:, 1:rk) * B;
end
end

% Simple home-made implementation of a stack -- if we really want to use
% this is still open for discussion. It does not matter much
% performance-wise in this context anyway.

function s = createStack()
s = {};
end

function s = hss_push(s, el)
s{end+1} = el;
end

function [el, s] = hss_pop(s)
el = s{end};
s = { s{1:end-1} };
end

function el = hss_elementAt(s, j)
el = s{j+1};
end

function s = hss_setElementAt(s, el, j)
s{j+1} = el;
end


