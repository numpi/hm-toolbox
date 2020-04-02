function [X,Q,Z, Y0] = hss_mldivide(A, B)
if A.leafnode == 1 % Base case
    X = B;
    X.D = full(A) \ full(B);
    return
end

if size(A, 1) == 0
    % This is a special case -- we use this just to copy the HSS tree
    % in place, since there is no computation to do.
    X = A' * B;
    return;
end

[L, Q, Z, ind, cind] = hss_mldivide_rec(A, {}, {}, {}, {});

QB = applyQ(B, Q, true);
Y0 = applyLinv(QB, L, ind, cind);
L = hss_project(L, cind, 'row');

% This performs QB - L*Y0
LY = L * Y0;
QB = hss_project(QB, cind, 'row');
QB = hss_sum(QB, -LY, true);
QB = hss_remove_leaf_level(QB);

% Project L on the columns as well
L = hss_project(L, cind, 'col');
L = hss_remove_leaf_level(L);
Y1 = hss_mldivide(L, QB);

Y1 = hss_unpack(Y0, Y1, ind, cind);

Y = hss_sum(Y0, Y1, false);

X = hss_compress(applyQ(Y, Z, false), hssoption('threshold'));

end

function [A, Q, Z, ind, cind] = hss_mldivide_rec(A, Q, Z, ind, cind)
if A.leafnode == 1  % LEAF
    k = size(A.U,2);
    [row, col] = size(A.D);
    
    if k >= row
        Q{length(Q) + 1} = eye(col);
        Z{length(Z) + 1} = eye(row);
        % Q = [Q, eye(col)]; Q = [Z, eye(row)];
        ind{length(ind)+1} = zeros(1,0);
        cind{length(cind)+1} = (1:col);
        return
    end
    
    % ind = [ind, [1:row-k]];  % Indices of computed entries in the solution, at the end
    ind{length(ind)+1} = [1:row-k];
    
    % cind = [cind, [row - k + 1: row]]; % Complementary of ind
    cind{length(cind)+1} = [row - k + 1 : row];
    
    QQ = zeros(row);
    % Compress off-diagonal block row
    [QQ(:, end:-1:1), A.U(end:-1:1, :)] = qr(A.U);
    Q = [Q, QQ];
    A.D = QQ' * A.D;
    
    % Triangularize part of diagonal block
    [ZZ, L] = qr(A.D(1:end-k,:)');
    A.D(end-k+1:end,:) = A.D(end-k+1:end,:) * ZZ;
    A.D(1:end-k,:) = L';
    
    % Update transformations on the right
    A.V = ZZ' * A.V;
    Z = [Z, ZZ];
else    % NOT A LEAF
    [A.A11, Q, Z, ind, cind] = hss_mldivide_rec(A.A11, Q, Z, ind, cind);  % recursive call on left  child
    [A.A22, Q, Z, ind, cind] = hss_mldivide_rec(A.A22, Q, Z, ind, cind);  % recursive call on right child
end
end

function  [X, Q, ind] = applyQ(X, Q, transp, ind)
if ~exist('ind', 'var')
    ind = {};
end
if X.leafnode == 1  % LEAF
    QQ = Q{1}; Q = { Q{2:end} };
    
    if ~isempty(ind)
        localind = ind{1};
        ind = { ind{2:end} };
        
        if transp
            X.D(localind,:) = QQ(localind,:)' * X.D(localind,:);
            X.U(localind,:) = QQ(localind,:)' * X.U(localind,:);
        else
            X.D(localind,:) = QQ(:,localind) * X.D(localind,:);
            X.U(localind,:) = QQ(:,localind) * X.U(localind,:);
        end
    else
        if transp
            X.D = QQ' * X.D;
            X.U = QQ' * X.U;
        else
            X.D = QQ * X.D;
            X.U = QQ * X.U;
        end
    end
else
    [X.A11, Q, ind] = applyQ(X.A11, Q, transp, ind);
    [X.A22, Q, ind] = applyQ(X.A22, Q, transp, ind);
end
end

function  [X, ind, cind] = applyLinv(X, L, ind, cind)
if ~exist('ind', 'var')
    ind = {};
end
if X.leafnode == 1  % LEAF
    LL = L.D;
    
    if ~isempty(ind)
        localind = ind{1};
        ind = { ind{2:end} };
        localcind = cind{1};
        cind = { cind{2:end} };
    else
        localind = 1 : size(X.D, 1);
        localcind = [];
    end
    
    X.D(localind,:) = LL(localind,localind)\X.D(localind,:);
    X.U(localind,:) = LL(localind,localind)\X.U(localind,:);
    X.D(localcind,:) = 0;
    X.U(localcind,:) = 0;
    
else
    [X.A11, ind, cind] = applyLinv(X.A11, L.A11, ind, cind);
    [X.A22, ind, cind] = applyLinv(X.A22, L.A22, ind, cind);
end
end

function [A, ind] = hss_project(A, ind, rc)
% Compute principal projection of an HSS matrix A(ind, :)
if A.leafnode == 1  % LEAF
    if ~isempty(ind)
        localind = ind{1};
        ind = { ind{2:end} };
    else
        localind = [];
    end
    
    if strcmp(rc, 'row')
        A.D = A.D(localind, :);
        A.U = A.U(localind, :);
    elseif strcmp(rc, 'col')
        A.D = A.D(:, localind);
        A.V = A.V(localind, :);
    else
        A.D = A.D(localind, localind);
        A.U = A.U(localind, :);
        A.V = A.V(localind, :);
    end
else
    [A.A11, ind] = hss_project(A.A11, ind, rc);
    [A.A22, ind] = hss_project(A.A22, ind, rc);
    [A.ml, A.nl] = size(A.A11);
    [A.mr, A.nr] = size(A.A22);
end
end

function [Y1, ind, cind] = hss_unpack(Y0, Y1, ind, cind)
if Y1.leafnode == 1 && Y0.leafnode == 1 
	ind1 = ind{1};
        lcind1 = cind{1};
        cind = { cind{2:end} };
        ind = { ind{2:end} };
        temp = zeros(size(Y0));
        temp(lcind1, :) = Y1.D; Y1.D = temp;
	if size(Y0, 1) > size(Y1.U, 1)
		temp = zeros(size(Y0, 1), size(Y1.U, 2));
		temp(lcind1, :) = Y1.U; 
		Y1.U = temp;
	end	
	%if size(Y0, 2) > size(Y1.V, 1)
	%	temp = zeros(size(Y0, 2), size(Y1.V, 2));
	%	temp(lcind1, :) = Y1.V; 
	%	Y1.V = temp;
	%end
	return
end

if Y1.leafnode == 1 
    Y1.leafnode = 0;
    Y1.A11 = hss();
    Y1.A22 = hss();
    Y1.A11.leafnode = 1;
    Y1.A22.leafnode = 1;
    Y1.A11.topnode = 0;
    Y1.A22.topnode = 0;
    
    lcind1 = cind{1};
    lcind2 = cind{2};
    ind1 = ind{1};
    ind2 = ind{2};
    
    cind = { cind{3:end} };
    ind = { ind{3:end} };
    Y1.A11.D = zeros(size(Y0.A11));
    Y1.A11.D(lcind1, :) = Y1.D(1:length(lcind1), 1:size(Y0.A11,2));
    Y1.A22.D = zeros(size(Y0.A22));
    Y1.A22.D(lcind2, :) = Y1.D(length(lcind1)+1:end, size(Y0.A11,2)+1:end);
    
    if Y1.topnode == 1
        Y1.U = zeros(length(lcind1) + length(lcind2), 0);
        Y1.V = zeros(size(Y0, 2), 0);
    end
    
    % A12 block
    Y1.A11.U = zeros(size(Y0.A11, 1), size(Y1.U,2) + length(lcind1));
    Y1.A11.U(end-length(lcind1)+1:end, :) = [ Y1.U(1:length(lcind1), :), eye(length(lcind1)) ];
    Y1.A22.V = [ Y1.V(size(Y0.A11.V,1)+1:end,:), Y1.D(1:length(lcind1), size(Y0.A11.V,1)+1:end)' ];
    Y1.B12 = zeros(size(Y1.A11.U, 2), size(Y1.A22.V, 2)); Y1.B12(end-length(lcind1)+1:end,end-length(lcind1)+1:end) = eye(length(lcind1));
    
    % A21 block
    Y1.A22.U = zeros(size(Y0.A22, 1), size(Y1.U, 2) + length(lcind2));
    Y1.A22.U(end-length(lcind2)+1:end, :) = [ Y1.U(length(lcind1)+1:end, :), eye(length(lcind2)) ];
    Y1.A11.V = [ Y1.V(1:size(Y0.A11.V,1),:), Y1.D(length(lcind1)+1:end, 1:size(Y0.A11.V,1))' ];
    Y1.B21 = zeros(size(Y1.A22.U, 2), size(Y1.A11.V, 2)); Y1.B21(end-length(lcind2)+1:end,end-length(lcind2)+1:end) = eye(length(lcind2));
    
    Y1.Rl = eye(size(Y1.A11.U, 2), size(Y1.U,2));
    Y1.Rr = eye(size(Y1.A22.U, 2), size(Y1.U,2));
    Y1.Wl = eye(size(Y1.A11.V, 2), size(Y1.V,2));
    Y1.Wr = eye(size(Y1.A22.V, 2), size(Y1.V,2));
    
    Y1.U = []; Y1.V = []; Y1.D = [];
    
    [Y1.ml, Y1.nl] = size(Y1.A11);
    [Y1.mr, Y1.nr] = size(Y1.A22);
else
    [Y1.A11, ind, cind] = hss_unpack(Y0.A11, Y1.A11, ind, cind);
    [Y1.A22, ind, cind] = hss_unpack(Y0.A22, Y1.A22, ind, cind);
    [Y1.ml, Y1.nl] = size(Y1.A11);
    [Y1.mr, Y1.nr] = size(Y1.A22);
end
end
