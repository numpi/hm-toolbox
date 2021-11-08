function H = halr_from_matvec_adaptive(Afun, ATfun, maxrank, m, n)
%---------------------------------------------------------------------
%
% Computes the halr representation of the matrix with assigned cluster 
% by means of matvecs
%
%	Afun, ATfun		handle functions that evaluate A*x and A'*x
%	maxrank			maximum rank allowed for the the low-rank blocks
%	cl				cluster tree
%
%---------------------------------------------------------------------
	tol = halroption('threshold');
	% H = halr('low-rank', zeros(size(cl, 1), 0), zeros(size(cl, 2), 0), 'cluster', cl);
	% matr = find_cluster_matrices(cl, {});
	% l = length(matr);
    
    
    cl = halr; cl.sz = [m n]; cl.admissible = 1;
    H = halr('low-rank', zeros(m,0), zeros(n,0), 'cluster', cl);
    l = 1;
   
    matr = cell(1, l);
    matr{1} = 0;
    
    % Generate the random vector for sampling the first level
    vectors = { randn(n, maxrank) };
    pos = { 1 };
    
    A0 = { Afun(vectors{1}) };
    
    [A0{1}, R0] = qr(A0{1}, 0);
    
    % FIXME: Make a sensible check
    if min(svd(R0)) <= tol
        % Compute the low-rank representation
        temp = ATfun(A0{1});
        temp = temp';
        
        [U, V] = compress_matrix(temp);
        U = A0{1} * U; 
        H.U = U; H.V = V;
    else
        matr{1} = 1; H.admissible = 0;
    end
    
	while ~can_stop(matr, l, H)
        H = add_one_level(H, l);
        l = l + 1;
        matr{l} = kron(matr{l-1}, zeros(2));
        
		[vectors, pos] = gen_random_sampling(matr, l, H, maxrank);
        
   		partialT = @(x) H'*x;
        partial  = @(x) H*x;
        
		AO = cell(1, length(vectors));
		for h = 1:length(vectors)
			AO{h} = Afun(vectors{h}) - partial(vectors{h});
			[I, J] = find(pos == h);
			for t = 1:length(I) % Orthogonalize the right products
				[rind, cind] = indices_from_path(dec2bin(I(t) - 1, l - 1), dec2bin(J(t) - 1, l - 1), H);
				[AO{h}(rind, :), R0] = qr(AO{h}(rind, :), 0);
                
                if min(svd(R0)) > tol
                    matr{l}(I(t), J(t)) = 1;
                    % Set the non-admissible flag
                    H = save_to_path(dec2bin(I(t) - 1, l - 1), dec2bin(J(t) - 1, l - 1), H);
                end
			end
        end
        
		[lvectors, lpos] = gen_basis_vectors(matr, l, H, pos, AO, maxrank);
		for h = 1:length(lvectors)
			temp = ATfun(lvectors{h}) - partialT(lvectors{h});	
			temp = temp';
			[J, I] = find(lpos == h);
			for t = 1:length(I)
				[rind, cind] = indices_from_path(dec2bin(I(t) - 1, l - 1), dec2bin(J(t) - 1, l - 1), H);
				[U, V] = compress_matrix(temp(:, cind));
				ht = pos(I(t), J(t));
				U = AO{ht}(rind, :) * U;
				H = save_to_path(dec2bin(I(t) - 1, l - 1), dec2bin(J(t) - 1, l - 1), H, U, V);
			end
        end        
    end
    
    H = add_dense_admissible_blocks(H);
    
	% Retrieve dense blocks at the end
	[vectors, pos] = gen_dense_sampling(matr, H);
	partial = @(x) H*x;	
	for h = 1:length(vectors)
		temp = Afun(vectors{h}) - partial(vectors{h});	
		[I, J] = find(pos == h);
		for t = 1:length(I)
				[rind, cind] = indices_from_path(dec2bin(I(t) - 1, l - 1), dec2bin(J(t) - 1, l - 1), H);
				F = temp(rind, 1:length(rind));
				H = save_to_path(dec2bin(I(t) - 1, l-1), dec2bin(J(t) - 1, l-1), H, F);
		end
	end
	
	
end

%------------------Auxiliary functions-----------------------------------
%----------------------------------------------------------------------------------------
function matr = find_cluster_matrices(cl, matr)
% Determine the auxiliary matrices that describe the partitioning at each level
% Ex: for a HODLR with 3 levels we return the cell array {1, eye(2), eye(4), eye(8)}
% because zeros correspond to low-rank blocks and ones to dense or recursively partioned ones
	if is_leafnode(cl)
		if cl.admissible
			matr = [matr, 0];
		else
			matr = [matr, 1];
		end
	else
		matr11 = find_cluster_matrices(cl.A11, {});
		matr12 = find_cluster_matrices(cl.A12, {});
		matr21 = find_cluster_matrices(cl.A21, {});
		matr22 = find_cluster_matrices(cl.A22, {});
		mx = max([length(matr11), length(matr12),length(matr21),length(matr22)]);
		matr11 = extend_matr(matr11, mx);
		matr12 = extend_matr(matr12, mx);
		matr21 = extend_matr(matr21, mx);
		matr22 = extend_matr(matr22, mx);
		for j = 1:mx
			matr{j} = [matr11{j}, matr12{j}; matr21{j}, matr22{j}];
		end
		matr = [1, matr];
	end
end
%------------------------------------------------------------------------
function matr = extend_matr(matr, mx)
% Generate the auxiliary matrices obtained by padding to mx levels 
	for j  = length(matr) + 1:mx
		matr{j} = kron(matr{j-1}, ones(2));
	end
end
%--------------------------------------------------------------------------
function [vectors, pos] = gen_random_sampling(matr, j, cl, maxrank)
% Generate the vectors for Martinsson's algorithm, for level j, i.e. vectors whose sparsity patterns is
% described with the auxiliary matrices 
	[vectors, pos] = compute_symbol_vectors(matr, j);
	% We reduce the number of vectors required by leveraging the free entries marked with 2
	% the correspondence in the matrix pos is updated accordingly
	[vectors, pos] = merge_vec(vectors, pos, 'r'); 
	for h = 1:length(vectors)
		w = zeros(size(cl, 1), maxrank);
		[I, J] = find(pos == h);
		for t = 1:length(I)
		% Applying the function dec2bin to the entries of I and J returns binary strings that 
		% encode the path to get the corresponding block starting from the root, 
			[rind, cind] = indices_from_path(dec2bin(I(t) - 1, j - 1), dec2bin(J(t) - 1, j - 1), cl);
			w(cind, :) = randn(length(cind), maxrank);
		end
		vectors{h} = w;
	end
	
end
%---------------------------------------------------------------------------
function [vectors, pos] = compute_symbol_vectors(matr, j)
	if j <= length(matr)
		P1 = kron(matr{j-1}, ones(2)); % P1 has ones for low-rank blocks of level j and for recursively partioned ones
		P2 = P1 - matr{j}; % P2 has ones only for the low-rank blocks of level j
	else
		P2 = matr{end}; % special case for retrieving the dense blocks
		P1 = P2;
	end
	[I, J] = find(P2);
	vectors = {};
	pos = P2;
	for h = 1:length(I) 
	% we generate the structure of the vector required for block I(h),J(h)
	% the symbolic vectors contains 0 where zeros are imposed, 1 where random Gaussian
	% entries are required and 2 where there are no constraints
		v = 2 * ones(size(P2,2), 1);
		v(find(P1(I(h), :))) = 0;
		v(J(h)) = 1; 		
		vectors= [vectors, v];
		pos(I(h), J(h)) = h;
	end
end
%------------------------------------------------------------------------------
function [vectors, pos] = merge_vec(vectors, pos, lr)
% Determine the 'minimal' number of vector structures we required for random sampling 
	j = 1;
	while j < length(vectors)
		h = j + 1;
		while h <= length(vectors)
			if (~(any(vectors{j} == 0 & vectors{h} == 1) || any(vectors{j} == 1 & vectors{h} == 0)) && lr == 'r' ) || (~(any(vectors{j} ~= 2 & vectors{h} == 1) || any(vectors{j} == 1 & vectors{h} ~= 2)) && lr == 'l') 
				w = vectors{j};
				w(vectors{h} == 0) = 0;
				w(vectors{h} == 1) = 1;
				vectors{j} = w;
				vectors(h) = [];
				pos(pos == h) = j;
				pos(pos > h) = pos(pos > h) - 1;
			else
				h = h + 1;
			end
		end
		j = j + 1;
	end
end
%-------------------------------------------------------------------------------
function [rind, cind] = indices_from_path(I, J, cl)
% Obtain the row and column indices of a block in cluster described with the path I,J
% I and J are binary strings describing the choices of the children at the various levels starting
% from the root
	if isempty(I)
		rind = 1:size(cl, 1);
		cind = 1:size(cl, 2);
	else
		roffset = 0;
		coffset = 0;
		if I(1) == '1'
			roffset = size(cl.A11, 1); 
		end
		if J(1) == '1'
			coffset = size(cl.A11, 2);
		end
		[rind, cind] = indices_from_path(I(2:end), J(2:end), ...
		subsref(cl, struct('type', '.', 'subs', ['A', char(I(1) + 1), char(J(1) +1)])));
		rind = rind + roffset;
		cind = cind + coffset;
	end

end
%---------------------------------------------
function U = colspan(A, tol)
	[U, S, ~] = svd(full(A), 'econ'); 
        k = sum(abs(diag(S)) > S(1,1) * tol);	
        U = U(:,1:k);
end
%---------------------------------------------------------------
function [vectors, lpos] = gen_basis_vectors(matr, j, cl, pos, AO, maxrank)
	% Generate the row vectors for Martinsson's algorithm, for level j, i.e. vectors whose sparsity patterns is
% described with the auxiliary matrices 
	for h = 1:length(matr)
		matr{h} = matr{h}';
	end
	[vectors, lpos] = compute_symbol_vectors(matr, j);
	% We reduce the number of vectors required by leveraging the free entries marked with 2
	% the correspondence in the matrix pos is updated accordingly
	[vectors, lpos] = merge_vec(vectors, lpos, 'l'); 
	for h = 1:length(vectors)
		w = zeros(size(cl, 2), maxrank);
		[J, I] = find(lpos == h);
		for t = 1:length(J)
		% Applying the function dec2bin to the entries of I and J returns binary strings that 
		% encode the path to get the corresponding block starting from the root, 
			[rind, cind] = indices_from_path(dec2bin(I(t) - 1, j - 1), dec2bin(J(t) - 1, j - 1), cl);
			ht = pos(I(t), J(t));
			w(rind, :) = AO{ht}(rind, :);
		end
		vectors{h} = w;
	end	
end
%----------------------------------------------
function H = save_to_path(I, J, H, varargin)
% Set the dense or the low-rank factorization U*V' into the block corresponding to the path I,J
	if isempty(I)
        switch length(varargin)
            case 0
                H.admissible = false;
            case 1
                H.F = varargin{1};
            case 2
                H.U = varargin{1};
                H.V = varargin{2};
                H.admissible = true;
            otherwise
                error('Unsupported call to save_to_path');
		end
	else
		res = subsref(H, struct('type', '.', 'subs', ['A', char(I(1) + 1), char(J(1) +1)]));
		res = save_to_path(I(2:end), J(2:end), res, varargin{:});	
		H = subsasgn(H, struct('type', '.', 'subs', ['A', char(I(1) + 1), char(J(1) +1)]), res);	
	end
end
%----------------------------------------------
function [vectors, pos] = gen_dense_sampling(matr, cl)
% Generate the vectors for extracting the dense blocks at the final level
% described with the auxiliary matrices 
	l = length(matr);
	[vectors, pos] = compute_symbol_vectors(matr, length(matr) + 1);
	% We reduce the number of vectors required by leveraging the free entries marked with 2
	% the correspondence in the matrix pos is updated accordingly
	[vectors, pos] = merge_vec(vectors, pos, 'r'); 
	for h = 1:length(vectors)
		w = zeros(size(cl, 1), 0);
		[I, J] = find(pos == h);
		for t = 1:length(I)
		% Applying the function dec2bin to the entries of I and J returns binary strings that 
		% encode the path to get the corresponding block starting from the root, 
			[rind, cind] = indices_from_path(dec2bin(I(t) - 1, l-1), dec2bin(J(t) - 1, l-1), cl);
			w(cind, 1:length(cind)) = eye(length(cind));
		end
		vectors{h} = w;
	end
	
end

function r = can_stop(matr, l, H)
    % FIXME: We should check if the blocks are smaller than the minimal
    % block size using H, instead of just estimating by recursive splitting
    if l > log2(min(size(H)) / halroption('block-size'))
        r = true;
    else
        r = all(matr{l} == 0);
    end
end

function H = add_one_level(H, l)
    if l == 1
        if ~H.admissible
            m = size(H, 1); 
            n = size(H, 2);
            m1 = ceil(m / 2);
            n1 = ceil(n / 2);
            
            H.A11 = halr; H.A11.sz = [m1, n1];     H.A11.admissible = true;
            H.A11.U = zeros(m1, 0); H.A11.V = zeros(n1, 0);
            H.A12 = halr; H.A12.sz = [m1, n-n1];   H.A12.admissible = true;
            H.A12.U = zeros(m1, 0); H.A12.V = zeros(n-n1, 0);
            H.A21 = halr; H.A21.sz = [m-m1, n1];   H.A21.admissible = true;
            H.A21.U = zeros(m-m1, 0); H.A21.V = zeros(n1, 0);
            H.A22 = halr; H.A22.sz = [m-m1, n-n1]; H.A22.admissible = true;
            H.A22.U = zeros(m-m1, 0); H.A22.V = zeros(n-n1, 0);
        end
    else
        H.U = []; H.V = [];
        
        if ~H.A11.admissible
            H.A11 = add_one_level(H.A11, l-1);
        end
        
        if ~H.A12.admissible
            H.A12 = add_one_level(H.A12, l-1);
        end
        
        if ~H.A21.admissible
            H.A21 = add_one_level(H.A21, l-1);
        end
        
        if ~H.A22.admissible
            H.A22 = add_one_level(H.A22, l-1);
        end
    end
end

function H = add_dense_admissible_blocks(H)
    if is_leafnode(H)
        if ~H.admissible && isempty(H.F)
            H.F = zeros(size(H));
        end
    else
        H.A11 = add_dense_admissible_blocks(H.A11);
        H.A12 = add_dense_admissible_blocks(H.A12);
        H.A21 = add_dense_admissible_blocks(H.A21);
        H.A22 = add_dense_admissible_blocks(H.A22);
    end
end

