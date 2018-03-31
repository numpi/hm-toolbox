function H = hss_from_full(A)
%HSS_FROM_FULL Builg an HSS representation of a dense matrix. 

m = size(A, 1);
n = size(A, 2);

if m ~= n
    error('Rectangular HSS matrices are not (yet) supported');
end

block_size = hssoption('block-size');
tol = hssoption('threshold');

% Prepare the tree for the HSS structure -- leaving all the blocks empty
H = build_hss_tree(m, n, block_size);
H.topnode = 1;

% Create the stack
if exist('java', 'file')
    % rs = java.util.Stack();
    % cs = java.util.Stack();
	rs = createStack();
	cs = createStack();
else
	rs = createStack();
	cs = createStack();
    % rs = javaObject('java.util.Stack');
    % cs = javaObject('java.util.Stack');
end

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
            Z = elementAt(cs, j);
            
            % FIXME: We should preallocate B
            B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+m,:) , B ];
        end
        
        [H.U, Z] = compress_hss_block(B, tol);
        
        % Push Z in the stack, and updates the elements in the as well
        counter = 0;
        for j = 0 : length(col_offsets) - 1
            W = elementAt(cs, j);
            
            cs = setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
                Z(counter+1:counter+size(W, 2), :)' ; ...
                W(col_offsets(j+1)+m+1:end,:) ], j);
            
            counter = counter + size(W, 2);
            col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
			% col_offsets(j+1) = size(Z, 2);
        end
        
        % Last but not least, push the new compressed stuff into the stack,
        % and add an offset entry for it
        rs = push(rs, Z(counter+1:end, :));
		rs1 = size(Z, 2);
        
        % And now do the columns 
        B = A(mh+m+1:end, nh+1:nh+n)';
        
        for j = length(row_offsets) - 1 : -1 : 0
            Z = elementAt(rs, j);
            
            % FIXME: We should preallocate B
            B = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+m,:) , B ];
        end
        
        [H.V, Z] = compress_hss_block(B, tol);
        
        % Push Z in the stack, and updates the elements in the as well
        counter = 0;
		for j = 0 : length(row_offsets) - 1
			W = elementAt(rs, j);

			rs = setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
				Z(counter+1:counter+size(W, 2), :)' ; ...
				W(row_offsets(j+1)+n+1:end,:) ], j);

			counter = counter + size(W, 2);
			row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
		end
        
        % Last but not least, push the new compressed stuff into the stack,
        % and add an offset entry for it
        cs = push(cs, Z(counter+1:end, :));
		cs1 = size(Z, 2);
		
        col_offsets = [ col_offsets , 0 ];
        row_offsets = [ row_offsets , 0 ];
        
        H.D = A(mh+1:mh+m, nh+1:nh+n);
	else
		% fprintf('> Non-leafnode, mh = %d, nh = %d, m = %d, n = %d\n', mh, nh, m, n);
        % Call the constructor recursively on the left and right childs 
        [H.hssl, rs, cs, row_offsets, col_offsets] = BuildHSS_iter(...
            H.hssl, A, tol, rs, cs, ...
            row_offsets, col_offsets, mh, nh);
        
        [H.hssr, rs, cs, row_offsets, col_offsets] = BuildHSS_iter(...
            H.hssr, A, tol, rs, cs, row_offsets, col_offsets, ...
            mh + size(H.hssl, 1), nh + size(H.hssl, 2));
		        
        % Extract Bl and Bu from the stacks, and merge the children
        [row1, rs] = pop(rs); [row2, rs] = pop(rs);
        [col1, cs] = pop(cs); [col2, cs] = pop(cs);
		
		rs1 = size(row2, 1) - size(row1, 1);
		cs1 = size(col2, 1) - size(col1, 1);
		
		[rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H);
		
		assert(rs1 == rrV)
		assert(cs1 == rrU)
		
		H.Bu = row2(1:rrV, 1:rlU)';
        H.Bl = col2(1:rrU, 1:rlV);
		
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
			Z = elementAt(cs, j);

			% FIXME: We should preallocate B
			B = [ Z(col_offsets(j+1)+1:col_offsets(j+1)+rU,:) , B ];
		end
		
		[U, Z] = compress_hss_block(B, tol);
		
		H.Rl = U(1:rlU, :);
		H.Rr = U(rlU+1:end, :);
		
		% Push Z in the stack, and updates the elements in the as well
        counter = 0;
		for j = 0 : length(col_offsets) - 1
			W = elementAt(cs, j);

			cs = setElementAt(cs, [ W(1:col_offsets(j+1), :) ; ...
				Z(counter+1:counter+size(W, 2), :)' ; ...
				W(col_offsets(j+1)+rU+1:end,:) ], j);

			counter = counter + size(W, 2);
			col_offsets(j+1) = col_offsets(j+1) + size(Z, 2);
		end
		
		rs = push(rs, Z(counter+1:end, :));
        
        B = [ col2(rrU+1:end, :), col1 ]';
		
		for j = length(row_offsets) - 1 : -1 : 0
			Z = elementAt(rs, j);

			% FIXME: We should preallocate B
			B = [ Z(row_offsets(j+1)+1:row_offsets(j+1)+rV,:) , B ];
		end
		
		[U, Z] = compress_hss_block(B, tol);
		
		H.Wl = U(1:rlV, :);
		H.Wr = U(rlV+1:end, :);
		
		% Push Z in the stack, and updates the elements in the as well
        counter = 0;
		for j = 0 : length(row_offsets) - 1
			W = elementAt(rs, j);

			rs = setElementAt(rs, [ W(1:row_offsets(j+1), :) ; ...
				Z(counter+1:counter+size(W, 2), :)' ; ...
				W(row_offsets(j+1)+rV+1:end,:) ], j);

			counter = counter + size(W, 2);
			row_offsets(j+1) = row_offsets(j+1) + size(Z, 2);
		end
		
		cs = push(cs, Z(counter+1:end, :));
        
		row_offsets = [ row_offsets, 0 ];
        col_offsets = [ col_offsets, 0 ];
    end
end

function [rU, rV, rlU, rlV, rrU, rrV] = generatorsRank(H)
	if H.hssl.leafnode
		rlU = size(H.hssl.U, 2);
		rlV = size(H.hssl.V, 2);
	else
		rlU = size(H.hssl.Rl, 2);
		rlV = size(H.hssl.Wl, 2);
	end
			
	if H.hssr.leafnode
		rrU = size(H.hssr.U, 2);
		rrV = size(H.hssr.V, 2);
	else
		rrU = size(H.hssr.Rl, 2);
		rrV = size(H.hssr.Wl, 2);
	end
	
	rU = rlU + rrU;
	rV = rlV + rrV;
end

function H = build_hss_tree(m, n, block_size)
    H = hss();
    
    if m ~= n
        error('Rectangular HSS matrices are not supported');
    end   

    H.topnode  = 0;
    
    if n > block_size
        [m1, m2] = split_indices(m);
        [n1, n2] = split_indices(n);
        
        H.ml = m1; H.mr = m2; 
        H.nl = n1; H.nr = n2;
        
        H.hssl = build_hss_tree(m1, n1, block_size);
        H.hssr = build_hss_tree(m2, n2, block_size);
        
        H.leafnode = 0;
    else
        H.leafnode = 1;
        H.D = zeros(m, n);
    end
end

function [n1, n2] = split_indices(n)
    n1 = ceil(n / 2);
    n2 = n - n1;
end

function [U, Z] = compress_hss_block(B, tol)
	if isempty(B)
		U = zeros(size(B, 1), 0);
		Z = zeros(0, size(B, 2));
		return
	end


   [U, B, Z] = svd(B, 'econ');
   
   if tol == 0
       rk = size(B, 1);
   else   
       rk = sum(diag(B) > tol * B(1,1));
   end
   
   U = U(:, 1:rk);
   B = B(1:rk, 1:rk);
   
   Z = Z(:, 1:rk) * B;
end

function s = createStack()
	s = {};
end

function s = push(s, el)
	s{end+1} = el;
end

function [el, s] = pop(s)
	el = s{end};
	s = { s{1:end-1} };
end

function el = elementAt(s, j)
	el = s{j+1};
end

function s = setElementAt(s, el, j)
	s{j+1} = el;
end


