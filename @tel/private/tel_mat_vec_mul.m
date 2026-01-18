function y = tel_mat_vec_mul(T, w)
%Computes T*w where T is in @tel
%INPUT: T in @tel format, w column vector
%OUTPUT: y = T*w

    % Base case T.M is dense
    if ~isa(T.M, 'tel')
        v_w = apply_blocks(T.V, w, true);
        y = apply_blocks(T.D, w, false) + apply_blocks(T.U, T.M * v_w, false);
    else
    % ---Recursion---
    v_w = apply_blocks(T.V, w, true);
    % Internal level
    inner_prod = tel_mat_vec_mul(T.M, v_w);
    % Riconstruction (D*w + U*inner_prod)
    y = apply_blocks(T.D, w, false) + apply_blocks(T.U, inner_prod, false);
    end
end

function out = apply_blocks(C, vec, is_hermitian)
    n_blocks = length(C);
    out_parts = cell(n_blocks, 1);
    curr_idx = 1;
    for i = 1:n_blocks
        [rows, cols] = size(C{i});
        if is_hermitian
            % Operation: C{i}' * v_part
            % Dimension: number of rows of C{i}
            sz_in = rows;
            v_part = vec(curr_idx : curr_idx + sz_in - 1);
            out_parts{i} = (C{i}') * v_part;
            curr_idx = curr_idx + sz_in;
        else
            % Operation: C{i} * v_part
            % Dimension: number of cols of C{i}
            sz_in = cols;
            v_part = vec(curr_idx : curr_idx + sz_in - 1);
            out_parts{i} = C{i} * v_part;
            curr_idx = curr_idx + sz_in;
        end
    end
    out = vertcat(out_parts{:});
end