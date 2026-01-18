function z = tel_vec_mat_mul(T, w)
%Computes w'*T where T is in @tel
%INPUT: T in @tel format, w column vector
%OUTPUT: z = w'*T

    % Base case T.M is dense
    if ~isa(T.M, 'tel')
        % z = D'*w + V * (M' * (U'*w))
        u_w = apply_blocks(T.U, w, true); % U' * w
        inner = (T.M') * u_w; % M' * (U'*w)
        z = apply_blocks(T.D, w, true) + apply_blocks(T.V, inner, false); % D'*w + V*inner
    else
    % ---Recursion---
    u_w = apply_blocks(T.U, w, true);
    % Internal level (trasposed)
    inner_prod = tel_vec_mat_mul(T.M, u_w);
    %Ricostruction (z = D'*w + V*inner_prod)
    z = apply_blocks_(T.D, w, true) + apply_blocks_(T.V, inner_prod, false);
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