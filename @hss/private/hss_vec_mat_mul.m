function x = hss_vec_mat_mul(v, A)
%X=HSS_MAT_VEC_MUL computes the fast Multiplication of a HSS matrix with a
%                  (block) vector x from right.
%                  The algorithodlr is taken from the paper:
%
% 		   S. Chandrasekaran, P. Dewilde, M. Gu, W. Lyons, and T. Pals, A fast
%                  solver for HSS representations via sparse matrices , SIAM J. Matrix Anal.
%                  Appl. 29 (2006/07), no. 1, 67--81 (electronic).

if (A.leafnode==0)            % not a leaf
    g = bottom_up(A, v);
    f = 0;
    x = up_bottom(A, v, g, f);
else
    x = v * A.D;
end


end

function g = bottom_up(A, v)
g = struct();

if (A.leafnode == 0)   % not a leaf
    
    g.gl = bottom_up(A.A11, v(:, 1:A.ml));
    g.gr = bottom_up(A.A22, v(:, A.ml+1:end));
    if A.topnode == 0    % not root
        g.vU = g.gl.vU * A.Rl + g.gr.vU * A.Rr;
    end
else        % leaf
    g.vU = v * A.U;
end

end

function x = up_bottom(A, v, g, f)

if (A.leafnode==0)   % not a leaf
    if A.topnode == 1 % root
        fl = g.gr.vU * A.B21;
        fr = g.gl.vU * A.B12;
    else % not top node
        fl = g.gr.vU * A.B21 + f * A.Wl';
        fr = g.gl.vU * A.B12 + f * A.Wr';
    end
    
    xl = up_bottom(A.A11, v(:, 1:A.ml), g.gl, fl);
    xr = up_bottom(A.A22, v(:, A.ml+1:end), g.gr, fr);
    x = [xl, xr];
    
else        %leaf node
    x = v * A.D + f * A.V';
end
end
