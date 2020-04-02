function x = hss_mat_vec_mul(A, v)
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
    x = A.D * v;
end


end

function g = bottom_up(A, v)
g = struct();

if (A.leafnode == 0)   % not a leaf    
    g.gl = bottom_up(A.A11, v(1:A.nl,:));
    g.gr = bottom_up(A.A22, v(A.nl+1:end,:));
    if A.topnode == 0    % not root
        g.Vv = A.Wl' * g.gl.Vv + A.Wr' * g.gr.Vv;
    end
else        % leaf
    g.Vv = A.V' * v;
end

end

function x = up_bottom(A, v, g, f)

if (A.leafnode==0)   % not a leaf
    if A.topnode == 1 % root
        fl = A.B12 * g.gr.Vv;
        fr = A.B21 * g.gl.Vv;
    else % not top node
        fl = A.B12 * g.gr.Vv + A.Rl * f;
        fr = A.B21 * g.gl.Vv + A.Rr * f;
    end
    
    xl = up_bottom(A.A11, v(1:A.nl,:), g.gl, fl);
    xr = up_bottom(A.A22, v(A.nl+1:end,:), g.gr, fr);
    x = [xl; xr];
    
else        %leaf node
    x = A.D * v + A.U * f;
end
end
