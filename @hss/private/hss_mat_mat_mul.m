function C = hss_mat_mat_mul(A, B)

if (A.leafnode==0)            % not a leaf
    g = bottom_up(A, B);
    f = 0;
    C = A;
    C = up_bottom(A, B, g, f, C);
else
    C = hss();
    C.leafnode = 1;
    C.topnode = 1;
    C.D = A.D * B.D;
end

C = hss_compress(C, hssoption('threshold'));

end

function g = bottom_up(A, B)
g = struct();

if (A.leafnode == 0)   % not a leaf
    
    g.gl = bottom_up(A.A11, B.A11);
    g.gr = bottom_up(A.A22, B.A22);
    if A.topnode == 0    % not root
        g.VU = A.Wl' * g.gl.VU * B.Rl + A.Wr' * g.gr.VU * B.Rr;
    end
else        % leaf
    g.VU = A.V' * B.U;
end
end


function C = up_bottom(A, B, g, f, C)

if (A.leafnode==0)   % not a leaf
    if A.topnode == 1 % root
        fl = A.B12 * g.gr.VU * B.B21;
        fr = A.B21 * g.gl.VU * B.B12;
    else % not top node
        fl = A.B12 * g.gr.VU * B.B21 + A.Rl * f * B.Wl';
        fr = A.B21 * g.gl.VU * B.B12 + A.Rr * f * B.Wr';
    end
    if A.topnode == 1
        C.B21 = blkdiag(A.B21, B.B21);
        C.B12 = blkdiag(A.B12, B.B12);
    else
        C.B21 = [A.B21 A.Rr * f * B.Wl'; zeros(size(B.B21,1),size(A.B21,2)) B.B21];
        C.B12 = [A.B12 A.Rl * f * B.Wr'; zeros(size(B.B12,1), size(A.B12,2)) B.B12];
        C.Wl = [A.Wl zeros(size(A.Wl,1),size(B.Wl,2)); B.B21' * g.gr.VU' * A.Wr B.Wl];
        C.Wr = [A.Wr zeros(size(A.Wr,1),size(B.Wr,2)); B.B12' * g.gl.VU' * A.Wl B.Wr];
        C.Rl = [A.Rl A.B12 * g.gr.VU * B.Rr; zeros(size(B.Rl,1),size(A.Rl,2)) B.Rl];
        C.Rr = [A.Rr A.B21 * g.gl.VU * B.Rl; zeros(size(B.Rr,1),size(A.Rr,2)) B.Rr];
    end
    C.A11 = up_bottom(A.A11, B.A11, g.gl, fl, C.A11);
    C.A22 = up_bottom(A.A22, B.A22, g.gr, fr, C.A22);
    
    C.ml = size(A.A11,1); C.nl = size(B.A11,2);
    C.mr = size(A.A22,1); C.nr = size(B.A22,2);
    
else        %leaf node
    C.D = A.D * B.D + A.U * f * B.V';
    C.U = [A.U, A.D * B.U];
    C.V = [B.D' * A.V, B.V];
    
    C.ml = []; C.nl = []; C.nr = []; C.mr = [];
end
end
