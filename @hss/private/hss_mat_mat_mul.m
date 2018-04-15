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

end

function g = bottom_up(A, B)
g = struct();

if (A.leafnode == 0)   % not a leaf
   
    g.gl = bottom_up(A.hssl, B.hssl);
    g.gr = bottom_up(A.hssr, B.hssr);
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
        fl = A.Bu * g.gr.VU * B.Bl;
        fr = A.Bl * g.gl.VU * B.Bu;
    else % not top node
        fl = A.Bu * g.gr.VU * B.Bl + A.Rl * f * B.Wl';
        fr = A.Bl * g.gl.VU * B.Bu + A.Rr * f * B.Wr';
    end
if A.topnode == 1
    C.Bl = blkdiag(A.Bl, B.Bl);
    C.Bu = blkdiag(A.Bu, B.Bu);
else
    C.Bl = [A.Bl A.Rr * f * B.Wl'; zeros(size(B.Bl,1),size(A.Bl,2)) B.Bl];
    C.Bu = [A.Bu A.Rl * f * B.Wr'; zeros(size(B.Bu,1), size(A.Bu,2)) B.Bu];
    C.Wl = [A.Wl zeros(size(A.Wl,1),size(B.Wl,2)); B.Bl' * g.gr.VU' * A.Wr B.Wl];
    C.Wr = [A.Wr zeros(size(A.Wr,1),size(B.Wr,2)); B.Bu' * g.gl.VU' * A.Wl B.Wr];
    C.Rl = [A.Rl A.Bu * g.gr.VU * B.Rr; zeros(size(B.Rl,1),size(A.Rl,2)) B.Rl];
    C.Rr = [A.Rr A.Bl * g.gl.VU * B.Rl; zeros(size(B.Rr,1),size(A.Rr,2)) B.Rr];
end
    C.hssl = up_bottom(A.hssl, B.hssl, g.gl, fl, C.hssl);
    C.hssr = up_bottom(A.hssr, B.hssr, g.gr, fr, C.hssr);
   
else        %leaf node
    C.D = A.D * B.D + A.U * f * B.V';
    C.U = [A.U, A.D * B.U];
    C.V = [B.D' * A.V, B.V];
end
end
