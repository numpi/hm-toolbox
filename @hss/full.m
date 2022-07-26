function [M, U, V] = full(H)
    if is_leafnode(H)
        M = H.D;
        U = H.U;
        V = H.V;
    else
        [M11, U1, V1] = full(H.A11);
        [M22, U2, V2] = full(H.A22);
        
        M = [ M11 , U1 * H.B12 * V2' ; U2 * H.B21 * V1', M22 ];
        
        if ~H.topnode
            U = blkdiag(U1, U2) * [ H.Rl ; H.Rr ];
            V = blkdiag(V1, V2) * [ H.Wl ; H.Wr ];
        else
            U = []; V = [];
        end
    end
end