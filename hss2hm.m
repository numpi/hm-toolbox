function B = hss2hm(A)
%HSS2HM Conversion from HSS to HM format. 
%
% B = HSS2HM(A) constructs an HM representation of the HSS matrix A. 

B = hss2hm_rec(A, hm());

end

function [B, U, V] = hss2hm_rec(A, B)
    if A.leafnode == 1
        B.F = A.D;
        B.sz = size(A);
        
        U = A.U;
        V = A.V;
    else
        [B.A11, U1, V1] = hss2hm_rec(A.A11, hm());
        [B.A22, U2, V2] = hss2hm_rec(A.A22, hm());
        
        B.U12 = U1 * A.B12;
        B.V12 = V2;
        
        B.U21 = U2 * A.B21;
        B.V21 = V1;
        
        B.sz = size(A);
        
        if ~A.topnode         
            U = [ U1 * A.Rl ; U2 * A.Rr ];
            V = [ V1 * A.Wl ; V2 * A.Wr ];
        end
    end
end

