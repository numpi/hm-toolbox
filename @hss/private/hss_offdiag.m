function [U,V] = hss_offdiag(A, ul)

[Ul,Vl] = hss_generators(A.A11);
[Ur,Vr] = hss_generators(A.A22);

if ~exist('ul', 'var')
    ul = 'all';
end

switch ul
    case 'upper'
        U = Ul * A.B12;
        V = Vr;
    case 'lower'
        U = Ur * A.B21;
        V = Vl;
    case 'all'
        U = [ zeros(size(Ul,1), size(A.B21,2)) , ...
            Ul * A.B12 ; Ur * A.B21 , zeros(size(Ur,1), size(A.B12,2)) ];
        V = [ Vl , zeros(size(Vl,1), size(Vr,2)) ; ...
            zeros(size(Vr,1), size(Vl,2)) , Vr ];
end

[U, V] = compress_factors(U, V, hssoption('threshold'));

end

function [U,V] = hss_generators(A)
if A.leafnode == 1
    U = A.U;
    V = A.V;
else
    [Ul,Vl] = hss_generators(A.A11);
    [Ur,Vr] = hss_generators(A.A22);
    U = [ Ul * A.Rl ; Ur * A.Rr ];
    V = [ Vl * A.Wl ; Vr * A.Wr ];
end
end
