function Ht = transpose(H)
%TRANSPOSE Transpose of the H-matrix H.

Ht = H;

if ~isempty(H.F)
    Ht.F = H.F.';
else
    Ht.A11 = H.A11.';
    Ht.A22 = H.A22.';
    Ht.U12 = conj(H.V21);
    Ht.V12 = conj(H.U21);
    Ht.U21 = conj(H.V12);
    Ht.V21 = conj(H.U12);
end

end

