function Ht = conj(H)
%CONJ Conjugate of HODLR matrix H.

Ht = H;

if ~isempty(H.F)
    Ht.F = conj(H.F);
else
    Ht.A11 = conj(H.A11);
    Ht.A22 = conj(H.A22);
    Ht.U12 = conj(H.U12);
    Ht.V12 = conj(H.V12);
    Ht.U21 = conj(H.U21);
    Ht.V21 = conj(H.V21);
end

end

