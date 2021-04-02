function Ht = conj(H)
%CONJ Conjugate of HALR matrix H.

Ht = H;

if is_leafnode(H)
    if H.admissible
        H.U = conj(H.U);
        H.V = conj(H.V);
    else
        Ht.F = conj(H.F);
    end
else
    Ht.A11 = conj(H.A11);
    Ht.A12 = conj(H.A12);
    Ht.A21 = conj(H.A21);
    Ht.A22 = conj(H.A22);
end

end

