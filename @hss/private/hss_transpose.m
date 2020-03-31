function B = hss_transpose(A)
B = A;
if B.leafnode == 1;
    B.D = B.D.';
else
    B.ml = A.nl;
    B.nl = A.ml;
    B.mr = A.nr;
    B.nr = A.mr;
    B.B12 = A.B21.';
    B.B21 = A.B12.';
    if B.A11.leafnode == 0
        B.A11.Rl = conj(A.A11.Wl);
        B.A11.Rr = conj(A.A11.Wr);
        B.A11.Wl = conj(A.A11.Rl);
        B.A11.Wr = conj(A.A11.Rr);
    else
	B.A11.U = conj(A.A11.V);
        B.A11.V = conj(A.A11.U);
    end
    if B.A22.leafnode == 0
        B.A22.Rl = conj(A.A22.Wl);
        B.A22.Rr = conj(A.A22.Wr);
        B.A22.Wl = conj(A.A22.Rl);
        B.A22.Wr = conj(A.A22.Rr);
    else
        B.A22.U = conj(A.A22.V);
        B.A22.V = conj(A.A22.U);
    end
    B.A11 = hss_transpose(B.A11);
    B.A22 = hss_transpose(B.A22);
end
end
