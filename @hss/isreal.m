function b = isreal(H)
if H.leafnode == 1
    b = isreal(H.D) && isreal(H.U) && isreal(H.V);
else
    b = isreal(H.A11) && isreal(H.A22) && isreal(H.B21) && isreal(H.B12) && isreal(H.Rl) && isreal(H.Rr) && isreal(H.Wr) && isreal(H.Wl);
end
end
