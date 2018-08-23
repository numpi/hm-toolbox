function B = hss_scalar_mul(a, B)
if B.leafnode == 1
    B.D = a * B.D;
else
    B.B12 = a * B.B12;
    B.B21 = a * B.B21;
    B.A11 = hss_scalar_mul(a, B.A11);
    B.A22 = hss_scalar_mul(a, B.A22);
end
end
