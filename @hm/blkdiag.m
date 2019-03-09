function H = blkdiag(H1, H2)
H = hm();

H.A11 = H1;
H.A22 = H2;
s1 = size(H1);
s2 = size(H2);
H.sz = [s1(1) + s2(1), s2(1) + s2(2)];
H.U21 = zeros(s2(1), 0);
H.U12 = zeros(s1(1), 0);
H.V21 = zeros(s1(2), 0);
H.V12 = zeros(s2(2), 0);

end


