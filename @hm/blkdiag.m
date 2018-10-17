function A = blkdiag(A1, A2)
A = hm();

A.A11 = A1;
A.A22 = A2;
s1 = size(A1);
s2 = size(A2);
A.sz = [s1(1) + s2(1), s2(1) + s2(2)];
A.U21 = zeros(s2(1), 0);
A.U12 = zeros(s1(1), 0);
A.V21 = zeros(s1(2), 0);
A.V12 = zeros(s2(2), 0);

end


