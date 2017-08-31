function H = plus(H1, H2)
%PLUS Sum two H matrices. 

hmatrix_pack_plus(H1, H2)

H = hmatrix_plus(H1, H2);
H = compress_hmatrix(H);

end

