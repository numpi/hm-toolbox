function H = minus(H1,H2)
%MINUS difference of matrices
if ~isa(H1,'hm') || ~isa(H2,'hm')
    error('Unsupported operation');
else
    H = hmatrix_minus(H1, H2);
    H = compress_hmatrix(H);
end
