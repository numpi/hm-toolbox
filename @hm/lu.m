function [HL, HU] = lu(H)
%LU computes LU factorization for HODLR matrices

[HL, HU] = hmatrix_lu(H);

HL = compress_hmatrix(HL);
HU = compress_hmatrix(HU);
