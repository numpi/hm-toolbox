function [HL, HU] = lu(H)
%LU computes LU factorization for HODLR matrices

[HL, HU] = hodlr_lu(H);

% HL = compress_hodlr(HL);
% HU = compress_hodlr(HU);
