function H = compress(H)
%COMPRESS Recompress the HODLR representation

Ht = H';
nrm = normest_Afun(@(x) H*x, @(x) Ht*x, size(H, 2), 1e-3);
H = compress_hmatrix(H, nrm);

end
