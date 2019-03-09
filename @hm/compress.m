function H = compress(H)
%COMPRESS Recompress the HM representation


H = compress_hmatrix(H, norm(H));

end
