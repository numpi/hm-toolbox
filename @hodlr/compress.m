function H = compress(H)
%COMPRESS Recompress the HODLR representation


H = compress_hodlr(H, norm(H));

end
