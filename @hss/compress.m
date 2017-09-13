function H = compress(H)
%COMPRESS Recompress the HSS representation

H = hss_compress(H, hssoption('threshold'));

end

