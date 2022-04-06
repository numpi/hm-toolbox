function H = compress(H, nrm)
%COMPRESS Recompress the HODLR representation

if ~exist('nrm', 'var')
    H = halr_compress(H);
else
    H = halr_compress(H, nrm);
end

end
