function H = compress(H, tol)
%COMPRESS Recompress the HSS representation

if ~exist('tol', 'var')
    tol = hssoption('threshold');
end

H = hss_compress(H, tol);

end

