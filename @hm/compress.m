function H = compress(H, tol)
%COMPRESS Recompress the HM representation

if ~exist('tol', 'var')
    tol = hmoption('threshold');
end

H = compress_hmatrix(H, tol * norm(H));

end
