function H = repartition(H, maxrank)
%REPARTITION

if ~exist('maxrank', 'var')
    maxrank = ceil(min(size(H)) / 10);
end

nrm = norm(H, hmatrixoption('norm'));

H = hmatrix_repartition(H, maxrank, nrm);
H = hmatrix_compress(H);

end

