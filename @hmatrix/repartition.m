function H = repartition(H, maxrank)
%REPARTITION
if ~exist('maxrank', 'var')
    maxrank = ceil(min(size(H)) / 10);
end

H = hmatrix_repartition(H, maxrank);
H = hmatrix_compress(H);

end

