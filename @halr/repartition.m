function H = repartition(H, maxrank)
%REPARTITION

if ~exist('maxrank', 'var')
    maxrank = ceil(min(size(H)) / 10);
end

nrm = norm(H, halroption('norm'));

H = halr_repartition(H, maxrank, nrm);
H = halr_compress(H);

end

