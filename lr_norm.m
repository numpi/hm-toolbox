function nrm = lr_norm(U,V,nrmtype)
if size(U,2) ~= size(V,2)
	error('uncompatible dimensions');
end
[~, RU] = qr(full(U),0);
[~, RV] = qr(full(V),0);

if ~exist('nrmtype', 'var') || ischar(nrmtype) && strcmp(nrmtype, 'fro')
    nrm = norm(RU * RV', 'fro');
else
    nrm = norm(RU * RV', nrmtype);
end

end
