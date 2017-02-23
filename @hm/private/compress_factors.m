function [ U, V ] = compress_factors(Uold, Vold)
%COMPRESS_FACTORS Compress a low-rank representation U*V'. 

threshold = 1e-12;

if size(Uold, 2) < 12
	U = Uold;
	V = Vold;
else
	[QU, RU] = qr(Uold, 0);
	[QV, RV] = qr(Vold, 0);

	[U,S,V] = svd(RU * RV');
	rk = sum(diag(S) > S(1,1) * threshold);
	U = QU * U(:,1:rk) * S(1:rk,1:rk);
	V = QV * V(:,1:rk);
end

end

