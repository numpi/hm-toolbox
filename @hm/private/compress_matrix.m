function [U, V] = compress_matrix(Uold, Vold)
%COMPRESS_MATRIX Get a low-rank representation for Uold * Vold'. 

threshold = eps;

if ~exist('Vold', 'var')
	[U,S,V] = svd(Uold);
	
	rk = sum(diag(S) > threshold);
	U = U(:,1:rk) * sqrt(S(1:rk,1:rk));
	V = V(:,1:rk) * sqrt(S(1:rk,1:rk));
else
	[QU,RU] = qr(Uold, 0);
	[QV,RV] = qr(Vold, 0);
	
	[U,V] = compress_matrix(RU * RV');
	U = QU * U;
	V = QV * V;
end


end

