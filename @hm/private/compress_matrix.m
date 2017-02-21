function [U, V] = compress_matrix(Uold, Vold)
%COMPRESS_MATRIX Get a low-rank representation for Uold * Vold'. 

threshold = eps;

if ~exist('Vold', 'var')
	[U,S,V] = svd(Uold);
	
	rk = sum(diag(S) > threshold);
	U = U(:,1:rk) * sqrt(S(1:rk,1:rk));
	V = V(:,1:rk) * sqrt(S(1:rk,1:rk));
else
	if size(Uold, 2) < 16
		U = Uold;
		V = Vold;
	else
		[QU,RU] = qr(Uold, 0);
		[QV,RV] = qr(Vold, 0);

		[U,S,V] = svd(RU * RV');
		rk = sum(diag(S) > threshold);
		U = QU * U(:,1:rk) * S(1:rk,1:rk);
		V = QV * V(:,1:rk);
	end
end

end

