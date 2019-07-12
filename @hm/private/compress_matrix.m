function [U,V] = compress_matrix(A)
%COMPRESS_MATRIX Compress a dense matrix into low rank format.

threshold = hmoption('threshold');
compression = hmoption('dense-compression');

if ~strcmp(compression, 'svd') && ~strcmp(compression, 'rrqr')
	error('Invalid value for the dense-compression option');
end

if max(size(A)) < 256 || strcmp(compression, 'svd')
    [U,S,V] = svd(A);

	k = sum(abs(diag(S)) > S(1,1) * threshold);

	U = U(:,1:k) * S(1:k,1:k);
	V = V(:,1:k);	
else
	[U, V] = prrqr(A, threshold);
	V = V';
end


end

