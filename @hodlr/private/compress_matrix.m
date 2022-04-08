function [U,V] = compress_matrix(A)
%COMPRESS_MATRIX Compress a dense matrix into low rank format.

threshold = hodlroption('threshold');
compression = hodlroption('compression');

if max(size(A)) < 256
	compression = 'svd';
elseif issparse(A)
	compression = 'lanczos';
end

switch compression
	case 'svd'
		[U,S,V] = svd(full(A), 'econ');

        switch hodlroption('norm')
            case 2
                k = sum(abs(diag(S)) > S(1,1) * threshold);
            case 'fro'
                t = diag(S);
                tt = sqrt(cumsum(t(end:-1:1).^2));
                k = sum(tt > tt(end) * threshold);
        end                

		U = U(:,1:k) * S(1:k,1:k);
		V = V(:,1:k);	
	case 'qr'
		[U, V] = prrqr(A, threshold);
		V = V';
	case 'lanczos'
		[U, S, V] = lanczos_svd(@(x,t) sparse_matvec(A, x, t), ...
			size(A, 1), size(A, 2), threshold);
		U = U * S;
	otherwise
		error('Unsupported compression method selected');
	end
end

function y = sparse_matvec(A, x, t)
	if strcmp(t, 'transp')
		y = A' * x;
	else
		y = A * x;
	end
end

