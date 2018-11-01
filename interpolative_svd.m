function [Xcol, Xrow, Jcol, Jrow] = interpolative_svd(A, tol)
	% compute interpolative low rank approximation of A
	if ~exist('tol', 'var')
		tol = hssoption('threshold');
	end
	[Qcol, Scol] = random_range_finder(A, size(A, 2), tol);
	[Qrow, Srow] = random_range_finder(A', size(A, 1), tol);
	[Xcol, Jcol] = interpolative(Qcol'); Xcol = Xcol';
	[Xrow, Jrow] = interpolative(Qrow'); Xrow = Xrow';
end
