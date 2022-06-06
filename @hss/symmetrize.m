function A = symmetrize(A, ud)
%
	if ~exist('ud', 'var')
		ud = 'down';
	end
A = hss_force_symmetry(A, ud);

end
