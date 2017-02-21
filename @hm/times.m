function H = times(H1, H2)
%TIMES Scalar multiplication

if (isfloat(H1) && isscalar(H1)) || (isfloat(H2) && isscalar(H2))
	H = H1 * H2;
else
	error('Unsupported operation');
end

end

