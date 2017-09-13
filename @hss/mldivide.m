function X = mldivide(A,B)
	if isa(B,'hss')
		error('unsupported');
	else
		[D,U,R,BB,W,V,tr] = hss2xia(A);
		X = hssulvsol(tr,D,U,R,BB,W,V,length(tr),B);
		%X = hss_solve(A, B);
	end
end
