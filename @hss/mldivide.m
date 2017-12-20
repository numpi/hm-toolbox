function X = mldivide(A,B)
	if isa(B,'hss')
		error('unsupported');
	else
		if size(A, 1) <= hssoption('block-size')
			X = full(A) \ B;
			return;
		end
		[D,U,R,BB,W,V,tr] = hss2xia(A);
		X = hssulvsol(tr,D,U,R,BB,W,V,length(tr),B);
		% X = hss_solve(A, B);
	end
end
