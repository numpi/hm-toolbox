function H = mrdivide(H1, H2)

if isa(H2, 'hss')
	if isscalar(H1)
		H = inv(H2) * H1;
	else
		H = hss_ulv_solve(H2', H1')'; 
	end	     
else
    if isscalar(H2)
        H = hss_scalar_mul(1 / H2, H1);
    elseif isa(H2, 'hm')
		H = hss2hm(H1)/H2;
    else
		H = hss_ulv_solve(H2', H1')';     
    end
end
