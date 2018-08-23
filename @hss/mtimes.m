function C = mtimes(A,B)
if isa(B,'hss')
    if isa(A,'hss')
        C = hss_mat_mat_mul(A, B);
    elseif isscalar(A)
        C = hss_scalar_mul(A,B);
    else
        C = hss_vec_mat_mul(A, B);
    end
else
    if isscalar(B)
        C = hss_scalar_mul(B,A);
    else
        C = hss_mat_vec_mul(A,B);
    end
end
end


