function X = mldivide(A, B)

if isa(B,'hss')
    if ~isa(A, 'hss') && isscalar(A)
        X = hss_scalar_mul(1 / A, B);
    else
        X = hss_mldivide(A, B);
    end
else
    if isa(A, 'hss')
        if size(A, 1) <= hssoption('block-size')
            X = full(A) \ B;
            return;
        end
        
        X = hss_ulv_solve(A,B);        
    else
        if isscalar(A)
            X = hss_scalar_mul(1 / A, B);
        end
    end

end
end
