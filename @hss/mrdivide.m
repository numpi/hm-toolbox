function X = mrdivide(A, B)
    if isscalar(B)
        X = hss_scalar_mul(1 ./ B, A);
    else
        X = mldivide(B', A')';        
    end
        


end
