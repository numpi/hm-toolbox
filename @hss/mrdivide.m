function X = mrdivide(A, B)

if isa(B, 'hss')
    X = (B' \ A')';        
else
    if isscalar(B)
        X = hss_scalar_mul(1 / B, A);
        return;
    end
    
    X = (B' \ A')';
end

end
