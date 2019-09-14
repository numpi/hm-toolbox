function H = mldivide(H1, H2)

if isa(H2,'hss')
    if isscalar(H1) && all(size(H1) == 1)
        H = hss_scalar_mul(1 / H1, H2);
    elseif isa(H1, 'hss')
        H = hss_mldivide(H1, H2);
    else
        H = H1\full(H2);
    end
else
    if size(H1, 1) <= hssoption('block-size') % We might remove this case, not sure is needed
            H = full(H1) \ H2;
            return;
    end    
    if isscalar(H2) && all(size(H2) == 1)
        H = inv(H1) * H2;
    elseif isa(H2, 'hodlr')
        H = hss2hodlr(H1)\H2;  
    else   
        H = hss_ulv_solve(H1, H2);     
    end

end
end



