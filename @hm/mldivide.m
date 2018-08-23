function H = mldivide(H1, H2)
%MLDIVIDE solve linear systems with HODLR-matrices
if isa(H1,'hm')
    if isempty(H1.U21) % Upper triangular system
        H = solve_upper_triangular(H1,H2);
    elseif isempty(H1.U12) % Lower triangular system
        H = solve_lower_triangular(H1,H2);
    else
        [HL, HU] = lu(H1);
        H = HU \ (HL\H2);
    end
elseif isscalar(H1)
    H = H2 *(1/H1);
end

