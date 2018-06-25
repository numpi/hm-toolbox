function y = ulv_solve_wrapper(nu, mu, ULV, x)

if nu > mu
    y = nu \ ulv_solve(ULV, x);
else
    y = -mu \ x;
end

