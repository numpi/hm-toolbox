function Y = solve_upper_triangular2(H, x)
% Compute X * H^(-1) with H upper triangular @hm object.

% Only dense matrices are supported as of now
if ~isempty(H.F)
    Y = x / H.F;
else
    mp = H.A11.sz(2);
    Y1 = solve_upper_triangular2(H.A11, x(:,1:mp));
    Y2 = solve_upper_triangular2(H.A22, x(:,mp+1:end) - (Y1 * H.U12) * H.V12');
    Y = [ Y1, Y2 ];
end

end

