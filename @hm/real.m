function rH = real(H)
%REAL Real part of an H-Matrix.
%
% rH = real(H) extracts the real part of an H-Matrix. Notice that in
%              general, the quasiseparable rank of rH can be twice as much
%              as the quasiseparable rank of H.

rH = H;

if is_leafnode(H)
    rH.F = real(H.F);
else
    rH.A11 = real(H.A11);
    rH.A22 = real(H.A22);
    
    % Notice that the signs are caused by the fact that we use the
    % ctranspose operator in MATLAB.
    rH.U21 = [ real(H.U21) , imag(H.U21) ];
    rH.V21 = [ real(H.V21) , imag(H.V21) ];
    rH.U12 = [ real(H.U12) , imag(H.U12) ];
    rH.V12 = [ real(H.V12) , imag(H.V12) ];
end


end

