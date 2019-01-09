function iH = imag(H)
%IMAG Imaginary part of an H-Matrix.
%
% iH = imag(H) extracts the imaginary part of an H-Matrix. Notice that in
%              general, the quasiseparable rank of rH can be twice as much
%              as the quasiseparable rank of H.

iH = H;

if is_leafnode(H)
    iH.F = imag(H.F);
else
    iH.A11 = imag(H.A11);
    iH.A22 = imag(H.A22);
    
    % Notice that the signs are caused by the fact that we use the
    % ctranspose operator in MATLAB.
    iH.U21 = [ real(H.U21) , imag(H.U21) ];
    iH.V21 = [ -imag(H.V21), real(H.V21) ];
    iH.U12 = [ real(H.U12) , imag(H.U12) ];
    iH.V12 = [ -imag(H.V12), real(H.V12) ];
end


end
