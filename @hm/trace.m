function t = trace(H)
%TRACE conputes the trace of a HODLR matrix
if ~isempty(H.F)
    t = trace(H.F);
else
    t = trace(H.A11) + trace(H.A22);
end
