function t = trace(H)
%TRACE computes the trace of a HODLR matrix

if size(H, 1) ~= size(H, 2)
    error('Trace is supported only for square matrices');
end

if ~isempty(H.F)
    t = trace(H.F);
else
    t = trace(H.A11) + trace(H.A22);
end
