function t = trace(H)
%TRACE computes the trace of a HODLR matrix

if ~check_cluster_equality(H)
    error('trace is supported only for square matrices with square diagonal blocks');
end

if is_leafnode(H)
    t = trace(H.F);
else
    t = trace(H.A11) + trace(H.A22);
end
