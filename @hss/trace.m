function t = trace(H)
%TRACE computes the trace of an HSS matrix

if ~check_cluster_equality(H)
    error('trace is supported only for square matrices with square diagonal blocks');
end

if H.leafnode,
    t = trace(H.D);
else
    t = trace(H.A11) + trace(H.A22);
end
