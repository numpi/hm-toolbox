function T = hss_trace(A)
if A.leafnode == 1
    T = trace(A.D);
else
    T = hss_trace(A.A11) + hss_trace(A.A22);
end
end
