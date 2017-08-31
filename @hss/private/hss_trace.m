function T = hss_trace(A)
if A.leafnode == 1
	T = trace(A.D);
else
	T = hss_trace(A.hssl) + hss_trace(A.hssr);
end
end
