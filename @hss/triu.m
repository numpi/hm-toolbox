function tH = triu(H)
if ~check_cluster_equality(H)
    error('triu is supported only for square matrices with square diagonal blocks');
end
tH = triu_ric(H);
tH = compress(tH);
end

function tH = triu_ric(H)
tH = H;
if H.leafnode == 1
    tH.D = triu(H.D);
else
    if H.topnode == 0
    	tH.B21 = zeros(size(tH.Rr, 1), size(tH.Wl, 1));
    elseif H.A11.leafnode == 1
	tH.B21 = zeros(size(tH.A22.U, 2), size(tH.A11.V, 2));
    else
	tH.B21 = zeros(size(tH.A22.Rr, 2), size(tH.A11.Wl, 2));
    end
    tH.A11 = triu_ric(H.A11);
    tH.A22 = triu_ric(H.A22);
end
end
