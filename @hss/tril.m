function tH = tril(H)
if ~check_cluster_equality(H)
    error('tril is supported only for square matrices with square diagonal blocks');
end
tH = tril_ric(H);
tH = compress(tH);
end

function tH = tril_ric(H)
tH = H;
if H.leafnode == 1
    tH.D = tril(H.D);
else
    if H.topnode == 0
        tH.B12 = zeros(size(tH.Rl, 1), size(tH.Wr, 1));
    elseif H.A11.leafnode == 1
        tH.B12 = zeros(size(tH.A11.U, 2), size(tH.A22.V, 2));
    else
        tH.B12 = zeros(size(tH.A11.Rl, 2), size(tH.A22.Wr, 2));
    end
    tH.A11 = tril_ric(H.A11);
    tH.A22 = tril_ric(H.A22);
end
end
