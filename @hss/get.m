function M = get(H, II, JJ)
% Equivalent to full(H(m,n))

[Is, Ip] = sort(II); 
[Js, Jp] = sort(JJ);
iIp(Ip) = [1:length(II)];
iJp(Jp) = [1:length(JJ)];

M = full(hss_sub(H, Is, Js));
M = M(iIp, iJp);
end

function H = hss_sub(H, mind, nind);

[m,n] = size(H);

if H.leafnode
    H.ml = []; H.nl = []; H.mr = []; H.nr = []; 
    H.D = H.D(mind, nind);
    H.U = H.U(mind, :);
    H.V = H.V(nind, :);
else
    [m1,n1] = size(H.A11);
    mind1 = intersect( 1:m1, mind );
    nind1 = intersect( 1:n1, nind );
    H.A11 = hss_sub(H.A11, mind1, nind1);
    H.ml = length(mind1); H.nl = length(nind1);
    
    [m2,n2] = size(H.A22);
    mind = mind - m1; nind = nind - n1;
    mind2 = intersect( 1:m2, mind );
    nind2 = intersect( 1:n2, nind );
    H.A22 = hss_sub(H.A22, mind2, nind2);
    H.mr = length(mind2); H.nr = length(nind2);
end

end

