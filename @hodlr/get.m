function M = get(H, I, J)
% Equivalent to full(H(m,n))

[Is, Ip] = sort(I); 
[Js, Jp] = sort(J);
iIp(Ip) = [1:length(I)];
iJp(Jp) = [1:length(J)];
M = full(hodlr_sub(H, Is, Js));
M = M(iIp, iJp);
end

function H = hodlr_sub(H, mind, nind);

[m,n] = size(H);

H.sz = [length(mind), length(nind)];
if is_leafnode(H)
    H.F = H.F(mind,nind);
else
    [m1,n1] = size(H.A11);
    mind1 = intersect( 1:m1, mind );
    nind1 = intersect( 1:n1, nind );
    H.A11 = hodlr_sub(H.A11, mind1, nind1);
    H.U12 = H.U12(mind1,:);
    H.V21 = H.V21(nind1,:);
    
    [m2,n2] = size(H.A22);
    mind = mind - m1; nind = nind - n1;
    mind2 = intersect( 1:m2, mind );
    nind2 = intersect( 1:n2, nind );
    H.A22 = hodlr_sub(H.A22, mind2, nind2);
    H.U21 = H.U21(mind2,:);
    H.V12 = H.V12(nind2,:);
end

end
