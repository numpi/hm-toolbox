function l = iszero(H)
% check whether the HALR H is a zero matrix
if is_leafnode(H)
    if H.admissible
        l = isempty(H.U);
    else
        l = all(all(H.F == 0));
    end
else
    l = iszero(H.A11) && iszero(H.A12) && iszero(H.A21) && iszero(H.A22);
end
end
