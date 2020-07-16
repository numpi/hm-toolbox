function t = trace(H)
%TRACE computes the trace of a HMATRIX matrix


if is_leafnode(H)
    if H.sz(1) ~= H.sz(2)
        error('TRACE:: non square diagonal block')
    end
    if H.admissible
        t = trace(H.V' * H.U);
    else
        t = trace(H.F);
    end
else
    t = trace(H.A11) + trace(H.A22);
end
end
