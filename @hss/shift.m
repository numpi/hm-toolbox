function H = shift(H, s)
%SHIFT Compute an HSS representation of H = H + sI (without recompression)

if is_leafnode(H)
    H.D = H.D + s*eye(size(H));
else
    H.A11 = shift(H.A11, s);
    H.A22 = shift(H.A22, s);
end

end

