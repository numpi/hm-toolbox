function leaf = is_leafnode(B)
%IS_LEAFNODE Returns true if the HALR B has no children

leaf = isempty(B.A11);

end

