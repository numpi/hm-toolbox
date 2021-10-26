function leaf = is_leafnode(B)
%IS_LEAFNODE Returns true if the HODLR B has no children

leaf = isempty(B.A11); %&& isempty(B.A22);

end

