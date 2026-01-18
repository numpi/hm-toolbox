function A = full(T)
    if isequal(class(T),'tel')
        if T.top %At the deepest level T.M is a dense matrix
            A=blkdiag(T.U{:})*T.M*blkdiag(T.V{:})'+blkdiag(T.D{:});
        else
            A=blkdiag(T.U{:})*full(T.M)*blkdiag(T.V{:})'+blkdiag(T.D{:});
        end
    else 
    error('Input must be an hm_toolbox @hss object.')
    end
end

