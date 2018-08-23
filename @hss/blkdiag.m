function A = blkdiag(A1, A2)
A = hss();

A.A11 = A1;
A.A22 = A2;

A.A11.topnode = 0;
A.A22.topnode = 0;

A.topnode = 1;
A.leafnode = 0;

A.ml = size(A1, 1);
A.nl = size(A1, 2);
A.mr = size(A2, 1);
A.nr = size(A2, 2);

A.B21 = [];
A.B12 = [];

if A.A11.leafnode == 1
    A.A11.U = zeros(A.ml, 0);
    A.A11.V = zeros(A.nl, 0);
else
    %A.A11.Rl = zeros(hss_generator_size(A.A11, 'left'), 0);
    %A.A11.Rr = zeros(hss_generator_size(A.A11, 'left'), 0);
    %A.A11.Wl = zeros(hss_generator_size(A.A11, 'right'), 0);
    %A.A11.Wr = zeros(hss_generator_size(A.A11, 'right'), 0);
    
    A.A11.Rl = zeros(size(A.A11.B12,1), 0);
    A.A11.Rr = zeros(size(A.A11.B21,1), 0);
    A.A11.Wl = zeros(size(A.A11.B21,2), 0);
    A.A11.Wr = zeros(size(A.A11.B12,2), 0);
end

if A.A22.leafnode == 1
    A.A22.U = zeros(A.mr, 0);
    A.A22.V = zeros(A.nr, 0);
else
    %A.A22.Rl = zeros(hss_generator_size(A.A22, 'left'), 0);
    %A.A22.Rr = zeros(hss_generator_size(A.A22, 'left'), 0);
    %A.A22.Wl = zeros(hss_generator_size(A.A22, 'right'), 0);
    %A.A22.Wr = zeros(hss_generator_size(A.A22, 'right'), 0);
    
    A.A22.Rl = zeros(size(A.A22.B12,1), 0);
    A.A22.Rr = zeros(size(A.A22.B21,1), 0);
    A.A22.Wl = zeros(size(A.A22.B21,2), 0);
    A.A22.Wr = zeros(size(A.A22.B12,2), 0);
end


end


