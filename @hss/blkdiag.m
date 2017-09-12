function A = blkdiag(A1, A2)
	A = hss();

	A.hssl = A1;
	A.hssr = A2;

	A.hssl.topnode = 0;
	A.hssr.topnode = 0;

	A.topnode = 1;
	A.leafnode = 0;

	A.ml = size(A1, 1);
	A.nl = size(A1, 2);
	A.mr = size(A2, 1);
	A.nr = size(A2, 2);

	A.Bl = [];
	A.Bu = [];

	if A.hssl.leafnode == 1
		A.hssl.U = zeros(A.ml, 0);
		A.hssl.V = zeros(A.nl, 0);
    else
		%A.hssl.Rl = zeros(hss_generator_size(A.hssl, 'left'), 0);
		%A.hssl.Rr = zeros(hss_generator_size(A.hssl, 'left'), 0);
		%A.hssl.Wl = zeros(hss_generator_size(A.hssl, 'right'), 0);
		%A.hssl.Wr = zeros(hss_generator_size(A.hssl, 'right'), 0);
		
		A.hssl.Rl = zeros(size(A.hssl.Bu,1), 0);
		A.hssl.Rr = zeros(size(A.hssl.Bl,1), 0);
		A.hssl.Wl = zeros(size(A.hssl.Bl,2), 0);
		A.hssl.Wr = zeros(size(A.hssl.Bu,2), 0);
	end

	if A.hssr.leafnode == 1
		A.hssr.U = zeros(A.mr, 0);
		A.hssr.V = zeros(A.nr, 0);
	else
		%A.hssr.Rl = zeros(hss_generator_size(A.hssr, 'left'), 0);
		%A.hssr.Rr = zeros(hss_generator_size(A.hssr, 'left'), 0);
		%A.hssr.Wl = zeros(hss_generator_size(A.hssr, 'right'), 0);
		%A.hssr.Wr = zeros(hss_generator_size(A.hssr, 'right'), 0);

		A.hssr.Rl = zeros(size(A.hssr.Bu,1), 0);
		A.hssr.Rr = zeros(size(A.hssr.Bl,1), 0);
		A.hssr.Wl = zeros(size(A.hssr.Bl,2), 0);
		A.hssr.Wr = zeros(size(A.hssr.Bu,2), 0);
    end


end


