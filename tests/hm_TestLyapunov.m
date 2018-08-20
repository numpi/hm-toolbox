function hm_TestLyapunov

for bs = [ 32, 256 ]

    hmoption('block-size', bs);

    % Simpler problem: 2D Laplacian discretized with finite difference, and
    % non-smooth RHS with singularity on the diagonal
    n = 4096; h = 1 / (n - 1);
    sA = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
    A  = hm('banded', sA, 1);
    C = h^2 * hm('chebfun2', @(x,y) log(1 + abs(x - y)), [0, 1], [0, 1], n);

    % Solver 1: completely HM
    X = lyap(A, -C);
    CheckTestResult(norm(A * X + X * A - C), '<', ...
        2 * norm(A) * norm(X) * hssoption('threshold'), ...
        [ 'Lyapunov equation with HSS coefficients ' ...
        sprintf('and RHS (block-size = %d)', bs) ]);

    % Sylvester equations are still unsupported by the HM interface -- tests
    % will be added here when they will become available. 

end

end