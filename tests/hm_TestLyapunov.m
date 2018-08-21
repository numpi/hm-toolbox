function hm_TestLyapunov

for bs = [ 32, 256 ]

    hmoption('block-size', bs);

    % Simpler problem: 2D Laplacian discretized with finite difference, and
    % non-smooth RHS with singularity on the diagonal
    n = 1024; h = 1 / (n - 1);
    sA = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
    A  = hm('banded', sA, 1);
    C = h^2 * hm('chebfun2', @(x,y) log(1 + abs(x - y)), [0, 1], [0, 1], n);

    % Solver 1: completely HM
    X = lyap(A, C);
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        2 * norm(A) * norm(X) * hmoption('threshold'), ...
        [ 'Lyapunov equation with HM coefficients ' ...
        sprintf('and RHS (block-size = %d)', bs) ]);
    
    % Solver 1: HM + sparsity
    X = lyap(A, C, sA);
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        2 * norm(A) * norm(X) * hmoption('threshold'), ...
        [ 'Lyapunov equation with HM coefficients ' ...
        sprintf('and RHS (exploiting sparsity, block-size = %d)', bs) ]);    
    
    % Solver 2: Sign function
    X = lyap(A, C, 'method', 'sign');
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        10 * 2 * norm(A) * norm(X) * hmoption('threshold'), ...
        [ 'Lyapunov equation with HM coefficients ' ...
        sprintf('and RHS (sign function, block-size = %d)', bs) ]);

    % Modify the problem to be unsymmetric, and test the Sylvester solvers.
    % Here we test a non-constant coefficients PDE
    g = @(x) (1 + x).^2; d = g(linspace(0, 1, n)); D = hm('diagonal', d);
    B = D * A; sB = spdiags(d.', 0, n, n) * sA;
    X = lyap(A, B, C);
    CheckTestResult(norm(A * X + X * B + C), '<', ...
     (norm(B) + norm(A)) * norm(X) * hmoption('threshold'), ...
     [ 'Sylvester equation with HM coefficients and RHS' ...
     sprintf(' (block-size = %d)', bs) ]);

    % Same problem, exploiting sparsity
    X = lyap(A, B, C, sA, sB);
    CheckTestResult(norm(A * X + X * B + C), '<', ...
     (norm(B) + norm(A)) * norm(X) * hssoption('threshold'), ...
     [ 'Sylvester equation with HSS coefficients and RHS' ...
     sprintf(' (sparse arithmetic used for EK, block-size = %d)', bs) ]);
 
    % Same problem, exploiting sparsity
    X = lyap(A, B, C, 'method', 'sign');
    CheckTestResult(norm(A * X + X * B + C), '<', ...
     10 * (norm(B) + norm(A)) * norm(X) * hssoption('threshold'), ...
     [ 'Sylvester equation with HSS coefficients and RHS' ...
     sprintf(' (sign function, block-size = %d)', bs) ]);
 
end

end
