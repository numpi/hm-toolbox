function hodlr_TestLyapunov

for bs = [ 32, 256 ]
    
    hodlroption('block-size', bs);
    
    % Simpler problem: 2D Laplacian discretized with finite difference, and
    % non-smooth RHS with singularity on the diagonal
    n = 1024; h = 1 / (n - 1);
    sA = spdiags(ones(n, 1) * [1 -2 1], -1:1, n, n);
    A  = hodlr('banded', sA, 1);
    t = linspace(0, 1, n);
    f = @(i,j) log(1 + abs(t(i).' - t(j)));
    C = h^2 * hodlr('handle', f, n, n);
    
    % Solver 1: completely hodlr
    X = lyap(A, C);
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        20 * norm(A) * norm(X) * hodlroption('threshold'), ...
        [ 'Lyapunov equation with hodlr coefficients ' ...
        sprintf('and RHS (block-size = %d)', bs) ]);
    
    % Solver 1: hodlr + sparsity
    X = lyap(A, C, sA);
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        20 * norm(A) * norm(X) * hodlroption('threshold'), ...
        [ 'Lyapunov equation with hodlr coefficients ' ...
        sprintf('and RHS (exploiting sparsity, block-size = %d)', bs) ]);
    
    % Solver 2: Sign function
    X = lyap(A, C, 'method', 'sign');
    CheckTestResult(norm(A * X + X * A + C), '<', ...
        20 * 2 * norm(A) * norm(X) * hodlroption('threshold'), ...
        [ 'Lyapunov equation with hodlr coefficients ' ...
        sprintf('and RHS (sign function, block-size = %d)', bs) ]);
    
    % Modify the problem to be unsymmetric, and test the Sylvester solvers.
    % Here we test a non-constant coefficients PDE
    g = @(x) (1 + x).^2; d = g(linspace(0, 1, n)); D = hodlr('diagonal', d);
    B = D * A; sB = spdiags(d.', 0, n, n) * sA;
    X = lyap(A, B, C);
    CheckTestResult(norm(A * X + X * B + C), '<', ...
        20 * (norm(B) + norm(A)) * norm(X) * hodlroption('threshold'), ...
        [ 'Sylvester equation with hodlr coefficients and RHS' ...
        sprintf(' (block-size = %d)', bs) ]);
    
    % Same problem, exploiting sparsity
    X = lyap(A, B, C, sA, sB);
    CheckTestResult(norm(A * X + X * B + C), '<', ...
        20 * (norm(B) + norm(A)) * norm(X) * hodlroption('threshold'), ...
        [ 'Sylvester equation with hodlr coefficients and RHS' ...
        sprintf(' (sparse arithmetic used for EK, block-size = %d)', bs) ]);
    
    % Same problem, exploiting sparsity
    X = lyap(A, B, C, 'method', 'sign');
    CheckTestResult(norm(A * X + X * B + C), '<', ...
        20 * (norm(B) + norm(A)) * norm(X) * hodlroption('threshold'), ...
        [ 'Sylvester equation with hodlr coefficients and RHS' ...
        sprintf(' (sign function, block-size = %d)', bs) ]);
    
end

end
