function [U, V] = chebfun2_low_rank(f, xdom, ydom, m, n)
%CHEBFUN2_LOW_RANK Create a low rank representation of the sampling of F.

threshold = hodlroption('threshold');

try
    % Chebapprox is a simplified chebfun drop-in replacement that only
    % performs the interpolation part --- and that we can use in Octave
    % where chebfun does not work. 
    cf = chebapprox2(f, [ xdom ydom ], threshold);
    cc = chebcoeffs(cf);
catch
    cf = chebfun2(f, [ xdom ydom ]);
    cc = chebcoeffs2(cf);
end

U = zeros(m, 0);
V = zeros(n, 0);

x = linspace(xdom(1), xdom(2), n);
y = linspace(ydom(1), ydom(2), m);

fim1 = ones(m, 1);
ly = linspace(-1, 1, m).';
fi = ly;

U = [ fim1, fi ];

for i = 1 : size(cc, 1)
    if sum(abs(cc(i,:))) > threshold
        if i > 2
            fn = 2 * fi .* ly - fim1;
            fim1 = fi;
            fi = fn;
            U = [ U, fi ];
        end
        
        % Local version of the polynomial X (recall we need to change the
        % variable inside the Chebyshev polynomials)
        lx = linspace(-1,1,n).';
        
        % Needed for the Clenshaw recursive formula
        gim1 = ones(n,1);
        gi = lx;
        g = cc(i,1) * ones(n,1) + cc(i,2) * linspace(-1,1,n).';
        for j = 3 : size(cc, 2)
            gn = 2 .* lx .* gi - gim1;
            gim1 = gi;
            gi = gn;
            g = g + cc(i,j) * gn;
        end
        V = [ V, g ];
    end
end

end

