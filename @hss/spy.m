function spy(H, variant)
%SPY Inspect the rank structure of H.
%
% SPY(H) draw a plot displaying the low-rank structure in the off-diagonal
%     blocks of H.
%
% SPY(H, 'svd') also plots bars representing the singular values of the
%     low-rank blocks, to indicate the decay present in those blocks. 

if ~exist('variant', 'var')
    variant = 'standard';
end

if strcmp('variant', 'svd')
    H = hss_proper(H);
end

m = size(H, 1);
n = size(H, 2);

xborders = [ 1 : n, n * ones(1, m), n:-1:1, zeros(1, m) ];
yborders = [ zeros(1,n), 1:m, m * ones(1, n), m:-1:1 ];

% plot(xborders, yborders, 'Color', diag_color);
set(gca, 'YDir', 'reverse');

axis([ 0 n 0 m ]);
hold on;

spy_draw_block(H, 0, 0, hssrank(H), variant);

hold off;

end

function spy_draw_block(H, xoffset, yoffset, rk, variant)
diag_color    = [ 0.1,  0.3,  0.95 ];
offdiag_color = [ 0.95,  0.95,  1.0  ];
svd_color = [1.0, 0.3, 0.4];

if H.leafnode == 1
    mm = size(H, 1);
    nn = size(H, 2);
    
    xb = [xoffset, xoffset + nn, xoffset + nn, xoffset ];
    yb = [yoffset, yoffset, yoffset + mm, yoffset + mm ];
    
    fill(xb, yb, diag_color);
else
    spy_draw_block(H.A11, xoffset, yoffset, rk, variant);
    mm = size(H.A11, 1);
    nn = size(H.A11, 2);
    spy_draw_block(H.A22, xoffset + nn, yoffset + mm, rk, variant);
    
    mm2 = size(H.A22, 1);
    nn2 = size(H.A22, 2);
    
    % Low rank block (1, 2)
    xb = [xoffset + nn, xoffset + nn + nn2, ...
        xoffset + nn + nn2, xoffset + nn ];
    yb = [yoffset, yoffset, yoffset + mm, yoffset + mm ];
    thisrk = rank(H.B12);
    
    fill(xb, yb, offdiag_color);

    if strcmp(variant, 'svd')
        % Draw the SVD plot of the block in the background
        s = log(svd(H.B12));
        nsvd = length(s);
        for j = 1 : nsvd      
            sjm = (s(1) - s(j)) / (s(1) - s(end)) * mm;
            xbj = [xoffset + nn + nn2 * (j-1)/nsvd, xoffset + nn + nn2 * j / nsvd, ...
                xoffset + nn + nn2 * j / nsvd, xoffset + nn + nn2 * (j-1)/nsvd ];
            ybj = [yoffset + sjm, yoffset + sjm, yoffset + mm, yoffset + mm ];
            fill(xbj, ybj, svd_color, 'LineStyle','none');
        end
    end

    text(xoffset + nn + nn2 / 2, yoffset + mm / 2, ...
        sprintf('%d', thisrk), ...
        'HorizontalAlignment', 'center');
    
    % Low rank block (2, 1)
    xb = [xoffset, xoffset + nn, xoffset + nn, xoffset ];
    yb = [yoffset + mm, yoffset + mm, ...
        yoffset + mm + mm2, yoffset + mm + mm2 ];
    thisrk = rank(H.B21);
    
    fill(xb, yb, offdiag_color);

    if strcmp(variant, 'svd')
        % Draw the SVD plot of the block in the background
        s = log(svd(H.B21));
        nsvd = length(s);
        for j = 1 : nsvd      
            sjm = (s(1) - s(j)) / (s(1) - s(end)) * mm2;
            xbj = [xoffset + nn * (j-1)/nsvd, xoffset + nn * j / nsvd, ...
                xoffset + nn * j / nsvd, xoffset + nn * (j-1)/nsvd ];
            ybj = [yoffset + mm + sjm, yoffset + mm + sjm, yoffset + mm + mm2, yoffset + mm + mm2 ];
            fill(xbj, ybj, svd_color, 'LineStyle','none');
        end
    end

    text(xoffset + nn / 2, yoffset + mm + mm2 / 2, ...
        sprintf('%d', thisrk), ...
        'HorizontalAlignment', 'center');
end
end
