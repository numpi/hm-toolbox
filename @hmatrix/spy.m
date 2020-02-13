function spy(H)
%SPY Inspect the rank structure of H.
%
% SPY(H) draw a plot displaying the low-rank structure in the off-diagonal
%     blocks of H.

m = size(H, 1);
n = size(H, 2);

xborders = [ 1 : n, n * ones(1, m), n:-1:1, zeros(1, m) ];
yborders = [ zeros(1,n), 1:m, m * ones(1, n), m:-1:1 ];

% plot(xborders, yborders, 'Color', diag_color);
set(gca, 'YDir', 'reverse');

axis([ 0 n 0 m ]);
hold on;

xd = zeros(4, 0);
yd = zeros(4, 0);
xo = zeros(4, 0);
yo = zeros(4, 0);
[xd,yd, xo, yo] = spy_draw_block(H, 0, 0, hmatrixrank(H), xd, yd, xo, yo);

diag_color    = [ 0.1,  0.3,  0.95 ];
offdiag_color = [ 0.95,  0.95,  1.0  ];

patch(xd, yd, diag_color);
patch(xo, yo, offdiag_color);

hold off;

end

function [xd, yd, xo, yo] = spy_draw_block(H, xoffset, yoffset, rk, xd, yd, xo, yo)


mm = H.sz(1);
nn = H.sz(2);

if is_leafnode(H)
    if ~H.admissible        
        xd = [ xd, [ xoffset ; xoffset ; xoffset + nn ; xoffset + nn ] ];
        yd = [ yd, [ yoffset ; yoffset + mm ; yoffset + mm ; yoffset ] ];
    else
        thisrk = size(H.U, 2);
        text(xoffset + nn / 2, yoffset + mm / 2, ...
            sprintf('%d', thisrk), ...
            'HorizontalAlignment', 'center', 'Interpreter', 'none');
        xo = [ xo, [ xoffset ; xoffset ; xoffset + nn ; xoffset + nn ] ];
        yo = [ yo, [ yoffset ; yoffset + mm ; yoffset + mm ; yoffset ] ];
    end
else
    mm = H.A11.sz(1);
    nn = H.A11.sz(2);
    [xd,yd,xo,yo] = spy_draw_block(H.A11, xoffset, yoffset, rk, xd, yd, xo, yo);
    [xd,yd,xo,yo] = spy_draw_block(H.A22, xoffset + nn, yoffset + mm, rk, xd, yd, xo, yo);
    [xd,yd,xo,yo] = spy_draw_block(H.A12, xoffset + nn, yoffset, rk, xd, yd, xo, yo);
    [xd,yd,xo,yo] = spy_draw_block(H.A21, xoffset, yoffset + mm, rk, xd, yd, xo, yo);
end
end
