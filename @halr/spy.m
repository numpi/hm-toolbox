function spy(H, levels)
%SPY Inspect the rank structure of H.
%
% SPY(H) draw a plot displaying the low-rank structure in the off-diagonal
%     blocks of H.
%
% SPY(H, LEVELS) shows the ranks of the off-diagonal blocks at depth
% LEVELS. If not specified, LEVELS = 6. 

m = size(H, 1);
n = size(H, 2);

xborders = [ 1 : n, n * ones(1, m), n:-1:1, zeros(1, m) ];
yborders = [ zeros(1,n), 1:m, m * ones(1, n), m:-1:1 ];

% plot(xborders, yborders, 'Color', diag_color);
set(gca, 'YDir', 'reverse');

axis([ 0 n 0 m ]);
hold on;

% Construct a struct to hold the global state of recursion in
% spy_draw_block, constructing the data for the patch. 
spy_state = struct();

% Coordinates of the diagonal blocks
spy_state.xd = zeros(4, 0);
spy_state.yd = zeros(4, 0);

% Coordinate of the off-diagonal blocks with text
spy_state.xo = zeros(4, 0);
spy_state.yo = zeros(4, 0);

% Coordinates of the text, and the rank to write
spy_state.xt = [];
spy_state.yt = [];
spy_state.rkt = [];

% Coordinates of the off-diagonal blocks without text
spy_state.xon = zeros(4, 0);
spy_state.yon = zeros(4, 0);

if ~exist('levels', 'var')
    levels = 6;
end

spy_state = spy_draw_block(H, 0, 0, spy_state, levels);

diag_color    = [ 0.1,  0.3,  0.95 ];
offdiag_color = [ 0.95,  0.95,  1.0  ];
offdiag_no_color = [ 0.9,  0.9,  1.0  ];

patch(spy_state.xd, spy_state.yd, diag_color);
patch(spy_state.xo, spy_state.yo, offdiag_color);
patch(spy_state.xon, spy_state.yon, offdiag_no_color);

for j = 1 : length(spy_state.xt)
    text(spy_state.xt(j), spy_state.yt(j), ...
        sprintf('%d', spy_state.rkt(j)), ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none');
end

hold off;

end

function spy_state = spy_draw_block(H, xoffset, yoffset, spy_state, levels)


mm = H.sz(1);
nn = H.sz(2);

if is_leafnode(H)
    if ~H.admissible        
        spy_state.xd = [ spy_state.xd, ...
            [ xoffset ; xoffset ; xoffset + nn ; xoffset + nn ] ];
        spy_state.yd = [ spy_state.yd, ...
            [ yoffset ; yoffset + mm ; yoffset + mm ; yoffset ] ];
    else
        if levels >= 0
            spy_state.xt = [ spy_state.xt, xoffset + nn / 2 ];
            spy_state.yt = [ spy_state.yt, yoffset + mm / 2 ];
            spy_state.rkt = [ spy_state.rkt, size(H.U, 2) ];
            
            spy_state.xo = [ spy_state.xo, ...
                [ xoffset ; xoffset ; xoffset + nn ; xoffset + nn ] ];
            spy_state.yo = [ spy_state.yo, ...
                [ yoffset ; yoffset + mm ; yoffset + mm ; yoffset ] ];
        else
            spy_state.xon = [ spy_state.xon, ...
                [ xoffset ; xoffset ; xoffset + nn ; xoffset + nn ] ];
            spy_state.yon = [ spy_state.yon, ...
                [ yoffset ; yoffset + mm ; yoffset + mm ; yoffset ] ];
        end
      end
else
    mm = H.A11.sz(1);
    nn = H.A11.sz(2);
    spy_state = spy_draw_block(H.A11, xoffset, yoffset, spy_state, levels - 1);
    spy_state = spy_draw_block(H.A22, xoffset + nn, yoffset + mm, spy_state, levels - 1);
    spy_state = spy_draw_block(H.A12, xoffset + nn, yoffset, spy_state, levels - 1);
    spy_state = spy_draw_block(H.A21, xoffset, yoffset + mm, spy_state, levels - 1);
end
end
