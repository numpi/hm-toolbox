function H = mrdivide(H1, H2)

% We rely on mldivide, that implements all the corner cases.
H = ( H2' \ H1' )';

end