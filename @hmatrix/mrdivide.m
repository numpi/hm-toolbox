function H = mrdivide(H1, H2)
%MRDIVIDE Compute H1 / H2

% We rely on mldivide for now
H = (H2' \ H1')';

end

