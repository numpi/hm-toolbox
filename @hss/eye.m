function I = eye(varargin)
%EYE Create an HM identity

sz = varargin{1};
I = hss('banded', speye(sz), 0, 0);

end
