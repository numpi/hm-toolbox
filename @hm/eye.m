function I = eye(varargin)
%EYE Create an HM identity

sz = varargin{1};
I = hm('diagonal', ones(sz(1),1));

end

