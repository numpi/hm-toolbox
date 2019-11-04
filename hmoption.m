function opt = hmoption(varargin)
%HMOPTION Deprecated syntax for HODLROPTION

warning('HMOPTION is deprecated. Please use HODLROPTION');

if nargin > 1
    hodlroption(varargin{:});
else
    opt = hodlroption(varargin{:});
end

end

