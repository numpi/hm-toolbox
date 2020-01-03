function H = hm(varargin)
%HM Constructor for HODLR matrices
%
% H = HM(...) is equivalent to H = HODLR(...), and uses an old deprecated
% syntax. Type help hodlr for further details.

warning('The HM constructor is deprecated. Please use HODLR instead');

H = hodlr(varargin{:});

end

