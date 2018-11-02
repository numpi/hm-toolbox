function [V, D] = eigs(varargin)
%EIGS Compute extremal eigenvalues of an HSS matrix.

p = inputParser;

p.addRequired('H');
p.addOptional('k', 5);
p.addOptional('opts', 'lm', @ischar);

p.parse(varargin{:});
args = p.Results;

n = size(args.H, 2);

switch args.opts
    case { 'lm', 'LM' }
        if nargout == 2
            [V, D] = eigs(@(x) args.H * x, n, args.k, args.opts);
        else
            V = eigs(@(x) args.H * x, n, args.k, args.opts);
        end
        
    case { 'sm', 'SM' }
        F = ulv(args.H);
        if nargout == 2
            [V, D] = eigs(@(x) ulv_solve(F, x), n, args.k, args.opts);
        else
            V = eigs(@(x) ulv_solve(F, x), n, args.k, args.opts);
        end
        
    otherwise
        error('Unsupported eigenvalue selection');
end

end

