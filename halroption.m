function opt = halroption(key, value)
%HALROPTION Set or get an option for the halr toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'compression': String indicating the strategy of compression (qr or svd)
%   'norm': Norm to use in truncations (2 or 'fro')

global halr_block_size
global halr_threshold
global halr_compression
global halr_norm

% Private options -- not exposed to users

% This option sets the maximum allowed rank for low-rank blocks with
% respect the dimension: for a block of size (m,n) the maximum allowed rank
% is halroption('rank-ratio') * min(m,n)
global halr_rank_ratio

if isempty(halr_rank_ratio)
    halr_rank_ratio = 0.5;
end

if isempty(halr_block_size)
    halr_block_size = 256;
end

if isempty(halr_compression)
    halr_compression = 'svd';
end

if isempty(halr_threshold)
    halr_threshold = 1e-12;
end

if isempty(halr_norm)
    halr_norm = 2;
end

if ~exist('key', 'var')
    error('Please specify a key');
end

if ~exist('value', 'var')
    switch key
        case 'block-size'
            opt = halr_block_size;
        case 'threshold'
            opt = halr_threshold;
        case 'compression'
            opt = halr_compression;
        case 'norm'
            opt = halr_norm;
        case 'rank-ratio'
            opt = halr_rank_ratio;
        otherwise
            error('Unsupported option specified');
    end
else
    switch key
        case 'block-size'
            if value <= 2
                error('minimum block size must be at least 3');
            else
                halr_block_size = value;
            end
        case 'threshold'
            if value < 0
                error('threshold has to be positive');
            else
                halr_threshold = max(eps, value);
            end
        case 'norm'
            if (isscalar(value) && value ~= 2) && ~strcmp(value, 'fro')
                error('Unsupported norm specified');
            else
                halr_norm = value;
            end
        case 'compression'
            if strcmp(value,'qr') || strcmp(value, 'svd')
                halr_compression = value;
            else
                error('Unsupported type of compression');
            end
        otherwise
            error('Unsupported option specified');
    end
    
end

