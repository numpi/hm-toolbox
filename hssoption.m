function opt = hssoption(key, value)
%HSSOPTION Set or get an option for the hss toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'compression': String indicating the strategy of compression (qr or svd)
%
% The special option HSSOPTION('clear') can be used to load all the
% default values, clearing all the previous HSSOPTION commands.

global hss_block_size
global hss_threshold
global hss_compression

if strcmp(key, 'clear')
    if exist('value', 'var')
        error('Specifying a value is unsupported for the special option "clear"');
    end
    
    clear hss_block_size;
    clear hss_threshold;
    clear hss_compression;
    
    return;
end

if isempty(hss_block_size)
    hss_block_size = 256;
end

if isempty(hss_compression)
    hss_compression = 'svd';
end

if isempty(hss_threshold)
    hss_threshold = 1e-12;
end

if ~exist('key', 'var')
    error('Please specify a key');
end

if ~exist('value', 'var')
    switch key
        case 'block-size'
            opt = hss_block_size;
        case 'threshold'
            opt = hss_threshold;
        case 'compression'
            opt = hss_compression;
        otherwise
            error('Unsupported option specified');
    end
else
    switch key
        case 'block-size'
            if value <= 2
                error('minimum block size must be at least 3');
            else
                hss_block_size = value;
            end
        case 'threshold'
            if value < 0
                error('threshold has to be positive');
            else
                hss_threshold = max(eps, value);
            end
        case 'compression'
            if strcmp(value,'qr') || strcmp(value, 'svd')
                hss_compression = value;
            else
                error('Unsupported type of compression');
            end
        otherwise
            error('Unsupported option specified');
    end
    
end

