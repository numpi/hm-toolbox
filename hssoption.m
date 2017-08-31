function opt = hssoption(key, value)
%HSSOPTION Set or get an option for the hss toolbox. 
%
% Valid options are: 
%   'block-size': Integer representing the minimum block size. 
%   'threshold': Value used for off-diagonal truncation. 

global hss_block_size
global hss_threshold

if isempty(hss_block_size)
	hss_block_size = 32;
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
		otherwise
			error('Unsupported option specified');
	end

end

