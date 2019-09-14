function html_file = CompileExamples(varargin)
%COMPILEGUIDE Compile the toolbox guide from .m files.

if exist('stylesheet.xsl', 'file')
    publish_cmd = @(file) publish(file, 'stylesheet', 'stylesheet.xsl');
else
    publish_cmd = @(file) publish(file);
end

if nargin == 0
	doc_files = { ...
		'hodlr_linear_system.m', ...
		'hodlr_lyapunov.m', ...
		'hss_lyapunov.m', ...
		'hm_expm.m'
		};
else
	doc_files = varargin;
end

main_file = 'examples.m';

for i = 1 : length(doc_files)
	fprintf('Publishing file %s ... ', doc_files{i});
    publish_cmd(doc_files{i});
	fprintf('done\n');
end

html_file = publish_cmd(main_file);

if exist('numpi.css', 'file')
    copyfile('numpi.css', 'html/numpi.css');
end


end

