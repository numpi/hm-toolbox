function html_file = CompileExamples
%COMPILEGUIDE Compile the toolbox guide from .m files.

% Create a parpool if available, so we don't skew timings later on
parfor i = 1 : 10
    % Nothing to do here
end

if exist('stylesheet.xsl', 'file')
    publish_cmd = @(file) publish(file, 'stylesheet', 'stylesheet.xsl');
else
    publish_cmd = @(file) publish(file);
end

doc_files = { ...
    'hodlr_linear_system.m', ...
    'hodlr_lyapunov.m', ...
    'hss_lyapunov.m'
    };

main_file = 'examples.m';

for i = 1 : length(doc_files)
    publish_cmd(doc_files{i});
end

html_file = publish_cmd(main_file);

if exist('numpi.css', 'file')
    copyfile('numpi.css', 'html/numpi.css');
end


end

