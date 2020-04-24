function CheckTestResult(value, op, bound, description)
%CHECKTESTRESULT Check if the value is in the expected range

switch op
    case '<'
        passed = value < bound;
    case '<='
        passed = value <= bound;
    case '>'
        passed = value > bound;
    case '>='
        passed = value >= bound;
    case '=='
        passed = value == bound;
    otherwise
        error('Unsupported operator');
end

if isstruct(description)
    start_was_called = true;
    time = toc(description.timer);
    text = description.description;
    timestring = sprintf('%6.2f s', time);
    description = sprintf('[ %s ] %s', ...
        timestring, description.description);
else
    start_was_called = false;
end

if start_was_called
    % Cleanup the line
    for j = 1 : length(text) + 11
        fprintf('\b');
    end
end

if passed
    fprintf('[   OK   ] %s\n', description);
else
    fprintf('[ FAILED ] %s\n', description);
    fprintf('           %e %s %e failed\n', value, op, bound);
end

assert(passed);


