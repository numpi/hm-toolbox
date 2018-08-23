function CheckTestResult(value, op, bound, description)
%CHECKTESTRESULT Check if the value is in the expected range

passed = false;

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

if passed
    fprintf('[   OK   ] %s\n', description);
else
    fprintf('[ FAILED ] %s\n', description);
    fprintf('           %e %s %e failed\n', value, op, bound);
end

assert(passed);


