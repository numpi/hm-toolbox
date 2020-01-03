function RunAllTests
%RUNALLTESTS Run all the unit tests

hodlroption('clear');
hssoption('clear');

hodlr_RunAllTests;
hss_RunAllTests;
hmatrix_RunAllTests;

hodlroption('clear');
hssoption('clear');

fprintf(' == HMTOOLBOX TEST SUITE == \n');
fprintf('  - All the tests have been completed successfully.\n');
fprintf('  - The options have been reset to their default values.\n');

end

