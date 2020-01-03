function hmatrix_RunAllTests
%hmatrix_RUNALLTESTS Run all the unit tests for @hmatrix

addpath ../

hmatrix_TestCreation;
hmatrix_TestOperations;

rmpath ../

end

