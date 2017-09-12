function hss_RunAllTests
%RUNALLTESTS Run all the unit tests.

addpath ../

hss_TestCreation;

hss_TestOperations;


rmpath ../

end
