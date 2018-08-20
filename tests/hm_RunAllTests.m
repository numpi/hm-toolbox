function hm_RunAllTests
%HM_RUNALLTESTS Run all the unit tests for @hm

addpath ../

hm_TestCreation;
hm_TestLyapunov;

rmpath ../

end

