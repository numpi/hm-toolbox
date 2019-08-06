function hodlr_RunAllTests
%hodlr_RUNALLTESTS Run all the unit tests for @hodlr

addpath ../

hodlr_TestCreation;
hodlr_TestOperations;
hodlr_TestLyapunov;

rmpath ../

end

