This is the entire Examples directory from the survival4
distribution. It contains tests for nearly all aspects of the package.

-----------------------
Running the tests in R
-----------------------
All the tests that work in R have been placed in Rtest.

Type

% R -v6 -nosave < Rtest >Rtest.out

to run them. You can compare the output to what I got using

% diff Rtest.out Rtest.out.supplied

Note that this requires write access to the directory.


