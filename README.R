Survival4 for R version 0.60 or so
-----------------------------------

This is a port of nearly all the survival4 package from S, with the
Changelog date 17/9/96 (making it contemporary with S-PLUS 3.4). It is
not GNU Public Licensed -- the author of the original version, Terry
Therneau, retains all the rights but has given permission for the R
port to be made and distributed.  The code was not "ported" in the
usual sense -- many of the changes were to R rather than survival4.

The mortality tables are in a separate library "ratetables" since they're
HUGE. They will be autoloaded on demand.

To install the library, untar this file and type
  R INSTALL survival4
at the shell prompt.

You will also need the splines library (for plot.cox.zph()) and, if
you want to run all the examples, the date library.

There are a lot of examples in the Examples directory. You will
probably want to move this directory to somewhere more accessible.
There is a test file that runs nearly all the Examples.

It has been slightly modified to remove some things that don't work and to
allow for the fact that different calculations that are exactly equal in S
are only equal to 14 or 15 digits in R. I'm not sure exactly why this is.

Things that don't work yet
--------------------------
predict(,type="terms") doesn't work until someone writes labels.lm()
