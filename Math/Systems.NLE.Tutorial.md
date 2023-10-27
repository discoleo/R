# Solving Non-Linear equations in R
**A brief Tutorial**

Solving systems of non-linear equations is often impeded by finding suitable starting values. The technique described here circumvents this problem. It is probably the simplest way to solve such systems.

> Note:
> Fitting non-linear models usually require a different set of techniques. This tutorial does NOT focus on such techniques.

## R Packages
This tutorial uses package *rootSolve*, but I presume that *nleqslv* works equally well.

## Problem

**SYSTEM:**\
	exp(x) = b * y + R;
	exp(y) = b * x + R;
	# where b, R = given parameters;

This system has the trivial solution x = y.

We want to find some non-trivial solutions with x != y, and possibly having complex values. Unfortunately, many of the starting values converge on the x = y roots. So we need a better approach.


## Approach

1. Step 1:
We will start with a different system:
- we **select** a set of values as the "solution" and basically "transform" the system to fit this "solution": i.e. we compute the 2 new "R" values for this new system;

2. Step 2:
We will "flow" the initial "solution", while changing the system back to the original system. Hopefully, the pseudo-solution will converge to a true/desired solution.


## R Code:

	# given parameters
	b = 1;
	R = 2;
	
	# used to solve the NLS
	library(rootSolve)
	# Helper functions
	source("Polynomials.Helper.Solvers.Num.R")
	
	### Step 1:
	
	# choosing some non-standard values may help;
	# (see also some tips later)
	x0 = 1i*sqrt(2)*c(-1, -6);
	R0 = exp(x0) - b*x0[c(2,1)]
	# create a seq from Xstart to Xend;
	path = expand.path(R0, R)
	
	### Step 2
	# - we need to "enhance" rootSolve
	#   to work also with complex roots:
	#   see file "Polynomials.Helper.Solvers.Num.R";

	# the actual NLS:
	solve.SExp = function(x, R, bb=1) {
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
		y = exp(x) - bb*x[c(2,1)] - R;
		y = rbind(Re(y), Im(y));
		return(y);
	}
	
	# and now lets solve the NLS:
	x = solve.path(solve.SExp, x0, path=path, bb=b)
	
	### Test
	exp(x) - b*x[c(2,1)]


### Alternative solution:
	x0 = c(2 + 1i * sqrt(7), 0.1 - 1i * sqrt(23));
	R0 = exp(x0) - b*x0[c(2,1)]
	path = expand.path(R0, R)
	x = solve.path(solve.SExp, x0, path=path, bb=b)


## Tips
- if the found solution is undesirable, then it is possible to choose different values for the "initial" solution;
- choosing c(0, ..., 0) is a very bad choice most of the time;
- choosing some irrational values (and including also negative values if negative roots are feasible) is possibly the best approach;
- if complex roots are not desirable: try a different set of values for the "initial" solution; "disabling" complex roots will probably fail, as the initial roots are likely to truly flow/converge to the complex values;

In case of multiple failed attempts to solve the system:
- various other methods are available, but they are much more complicated and messier;

I hope this helps solving a few more non-linear systems.


## GitHub

The code is also available on GitHub:
1. **Helper Functions:**
> https://github.com/discoleo/R/blob/master/Math/Polynomials.Helper.Solvers.Num.R
2. **Exponential example:**
> - see Section "Supplementary Info";
> https://github.com/discoleo/R/blob/master/Math/DE.Systems.Lambert.R


## References

### Varia

1. Fields Institute: Riemannian Optimization on Embedded Manifolds Using Homotopy Continuation
> https://www.youtube.com/watch?v=YT-2LvTeud4&list=PLArBKNfJxuun0G_hi4HMq4dH5a9v2jFdB&index=5
> see from around 14:00;

