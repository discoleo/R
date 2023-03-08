

### Systems: Lambert


####################

### Helper Functions

# used to solve the NLS
library(rootSolve)
# Helper functions
source("Polynomials.Helper.Solvers.Num.R")

# Note:
# IF solve.all FAILS:
# - verify that FUN combines the Re & Im parts of x!


###################

### x^x = k
k = 2
x = exp(lambertWp(log(k)))
#
x^x


### x^n * log(x/b) = 1
n = 3
b = log(2)
#
x = b * exp(lambertWp(n/b^n) / n)
x^n * log(x/b) # = 1


###
x = exp(lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


#####################

### exp(x^2) + b1*x + b0 = 0

solve.exp2 = function(b, x0, ...) {
	b0 = b[2]; b1 = b[1];
	FUN = function(x) {
		x = x[1] + x[2]*1i;
		err = test.exp2(x, b=b);
		c(Re(err), Im(err));
	}
	if(b0 == 0) {
		# exp(2*x^2) = b1^2*x^2;
		# W(-2): NO real solution!
		if(any(Im(x0) == 0)) warning("Unlikely to find solution!");
		sol = solve.all(FUN, x0=x0, ...);
		return(sol);
	}
	# TODO
}
test.exp2 = function(x, b) {
	exp(x^2) + b[1]*x + b[2];
}

###
b = c(2, 3)
# TODO


###
b = c(2, 0)
# x0 = - sqrt(W(-2/b1^2) / -2);
# only the negative root seems valid;
x0 = - sqrt((-0.79402363 + 0.77011175i) / -2);
x0 = rbind(-1 - 1i, -1 + 1i); # works as well
sol = solve.exp2(b, x0=x0)
test.exp2(sol, b=b)

