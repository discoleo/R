

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
		# W(-2/b1^2): NO real solution!
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


###################
###################

### System:
# exp(x) = b*y + R
# exp(y) = b*z + R
# exp(z) = b*x + R

# the actual NLS:
solve.SExp = function(x, R, bb=1) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	y = exp(x) - bb*x[c(2,3,1)] - R;
	y = rbind(Re(y), Im(y));
	return(y);
}

# Parameters:
b = 1;
R = 2;
	
### Step 1:
	
# - choosing some non-standard values may help;
# - Note: exponentials may easily blow up;
x0 = c(2,2,1/3) + 1i*sqrt(2)*c(-2, 2, -1/4);
R0 = exp(x0) - b*x0[c(2,3,1)]
# create a seq from Rstart to Rend;
path = expand.path(R0, R)

### Step 2:
x = solve.path(solve.SExp, x0, path=path, bb=b)

### Test
exp(x) - b*x[c(2,3,1)]

# Non-Trivial Solution:
print(x)

