
####################

### Helper Functions

source("Polynomials.Helper.R")

library(deSolve)


subst.part = function(x, param) {
	if(! inherits(param, "list")) param = as.list(param);
	if(inherits(x, "expression")) x = x[[1]];
	do.call(substitute, list(x, param));
}


######################

### Order 1 => Order 2

# dy = b*x*y^2 + fx

# D =>
d2y - 2*b*x*y*dy - b*y^2 - dfx # = 0
# Subst =>
d2y - 2*b*x*y*(b*x*y^2 + fx) - b*y^2 - dfx # = 0

### ODE:
d2y - 2*b^2*x^2*y^3 - b*y^2 - 2*b*x*fx*y - dfx # = 0

### Variants:
d2y - b*x*y*dy - b^2*x^2*y^3 - b*y^2- b*x*fx*y - dfx # = 0

# Note:
# - Derived ODE: is 1 order higher and power of y is also higher;
# - Additional solution may be possible;


##########
### Gen 1:

# dy = b2*x*y^2 + b1*x*y + b0*y + fx

# D =>
d2y - 2*b2*x*y*dy - b1*x*dy - b0*dy - b2*y^2 - b1*y - dfx # = 0
# Subst =>
d2y - (2*b2*x*y + b1*x + b0)*(b2*x*y^2 + b1*x*y + b0*y + fx) +
	- b2*y^2 - b1*y - dfx # = 0

### ODE:
d2y - 2*b2^2*x^2*y^3 - (3*b1*b2*x^2 + 3*b0*b2*x + b2)*y^2 +
	- (b1^2*x^2 + 2*b2*x*fx + 2*b1*b0*x + b1 + b0^2)*y +
	- b1*x*fx - b0*fx - dfx # = 0

# p1 = as.pm(...)
# p1 = sort.dpm(p1)


### Example:
b = c(1,2,3)
names(b) = paste0("b", 2:0)

### ODE:
d2y - 2*x^2 * y^3 - (6*x^2 + 9*x + 1) * y^2 +
	- (4*x^2 + 2*x*fx + 12*x + 11) * y +
	- 2*x*fx - 3*fx - dfx # = 0

# Initial ODE:
dy - x*y^2 - 2*x*y - 3*y - fx # = 0

### Solve:
ord1.f = function(t, y, parms, ...) {
	x  = t;
	fx = parms$FUN(x);
	dy = x*y^2 + 2*x*y + 3*y + fx;
	return(list(dy));
}
ord2.f = function(t, y, parms, ...) {
	x  = t; dy = y[2]; y = y[1];
	fx = parms$FUN(x); dfx = params$DFUN(x);
	d2y = 2*x^2 * y^3 + (6*x^2 + 9*x + 1) * y^2 +
		+ (4*x^2 + 2*x*fx + 12*x + 11) * y +
		+ 2*x*fx + 3*fx + dfx;
	return(list(c(dy, d2y)));
}

### Test:
params = list(
	FUN  = \(x) - x^2/3 - 3*x - 1,
	DFUN = \(x) - 2/3*x - 3);
x  = seq(0, 4, by = 1/4096 / 4);
y  = c(y = -1/2); # y = c(y = -2/7);
y2 = c(y = y, dy = 3*y + params$FUN(0))

sol1 = ode(y,  x, func = ord1.f, parms = params)
sol2 = ode(y2, x, func = ord2.f, parms = params, type="rk4")
col = "#FF6432A0"
matplot(sol1[,1], sol1[,-1], type = "l", lwd=2)
matplot(sol2[,1], sol2[, 2], type = "l", lwd=3, col=col, add=T)


# TODO:
# Check what happens for x > 2?
# Just numerical instabilities?



