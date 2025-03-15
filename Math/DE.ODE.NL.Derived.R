
####################

### Helper Functions

source("Polynomials.Helper.R")


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

# TODO: check;

