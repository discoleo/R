########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Log w. Radicals
##
## draft v.0.1e



####################

### Helper Functions

library(deSolve)

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


#######################
#######################

### Simple: SQRT

### y = log(sqrt(x + b0) + d)

# Check:
x = sqrt(3); b0 = -sqrt(2); d = 1/3; params = list(x=x, b0=b0, d=d);
e = expression(log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*dy - 1/(sqrt(x + b0) + d) /sqrt(x + b0) # = 0
2*dy - (x + b0 - d*sqrt(x + b0))/(x + b0 - d^2) / (x + b0) # = 0
2*(x + b0)*(x + b0 - d^2)*dy - (x + b0 - d*sqrt(x + b0)) # = 0

# D2 =>
4*(x + b0)*(x + b0 - d^2)*d2y + 4*(2*x + 2*b0 - d^2)*dy +
	- (2 - d/sqrt(x + b0)) # = 0
4*(x + b0)^2*(x + b0 - d^2)*d2y + 4*(x + b0)*(2*x + 2*b0 - d^2)*dy +
	- (2*(x + b0) - d*sqrt(x + b0)) # = 0

### ODE:
4*(x + b0)*(x + b0 - d^2)*d2y + 2*(3*x + 3*b0 - d^2)*dy - 1 # = 0


#####################

### y = x * log(sqrt(x + b0) + d)

# Check:
x = sqrt(3); b0 = - sqrt(2); d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(x * log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 2*y - x^2/(sqrt(x + b0) + d) /sqrt(x + b0) # = 0
2*x*(x + b0)*(x + b0 - d^2)*dy +
	- 2*(x + b0)*(x + b0 - d^2)*y - x^2*(x + b0 - d*sqrt(x + b0)) # = 0

# D2 =>
4*x*(x + b0)^2*(x + b0 - d^2)*d2y +
	+ 4*(x + b0)*(2*x^2 + 2*b0*x - d^2*x) * dy +
	- 4*(x + b0)*(2*x + 2*b0 - d^2)*y +
	+ (4*(x + b0) + x) * d*x*sqrt(x + b0) +
	- (x + b0)*(6*x^2 + 4*b0*x) # = 0
4*x^2*(x + b0)^2*(x + b0 - d^2)*d2y +
	- 2*x*(x + b0)*(x^2 + (5*b0 - 3*d^2)*x + 4*b0^2 - 4*d^2*b0) * dy +
	- 4*x*(x + b0)*(2*x + 2*b0 - d^2)*y +
	+ 2*(x + b0)*(4*(x + b0) + x) * (x + b0 - d^2)*y +
	- x^3*(x + b0) # = 0

### ODE:
4*x^2*(x + b0)*(x + b0 - d^2)*d2y +
	- 2*x*(x^2 + (5*b0 - 3*d^2)*x + 4*b0^2 - 4*d^2*b0) * dy +
	+ 2*(x^2 + (5*b0 - 3*d^2)*x + 4*b0^2 - 4*d^2*b0) * y - x^3 # = 0


### Special Cases:

### b0 = d^2
x = sqrt(3); b0 = sqrt(2); d = - 2^(1/4);
params = list(x=x, b0=b0, d=d);
# Re-run assignments above:
4*x^2*(x + b0)*d2y - 2*x*(x + 2*b0)*dy + 2*(x + 2*b0)*y - x^2 # = 0


####################

### y = sqrt(x + b0) * log(sqrt(x + b0) + d)

# Check:
x = sqrt(3); b0 = - sqrt(2); d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(sqrt(x + b0) * log(sqrt(x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*dy - y / (x+b0) - 1/(sqrt(x + b0) + d) # = 0
2*(x+b0)*(x + b0 - d^2)*dy - (x + b0 - d^2)*y - (x+b0)*(sqrt(x + b0) - d) # = 0

# D2 =>
4*(x+b0)*(x + b0 - d^2)*d2y + 2*(3*x + 3*b0 - d^2)*dy +
	- 2*y - 3*sqrt(x + b0) + 2*d # = 0


### ODE:
4*(x+b0)^2*(x + b0 - d^2)*d2y + 4*d^2*(x+b0)*dy +
	+ (x + b0 - 3*d^2)*y - d*(x+b0) # = 0


### Special Cases:

### b0 = d = 1
b0 = 1; d = 1;
params = list(x=x, b0=b0, d=d);
# Re-run assignments above:
4*x*(x+1)^2*d2y + 4*(x+1)*dy + (x - 2)*y - (x+1) # = 0


####################
####################

### y = log(sqrt(x^2 + b0) + d)

# Check:
x = sqrt(3); b0 = -sqrt(2); d = 1/3; params = list(x=x, b0=b0, d=d);
e = expression(log(sqrt(x^2 + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
dy - x/(sqrt(x^2 + b0) + d) /sqrt(x^2 + b0) # = 0
dy - x*(x^2 + b0 - d*sqrt(x^2 + b0))/(x^2 + b0 - d^2) / (x^2 + b0) # = 0
(x^2 + b0)*(x^2 + b0 - d^2)*dy - x*(x^2 + b0 - d*sqrt(x^2 + b0)) # = 0

# D2 =>
(x^2 + b0)*(x^2 + b0 - d^2)*d2y +
	+ 2*x*(2*x^2 + 2*b0 - d^2)*dy +
	+ 2*d*sqrt(x^2 + b0) - b0*d*sqrt(x^2 + b0) / (x^2 + b0) +
	- (3*x^2 + b0) # = 0
(x^2 + b0)^2*(x^2 + b0 - d^2)*d2y +
	+ 2*x*(x^2 + b0)*(2*x^2 + 2*b0 - d^2)*dy +
	+ (2*x^2 + b0) * d*sqrt(x^2 + b0) +
	- (x^2 + b0)*(3*x^2 + b0) # = 0


### ODE:
x*(x^2 + b0)*(x^2 + b0 - d^2)*d2y +
	+ (2*x^4 + b0*x^2 - b0^2 + d^2*b0)*dy - x^3 # = 0


### Special Cases:

### b0 = d^2
x = sqrt(3); b0 = sqrt(2); d = - 2^(1/4);
params = list(x=x, b0=b0, d=d);
# Re-run assignments above:
x*(x^2 + b0)*d2y + (2*x^2 + b0)*dy - x # = 0


####################

### y = x * log(sqrt(x^2 + b0) + d)

# Check:
x = sqrt(3); b0 = -sqrt(2); d = 1/3; params = list(x=x, b0=b0, d=d);
e = expression(x * log(sqrt(x^2 + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*dy - y - x^3/(sqrt(x^2 + b0) + d) /sqrt(x^2 + b0) # = 0
x*dy - y - x^3*(x^2 + b0 - d*sqrt(x^2 + b0))/(x^2 + b0 - d^2) / (x^2 + b0) # = 0
x*(x^2 + b0)*(x^2 + b0 - d^2)*dy +
	- (x^2 + b0)*(x^2 + b0 - d^2)*y +
	- x^3*(x^2 + b0 - d*sqrt(x^2 + b0)) # = 0

# D2 =>
x*(x^2 + b0)*(x^2 + b0 - d^2)*d2y +
	+ 2*x^2*(2*x^2 + 2*b0 - d^2) * dy +
	- 2*x*(2*x^2 + 2*b0 - d^2)*y +
	+ (4 - b0/(x^2 + b0)) * x^2 * d*sqrt(x^2 + b0) +
	- x^2*(5*x^2 + 3*b0) # = 0


### ODE:
x^2*(x^2 + b0)*(x^2 + b0 - d^2)*d2y +
	- x*((3*b0 - 2*d^2)*x^2 + 3*b0^2 - 3*d^2*b0) * dy +
	+ ((3*b0 - 2*d^2)*x^2 + 3*b0^2 - 3*d^2*b0) * y - x^5 # = 0


### Special Cases:

### b0 = d^2
x = sqrt(3); b0 = sqrt(2); d = - sqrt(b0);
params = list(x=x, b0=b0, d=d);
# Re-run assignments above:
x^2*(x^2 + b0)*d2y - b0*x*dy + b0*y - x^3 # = 0

