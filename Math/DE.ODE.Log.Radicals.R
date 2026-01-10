########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Linear ODEs - Log w. Radicals
##
## draft v.0.1g

### Type
# Simple:  y = Log(Radical + P0(x))
# Coupled: y = Radical * Log(Radical + P0(x))

# Note:
# - for Radicals combined with simple Logs:
#   see file DE.ODE.Radicals.R;

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

### ODE:
4*(x + b0)*(x + b0 - d^2)*d2y + 2*(3*x + 3*b0 - d^2)*dy - 1 # = 0


# D =>
2*dy - 1/(sqrt(x + b0) + d) /sqrt(x + b0) # = 0
2*dy - (x + b0 - d*sqrt(x + b0))/(x + b0 - d^2) / (x + b0) # = 0
2*(x + b0)*(x + b0 - d^2)*dy - (x + b0 - d*sqrt(x + b0)) # = 0

# D2 =>
4*(x + b0)*(x + b0 - d^2)*d2y + 4*(2*x + 2*b0 - d^2)*dy +
	- (2 - d/sqrt(x + b0)) # = 0
4*(x + b0)^2*(x + b0 - d^2)*d2y + 4*(x + b0)*(2*x + 2*b0 - d^2)*dy +
	- (2*(x + b0) - d*sqrt(x + b0)) # = 0


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


#############

### y = (sqrt(x + b0) + a0) * log(sqrt(x + b0) + d)

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = -sqrt(2); a0 = 2/5; d = 1/3; c0 = -2/3;
params = list(x=x, b0=b0, a0=a0, d=d, c0=c0);
e = expression((sqrt(x + b0) + a0) * log(sqrt(x + b0) + d) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*(x+b0)*dy - sqrt(x+b0) * log(sqrt(x+b0) + d) - (sqrt(x+b0) + a0) * sqrt(x+b0) / (sqrt(x+b0) + d) # = 0
2*(x+b0)*(x+b0 - d^2)*(x+b0 - a0^2)*dy - (x+b0-d^2)*(x + b0 - a0*sqrt(x+b0)) * (y - c0) +
	- (x + b0 - a0^2) * ((x + b0 - a0*d)*sqrt(x+b0) + (a0-d)*(x+b0)) # = 0

# D2 =>
4*(x+b0)^2*(x+b0 - d^2)*(x+b0 - a0^2)*d2y +
	+ 2*(x+b0)*((x+b0)*(5*x+5*b0 - 4*a0^2 - 3*d^2) + 2*a0^2*d^2)*dy +
	+ a0*sqrt(x+b0) * 2*(x+b0)*(x+b0-d^2)*dy +
	- (2*(x+b0)*(2*x+2*b0-d^2) - a0*(3*x+3*b0-d^2)*sqrt(x+b0)) * (y - c0) +
	- (5*(x+b0)^2 - (3*a0^2+3*a0*d)*(x+b0) + a0^3*d) * sqrt(x+b0) +
	- 2*(2*x+2*b0 - a0^2) * (a0-d)*(x+b0) # = 0
4*(x+b0)^2*(x+b0 - d^2)*(x+b0 - a0^2)^2*d2y +
	+ 2*(x+b0)*(x+b0 - a0^2)*((x+b0)*(5*x+5*b0 - 4*a0^2 - 3*d^2) + 2*a0^2*d^2)*dy +
	- (x+b0 - a0*sqrt(x+b0)) * ((2*x+2*b0 - 2*a0^2)*(2*x+2*b0-d^2) + a0^2 * (x+b0-d^2)) * (y - c0) +
	- (x+b0 - a0^2) * (5*(x+b0)^2 - (3*a0^2+3*a0*d)*(x+b0) - a0*(a0-d)*(x+b0) + a0^3*d) * sqrt(x+b0) +
	- (x+b0 - a0^2) * (2*(a0-d)*(x+b0)*(2*x+2*b0 - a0^2) - a0*(x+b0)*(x + b0 - a0*d)) # = 0
4*(x+b0)^2*(x+b0 - a0^2)*(x+b0 - d^2)^2 * d2y +
	+ 2*(x+b0)*(x+b0 - a0^2)*(x+b0 - d^2)^2 * dy +
	- (x+b0)*(x+b0 - a0^2) * ((x+b0) - (3*d^2- 2*a0*d)) * sqrt(x+b0) +
	+ (x+b0)*(x+b0 - a0^2) * (a0*(x+b0) + a0*d^2 - 2*d^3) # = 0

# Analysis:

# SQRT:
((x+b0 - a0^2)*(x+b0 - a0*d) - a0*(x+b0-d^2)*(y - c0)) * sqrt(x+b0) # ==
2*(x+b0)*(x+b0 - d^2)*(x+b0 - a0^2)*dy - (x+b0-d^2)*(x + b0)*(y - c0) - (x+b0 - a0^2)*(a0-d)*(x+b0);

# Linear ODE: only when a0 = 0;

# Case: a0 = 0 =>
(x+b0) * sqrt(x+b0) # ==
2*(x+b0)*(x+b0 - d^2)*dy - (x+b0 - d^2)*(y - c0) + d*(x+b0);

# ODE:
4*(x+b0)^2*(x+b0 - d^2) * d2y + 4*d^2*(x+b0) * dy +
	+ (x+b0 - 3*d^2)*(y - c0) - d*(x+b0) # = 0


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


######################

### y = sqrt(x^2 + b0) * log(sqrt(x^2 + b0) + d)

# Check:
x = 4^(1/3); b0 = - sqrt(2); d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(sqrt(x^2 + b0) * log(sqrt(x^2 + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
dy - x*y / (x^2 + b0) - x/(sqrt(x^2 + b0) + d) # = 0
(x^2 + b0)*(x^2 + b0 - d^2)*dy - x*(x^2 + b0 - d^2)*y +
	- x*(x^2 + b0)*(sqrt(x^2 + b0) - d) # = 0

# D2 =>
(x^2 + b0)*(x^2 + b0 - d^2)*d2y +
	+ x*(3*x^2 + 3*b0 - d^2)*dy - (3*x^2 + b0 - d^2)*y +
	- (4*x^2 + b0)*sqrt(x^2 + b0) + d*(3*x^2 + b0) # = 0


### ODE:
x*(x^2 + b0)^2*(x^2 + b0 - d^2)*d2y +
	- (x^2 + b0)*(x^4 + (2*b0 - 3*d^2)*x^2 + b0^2 - d^2*b0) * dy +
	+ x^3*(x^2 + b0 - 3*d^2)*y - d*x^3*(x^2 + b0) # = 0


### Special Cases:

### b0 = d^2
x = sqrt(3); d = - 2^(1/4); b0 = d^2;
params = list(x=x, b0=b0, d=d);
# Note: value of d (+ or - sqrt(b0)) affects y;
# Re-run assignments above:
x*(x^2 + d^2)^2*d2y - (x^2 + d^2)*(x^2 - d^2)*dy +
	+ x*(x^2 - 2*d^2)*y - d*x*(x^2 + d^2) # = 0

##############

### y = (sqrt(x^2 + b0) + a0) * log(sqrt(x^2 + b0) - x)

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = sqrt(2); a0 = 2/5; c0 = -2/3;
params = list(x=x, b0=b0, a0=a0, c0=c0);
e = expression((sqrt(x^2 + b0) + a0) * log(sqrt(x^2 + b0) - x) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + b0)*dy - x*sqrt(x^2 + b0) * log(sqrt(x^2 + b0) - x) +
	+ (x^2 + b0 + a0*sqrt(x^2 + b0)) # = 0
# TODO


##############

### y = (sqrt(x^2 + b0) + x) * log(sqrt(x^2 + b0) - x)

# Check:
# for Quasi-Homogenous: c0 = 0;
x = sqrt(3); b0 = sqrt(2) + 1; c0 = -5/7;
params = list(x=x, b0=b0, c0=c0);
e = expression((sqrt(x^2 + b0) + x) * log(sqrt(x^2 + b0) - x) + c0)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### ODE:
# Note: has coefficient with SQRT;
(x^2 + b0)*d2y + x*dy - (y - c0) + 2*x + 2*sqrt(x^2 + b0) # = 0


# D =>
(x^2 + b0)*dy - sqrt(x^2 + b0) * (y - c0) +
	+ x^2 + b0 + x*sqrt(x^2 + b0) # = 0

# D2 =>
(x^2 + b0)*d2y + 2*x*dy - sqrt(x^2 + b0) * dy - x/sqrt(x^2 + b0) * (y - c0) +
	+ 2*x + sqrt(x^2 + b0) + x^2/sqrt(x^2 + b0) # = 0
(x^2 + b0)*d2y + 2*x*dy - (x/sqrt(x^2 + b0) + 1) * (y - c0) +
	+ 3*x + 2*sqrt(x^2 + b0) + x^2/sqrt(x^2 + b0) # = 0
(x^2 + b0)*d2y + x*dy - (y - c0) + 2*x + 2*sqrt(x^2 + b0) # = 0


####################

### Gen: y = x * log(sqrt(x^2 + b1*x + b0) + d)

# Check:
x = sqrt(5); b1 = - sqrt(3); b0 = sqrt(2); d = 1/3;
params = list(x=x, b0=b0, d=d);
e = expression(x * log(sqrt(x^2 + b1*x + b0) + d))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
2*x*dy - 2*y - x^2*(2*x + b1)/(sqrt(x^2 + b1*x + b0) + d) /sqrt(x^2 + b1*x + b0) # = 0
2*x*dy - 2*y - x^2*(2*x + b1)*
	(x^2 + b1*x + b0 - d*sqrt(x^2 + b1*x + b0))/(x^2 + b1*x + b0 - d^2) / (x^2 + b1*x + b0) # = 0
2*x*(x^2 + b1*x + b0)*(x^2 + b1*x + b0 - d^2)*dy +
	- 2*(x^2 + b1*x + b0)*(x^2 + b1*x + b0 - d^2)*y +
	- x^2*(2*x + b1)*(x^2 + b1*x + b0 - d*sqrt(x^2 + b1*x + b0)) # = 0

# D2 =>
4*x*(x^2 + b1*x + b0)*(x^2 + b1*x + b0 - d^2)*d2y +
	+ 4*x*(2*x + b1)*(2*x^2 + 2*b1*x + 2*b0 - d^2)*dy +
	- 4*(2*x + b1)*(2*x^2 + 2*b1*x + 2*b0 - d^2)*y +
	+ 4*x*(3*x + b1) * d*sqrt(x^2 + b1*x + b0) +
		+ x^2*(2*x + b1)^2 * d*sqrt(x^2 + b1*x + b0) / (x^2 + b1*x + b0) +
	- 2*(10*x^4 + 12*b1*x^3 + 6*b0*x^2 + 3*b1^2*x^2 + 2*b0*b1*x) # = 0

### ODE:
4*x^2*(2*x + b1)*(x^2 + b1*x + b0)*(x^2 + b1*x + b0 - d^2)*d2y +
	- 2*x*(4*b1*x^4 + (12*b0 + 5*b1^2 - 8*d^2)*x^3 +
		+ b1*(20*b0 + b1^2 - 12*d^2)*x^2 +
		+ (12*b0^2 - 12*d^2*b0 + 5*b0*b1^2 - 3*d^2*b1^2)*x +
		+ 4*b0*b1*(b0 - d^2)) * dy +
	+ 2*(4*b1*x^4 + (5*b1^2 + 12*b0 - 8*d^2)*x^3 + b1*(b1^2 + 20*b0 - 12*d^2)*x^2 +
		+ (12*b0*(b0 - d^2) + b1^2*(5*b0 - 3*d^2))*x +
		+ 4*b0*b1*(b0 - d^2)) * y +
	- x^3*(2*x + b1)^3 # = 0


### Special Cases:

### b0 = d^2
d = - 2^(1/4); b0 = d^2;
params = list(x=x, b0=b0, d=d);
# Re-run assignments above:
4*x^2*(x + b1)*(2*x + b1)*(x^2 + b1*x + b0)*d2y +
	- 2*x*(4*b1*x^3 + (5*b1^2 + 4*b0)*x^2 + b1*(b1^2 + 8*b0)*x + 2*b0*b1^2) * dy +
	+ 2 * (4*b1*x^3 + (5*b1^2 + 4*b0)*x^2 + b1*(b1^2 + 8*b0)*x + 2*b0*b1^2) * y +
	- x^2*(2*x + b1)^3 # = 0

