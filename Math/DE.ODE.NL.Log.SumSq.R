########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Sum of Log^2
##
## draft v.0.1c


####################
### Logarithmic  ###
### Higher Power ###
####################


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


######################
######################

### Type: 2 Components

### Derived from:
### y = (log(P1(x)))^2 + (log(P2(x)))^2
# ODE Type: d2y^2


### Example:
y = x^p * (log(x + a)^2 + log(x + b))^2
# - Simpler example;
# Note:
# - Coefficients of 2 Logs have to be identical,
#   except for a constant scaling;

# Check:
a = sqrt(2); b = sqrt(3); # b = -a;
p = 2/3;
x = 5^(2/3);
params = list(x=x, a=a, b=b, p=p);
e = expression(x^p * (log(x+a)^2 + log(x+b)^2))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x+a)*(x+b)*dy - p*(x+a)*(x+b)*y - 2*x^(p+1) * ((x+b) * log(x+a) + (x+a) * log(x+b)) # = 0

# D2 =>
x*(x+a)*(x+b) * d2y - ((p-3)*x^2 + (p-2)*(a+b)*x + (p-1)*a*b) * dy +
	- p*(2*x + a+b) * y +
	- 2*(p+1)*x^p * ((x+b) * log(x+a) + (x+a) * log(x+b)) +
	- 2*x^(p+1) * (log(x+a) + log(x+b)) +
	- 2*x^(p+1) * ((x+b)/(x+a) + (x+a)/(x+b)) # = 0
x^2*(x+a)*(x+b) * d2y - x * ((2*p-2)*x^2 + (2*p-1)*(a+b)*x + 2*p*a*b) * dy +
	+ p*((p-1)*x^2 + p*(a+b)*x + (p+1)*a*b) * y +
	- 2*x^(p+2) * (log(x+a) + log(x+b)) +
	- 2*x^(p+2) * ((x+b)/(x+a) + (x+a)/(x+b)) # = 0

# Linear System:
F = x^2*(x+a)*(x+b) * d2y - x * ((2*p-2)*x^2 + (2*p-1)*(a+b)*x + 2*p*a*b) * dy +
	+ p*((p-1)*x^2 + p*(a+b)*x + (p+1)*a*b) * y +
	- 2*x^(p+2) * ((x+b)/(x+a) + (x+a)/(x+b));
log(x+a) # ==
(x*(x*(x+a)*(x+b)*dy - p*(x+a)*(x+b)*y) - (x+a)*F) / (2*(b-a)*x^(p+2));
(x+a)*(x^2*(x+b)*dy - p*x*(x+b)*y - F) / (2*(b-a)*x^(p+2));
#
log(x+b) # ==
(x*(x*(x+a)*(x+b)*dy - p*(x+a)*(x+b)*y) - (x+b)*F) / -(2*(b-a)*x^(p+2));
(x+b)*(x^2*(x+a)*dy - p*x*(x+a)*y - F) / -(2*(b-a)*x^(p+2));

# Substitution in Eq for y:
4*(b-a)^2*x^(p+4) * y +
	- (x+a)^2*(F - x^2*(x+b)*dy + p*x*(x+b)*y)^2 +
	- (x+b)^2*(F - x^2*(x+a)*dy + p*x*(x+a)*y)^2 # = 0

# TODO: simplify;


######################

### Derived from:
### y = (log(P1(x)))^2 + (log(P2(x)))^2

### Example:
y = (log(x^2 + a))^2 + (log(x^2 + b))^2

# Check:
a = sqrt(2); b = sqrt(3); # b = -a;
x = 5^(2/3);
params = list(x=x, a=a, b=b);
e = expression((log(x^2 + a))^2 + (log(x^2 + b))^2)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D(y):
dy - 4*x*log(x^2 + a) / (x^2 + a) - 4*x*log(x^2 + b) / (x^2 + b) # = 0
#
(x^2+a)*(x^2+b)*dy +
	- 4*x*(x^2+b)*log(x^2 + a) - 4*x*(x^2+a)*log(x^2 + b) # = 0

### D2(y):
(x^2+a)*(x^2+b) * d2y + (4*x^3 + 2*(a+b)*x) * dy +
	- 4*(3*x^2 + b)*log(x^2 + a) - 4*(3*x^2 + a)*log(x^2 + b) +
	- 8*x^2*(x^2+b) / (x^2 + a) - 8*x^2*(x^2+a) / (x^2 + b) # = 0

### Solve Linear system:
# T = 8*x^2*(x^2+b) / (x^2 + a) + 8*x^2*(x^2+a) / (x^2 + b) - (4*x^3 + 2*(a+b)*x)*dy;
T = 8*x^2*(2*x^4 + 2*(a+b)*x^2 + a^2 + b^2) / ((x^2 + a)*(x^2 + b)) - (4*x^3 + 2*(a+b)*x)*dy;
Z = 8*x^2*(2*x^4 + 2*(a+b)*x^2 + a^2 + b^2) / ((x^2 + a)*(x^2 + b));
M = (x^2+a)*(x^2+b);
#
log(x^2 + a) # ==
	(x*(x^2+a)*M*d2y - (3*x^2 + a)*M*dy - x*(x^2+a)*T) / (8*(a-b)*x^3);
	(x^2+a)*(x*M*d2y + (x^4 + (a-b)*x^2 - a*b)*dy - x*Z) / (8*(a-b)*x^3);
log(x^2 + b) # ==
	(x*(x^2+b)*M*d2y - (3*x^2 + b)*M*dy - x*(x^2+b)*T) / -(8*(a-b)*x^3);
	(x^2+b)*(x*M*d2y + (x^4 - (a-b)*x^2 - a*b)*dy - x*Z) / -(8*(a-b)*x^3);


### ODE:
64*(a-b)^2 * x^6 * y +
	- (x^2+a)^2 * (x*M*d2y + (x^4 + (a-b)*x^2 - a*b)*dy - x*Z)^2 +
	- (x^2+b)^2 * (x*M*d2y + (x^4 - (a-b)*x^2 - a*b)*dy - x*Z)^2 # = 0

# TODO: simplify;

### Special Cases:

### Case: b = -a;
Z = 16*x^2*(x^4 + a^2) / (x^4 - a^2);
M = x^4 - a^2;
256*a^2 * x^6 * y +
	- (x^2+a)^2 * (x*M*d2y + (x^2 + a)^2*dy - x*Z)^2 +
	- (x^2-a)^2 * (x*M*d2y + (x^2 - a)^2*dy - x*Z)^2 # = 0
256*a^2 * x^6 * y +
	- (x^2+a)^2 * (x*M*d2y + (x^2 + a)^2*dy)^2 +
	- (x^2-a)^2 * (x*M*d2y + (x^2 - a)^2*dy)^2 +
	+ 64*x^4*(x^4 + a^2)^2 * d2y +
	+ 64*x^3*(x^4 + a^2) * (x^8 + 6*a^2*x^4 + a^4) / (x^4-a^2) * dy +
	- 256*x^6*(x^4 + a^2)^2 * (1/(x^2-a)^2 + 1/(x^2+a)^2) # = 0

