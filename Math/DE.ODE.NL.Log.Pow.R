########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Powers of Log
##
## draft v.0.1b


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

### Example:
y = (log(x^2 + a))^2 + (log(x^2 + b))^2

# Check:
a = sqrt(2); b = sqrt(3); # b = -a;
x = 5^(4/3);
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

