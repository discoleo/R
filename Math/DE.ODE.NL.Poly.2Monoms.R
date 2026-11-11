#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Polynomial Types
## Encoded: 2 Monomials
##
## draft v.0.1a


### Polynomial Type
# w. 2 encoded Monomials;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")

#########################


### y^5 - 5*x*y^3 + 5*x^2*y - 2*b*x = 0

# Check:
x = 1/sqrt(3) + 0i;
b = sqrt(5);
params = list(x=x, b=b);
ye = expression((b*x + sqrt(b^2*x^2 - x^5))^(1/5) + (b*x - sqrt(b^2*x^2 - x^5))^(1/5));
y  = eval(ye, params);
dye = D(ye, "x"); dy = eval(dye, params);
d2y = eval(D(dye, "x"), params);

# D =>
5*y^4*dy - 15*x*y^2*dy - 5*y^3 + 5*x^2*dy + 10*x*y - 2*b # = 0

# System: in y^5 & y^3
y^5 - 5*x*y^3 + 5*x^2*y - 2*b*x # = 0
5*y^5*dy - 5*(3*x*dy + y)*y^3 + 5*x^2*y*dy + 10*x*y^2 - 2*b*y # = 0
# =>
y^5 # ==
(- (3*x*dy + y)*(5*x^2*y - 2*b*x) + x*(5*x^2*y*dy + 10*x*y^2 - 2*b*y)) / (-2*x*dy + y);
(10*x^3*y*dy - 6*b*x^2*dy - 5*x^2*y^2) / (2*x*dy - y);
#
y^3 # ==
(- 5*dy*(5*x^2*y - 2*b*x) + (5*x^2*y*dy + 10*x*y^2 - 2*b*y)) / (-25*x*dy + 5*(3*x*dy + y));
(20*x^2*y*dy - 10*b*x*dy - 10*x*y^2 + 2*b*y) / (10*x*dy - 5*y);

# D2 =>
y^4*d2y - 3*x*y^2*d2y + x^2*d2y + 4*y^3*dy^2 - 6*x*y*dy^2 +
	- 6*y^2*dy + 4*x*dy + 2*y # = 0
y^5*d2y - 3*x*y^3*d2y + x^2*y*d2y + 4*y^4*dy^2 - 6*x*y^2*dy^2 +
	- 6*y^3*dy + 4*x*y*dy + 2*y^2 # = 0
5*(10*x^3*y*dy - 6*b*x^2*dy - 5*x^2*y^2) * d2y +
	- (3*x*d2y - 4*y*dy^2 + 6*dy) * (20*x^2*y*dy - 10*b*x*dy - 10*x*y^2 + 2*b*y) +
	+ (5*x^2*y*d2y - 30*x*y^2*dy^2 + 20*x*y*dy + 10*y^2) * (2*x*dy - y) # = 0

### ODE:
3*b*x*y*d2y - 10*(x^2*y^2 - 2*b*x*y) * dy^3 + (5*x*y^3 - 4*b*y^2 + 40*x^2*y - 30*b*x) * dy^2 +
	- 30*x*y^2*dy + 6*b*y*dy + 5*y^3 # = 0

