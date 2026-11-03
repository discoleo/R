########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Atan( Log )
##
## draft v.0.1a


### NL ODEs: Atan-Type
#   Subtype: Atan( Log )


### Examples:

# x^2 * d2y^2 + 2*x * dy*d2y + 4*x^2 * dy^4 - 4*k*x*dy^3 + dy^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


########################

### Simple Components
# (Non-Coupled)

### y = atan(log(x)/k)

# Check:
e = expression(atan(log(x)/k));
k = sqrt(5);
x = sqrt(3);
params = list(x=x, k=k);
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


# D =>
x*(log(x)^2 + k^2) * dy - k # = 0

# D2 =>
k*x*d2y + k*dy + 2*x*log(x) * dy^2 # = 0
(k*x*d2y + k*dy)^2 - 4*x^2*log(x)^2 * dy^4 # = 0
k*(x*d2y + dy)^2 + 4*x * (k*x*dy - 1) * dy^3 # = 0

### ODE:
k*x^2 * d2y^2 + 2*k*x * dy*d2y + 4*k*x^2 * dy^4 - 4*x*dy^3 + k*dy^2 # = 0
# Re-formulated:
K = 1/k;
x^2 * d2y^2 + 2*x * dy*d2y + 4*x^2 * dy^4 - 4*K*x*dy^3 + dy^2 # = 0

