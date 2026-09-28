#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log * Exp
##
## draft v.0.1d

### Log to Power:
# y = Log(P1(x))^2 * Exp(P2(x)})

### Motto:
# The new AI proofs in Maths look only more and more ridiculous!


### Examples:

# 2*x*y*d2y - x*dy^2 - 2*(k*x - 1) * y*dy + k*(k*x - 2) * y^2 = 0;
# 3*x*y*d2y - 2*x*dy^2 - (2*k*x - 3) * y*dy + k*(k*x - 3) * y^2 = 0;
# n*x*y*d2y - (n-1)*x*dy^2 - (2*k*x - n) * y*dy + k*(k*x - n) * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Log(x)^2 * Exp(k*x)

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(log(x)^2 * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 2*log(x)*exp(k*x) # = 0

# D2 =>
x*d2y + dy - k*x*dy - k*y +
	- 2*k*log(x)*exp(k*x) - 2*exp(k*x)/x # = 0
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y - 2*exp(k*x) # = 0
2*x^2*d2y - 2*x*(2*k*x - 1)*dy + 2*k*x*(k*x - 1)*y - (x*dy - k*x*y)^2 / y # = 0
2*x^2 * y*d2y - 2*x*(2*k*x - 1) * y*dy + 2*k*x*(k*x - 1) * y^2 +
	- (x*dy - k*x*y)^2 # = 0

### ODE:
2*x*y*d2y - x*dy^2 - 2*(k*x - 1) * y*dy + k*(k*x - 2) * y^2 # = 0


#########################

### Type: 2 Components w. Exp
# y = Log(x)^2 * Exp(k*x) + c1*Exp(k*x)

# Check:
k  = 1/sqrt(5);
c1 = 2^(1/3);
x = sqrt(3); params = list(x=x, k=k, c1=c1);
e = expression(log(x)^2 * exp(k*x) + c1*exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 2*log(x)*exp(k*x) # = 0

# D2 =>
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y - 2*exp(k*x) # = 0

# =>
2*log(x)^2 * exp(k*x) # ==
- c1*x^2*d2y + c1*x*(2*k*x - 1)*dy - c1*k*x*(k*x - 1)*y + 2*y;
# =>
x^2*d2y - x*(2*k*x - 1)*dy + k*x*(k*x - 1)*y +
	+ (x*dy - k*x*y)^2 / (c1*x^2*d2y - c1*x*(2*k*x - 1)*dy + c1*k*x*(k*x - 1)*y - 2*y) # = 0

# TODO: simplify;


#####################
#####################

#####################
### Higher Powers ###

### y = Log(x)^3 * Exp(k*x)

# Check:
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, k=k);
e = expression(log(x)^3 * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - 3*log(x)^2 * exp(k*x) # = 0

# D2 =>
x^2 * d2y - x*(k*x - 1)*dy - k*x*y +
	- 3*k*x*log(x)^2 * exp(k*x) - 6*log(x) * exp(k*x) # = 0
x^2 * d2y - x*(k*x - 1)*dy - k*x*y +
	- k*x * (x*dy - k*x*y) - 6*log(x) * exp(k*x) # = 0
3*x^2 * y*d2y - 3*x*(2*k*x - 1) * y*dy + 3*k*x*(k*x - 1) * y^2 +
	- 2 * (x*dy - k*x*y)^2 # = 0

### ODE:
3*x*y*d2y - 2*x*dy^2 - (2*k*x - 3) * y*dy + k*(k*x - 3) * y^2 # = 0


################################

### Gen: y = Log(x)^n * Exp(k*x)

# Check:
n = sqrt(pi); # n = exp(1/pi);
k = 1/sqrt(5);
x = sqrt(3); params = list(x=x, n=n, k=k);
e = expression(log(x)^n * exp(k*x))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*dy - k*x*y - n*log(x)^(n-1) * exp(k*x) # = 0

### ODE:
n*x*y*d2y - (n-1)*x*dy^2 - (2*k*x - n) * y*dy + k*(k*x - n) * y^2 # = 0

