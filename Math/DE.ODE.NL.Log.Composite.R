########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Composite Logarithms
##
## draft v.0.2b

###########################
### Compositions of Log ###
###########################

### y = log(exp(P1(x)) + P2(x)) + F0(x)

### y = log(exp(x^n) + b0) + f0

# Note: b0 has NO impact on ODE;
x = sqrt(3); n = -1/5; b0 = 4/5; bf = sqrt(2);
params = list(x=x, n=n, b0=b0, bf=bf);
e = expression(log(exp(x^n) + b0) + bf*x)[[1]];
f = bf*x; df = bf; d2f = 0;
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### D =>
(exp(x^n)+b0)*dy - n*x^(n-1)*exp(x^n) - df*(exp(x^n)+b0) # = 0
# exp(x^n) = - b0*(dy - df) / (dy - n*x^(n-1) - df);

### D2 =>
(exp(x^n)+b0)*d2y + n*x^(n-1)*exp(x^n)*dy +
	- n*(n-1)*x^(n-2)*exp(x^n) - n^2*x^(2*n-2)*exp(x^n) +
	- d2f*(exp(x^n)+b0) - n*x^(n-1)*df*exp(x^n) # = 0
(d2y + n*x^(n-1)*dy +
		- n*(n-1)*x^(n-2) - n^2*x^(2*n-2) - n*x^(n-1)*df - d2f) * exp(x^n) +
	+ b0*d2y - d2f*b0 # = 0
# Subst =>
b0*(d2y + n*x^(n-1)*dy +
		- n*(n-1)*x^(n-2) - n^2*x^(2*n-2) - n*x^(n-1)*df - d2f) * (dy - df) +
	- (b0*d2y - b0*d2f)*(dy - n*x^(n-1) - df) # = 0
b0*n*x^(n-1)*d2y + b0*n*x^(n-1)*dy^2 - 2*b0*n*x^(n-1)*df*dy +
	- b0*(n*(n-1)*x^(n-2) + n^2*x^(2*n-2)) * (dy - df) +
	+ b0*n*x^(n-1)*df^2 - b0*n*x^(n-1)*d2f # = 0

### ODE:
x*d2y + x*dy^2 - (2*x*df + n*x^n + (n-1))*dy +
	+ (n*x^n + (n-1))*df + x*df^2 - x*d2f # = 0


####################

### y = x * log(exp(x^n) + b0) + f0

# Note: b0 has NO impact on ODE;
x = sqrt(3); n = -1/5; b0 = 4/5; bf = sqrt(2);
params = list(x=x, n=n, b0=b0, bf=bf);
e = expression(x * log(exp(x^n) + b0) + bf*x)[[1]];
f = bf*x; df = bf; d2f = 0;
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

### D =>
x*(exp(x^n)+b0)*dy - (exp(x^n)+b0)*(y-f) - n*x^(n+1)*exp(x^n) - x*df*(exp(x^n)+b0) # = 0
# exp(x^n) = - b0*(x*dy - (y-f) - x*df) / (x*dy - (y-f) - n*x^(n+1) - x*df);

### D2 =>
# TODO


#######################
#######################

### y = log(log(P(x))) + F0(x)

### y = log(log(x^2 + k)) + f

### D =>
(x^2+k)*log(x^2 + k)*dy - (x^2+k)*df*log(x^2 + k) - 2*x # = 0
# (x^2+k)*log(x^2 + k) = 2*x / (dy - df);

### D2 =>
(x^2+k)*log(x^2 + k)*d2y + 2*x*dy + 2*x*log(x^2 + k)*dy +
	- 2*x*(x^2+k)*df - 2*x*df*log(x^2 + k) - (x^2+k)*d2f*log(x^2 + k) - 2 # = 0
2*x*(x^2+k)/(dy - df) * d2y + 2*x*(x^2+k)*dy + 4*x^2/(dy - df) * dy +
	- 2*x*(x^2+k)^2*df - 4*x^2*df/(dy - df) - 2*x*(x^2+k)*d2f/(dy - df) - 2*(x^2+k) # = 0

# TODO

