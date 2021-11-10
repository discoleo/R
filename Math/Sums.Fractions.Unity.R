########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Roots of Unity
###
### draft v.0.1f


### Infinite Sums
### A.) Sum( (-1)^n / (k1*n + k0))
### B.) Based on I(x^k * log(x^n + 1))

### Theory:
# - can be transformed into:
#   Integral( x^(k0 - 1) / (x^n + 1) );


### Integral( P(x) / (x^n + 1) ):
# - see file:
#   Integrals.Fractions.Unity;


#######################

### Helper Functions

# required:
source("Polynomials.Helper.R")


### Other

sum.frn = function(x, n, k0=1) {
	x0 = if(k0 == 1) 1 else x^(k0-1);
	x0 / (x^n + 1)
}
# converges slowly:
sum.basicFr = function(n, iter=8000, k0=1) {
	sign.fr = rep(c(1,-1), iter %/% 2);
	if(iter %% 2 == 1) sign.fr = c(sign.fr, 1)
	1/k0 - sum(sign.fr / (n*seq(iter) + k0))
}
coeffs.frn = function(n=5) {
	# Roots of unity
	len = (n-1)/2; # only half the roots of unity
	# Simplification:
	m.sum = 2*cos(2*pi*seq(len) / n);
	m.shift = m.sum/2;
	# Coefficients
	b0 = 1/n; b = -2*b0;
	a = b0 * m.sum;
	D = b + a*m.shift;
	m.sq = sqrt(1 - m.shift^2);
	coeffs = list(a=a, b0=b0, b=b, D=D, m.sum=m.sum, m.shift=m.shift, m.sq=m.sq);
	return(coeffs)
}

### Exact formula: k0 = 1
sumExact.frn = function(n, lower=0, upper=1) {
	# only odd n;
	# for even n: see Integrals.Fractions.Unity;
	cf = coeffs.frn(n=n)
	a = cf$a; b0 = cf$b0; b = cf$b; # D = cf$D;
	m.sum = cf$m.sum; m.shift = cf$m.shift; m.sq = cf$m.sq;
	D = b + a*m.shift;
	r = b0*log(2) + sum( a/2 * log(upper^2 + 1 + m.sum*upper) ) +
		- sum(D / m.sq * atan((upper + m.shift)/m.sq)) +
		+ sum(D / m.sq * atan((lower + m.shift)/m.sq));
	### Note:
	# sum( a/2 * log(lower^2 + 1 + m.sum*lower) ) == 0 for lower = 0!
	if(lower != 0) r = r - sum( a/2 * log(lower^2 + 1 + m.sum*lower) );
	return(r);
}


########################
########################

### Section A:
### Simple Sums

### Sum( (-1)^n / (5*n + k0) )

### Integral:
# I = I[from=0, to=1] (-1)^n * x^(k0 - 1) * sum( x^(5*n) ) dx
# for k0 = 1 =>
# I = Sum( (-1)^n / (5*n + k0) )


### Examples:

n = 5;
# - any odd n!
# - formula is slightly different for even n;

### Sum vs Integral:
integrate(sum.frn, lower=0, upper=1, n=n)
sum.basicFr(n)
### Exact formula:
sumExact.frn(n)


### Ex 2:
n = 9
integrate(sum.frn, lower=0, upper=1, n=n)
sum.basicFr(n)
### Exact formula:
sumExact.frn(n)


#################

### Larger k0:
# Sum( (-1)^n / (5*n + k0) )
# where k0 = 2
n = 9
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b;
m.sum = cf$m.sum; m.shift = cf$m.shift; m.sq = cf$m.sq;

# Formula differs based on k0
k0 = 2 # fixed
#
integrate(sum.frn, lower=1E-8, upper=1, n=n, k0=k0)
sum.basicFr(n, k0=k0)

### Exact formula:
b0 + sum(a) - b0*log(2) +
	- sum( (a*m.sum + b)/2 * log(1 + 1 + m.sum) ) +
	+ sum( ((a*m.sum + b)*m.sum/2 - a) / m.sq * atan((1 + m.shift)/m.sq) ) +
	- sum( ((a*m.sum + b)*m.sum/2 - a) / m.sq * atan((0 + m.shift)/m.sq) )



###
k0 = 3
#
integrate(sum.frn, lower=1E-8, upper=1, n=n, k0=k0)
sum.basicFr(n, k0=k0)

# TODO

##################

### Test
### Fraction decompositions:

n = 5;
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b; D = cf$D;
m.sum = cf$m.sum; m.shift = cf$m.shift; m.sq = cf$m.sq;

# only a test value
x = 3;

###########
### k0 = 1:
1/(x^n + 1) # ==
b0/(x + 1) + sum( (a*x - b) / (x^2 + m.sum*x + 1) )
# decomposed into: LOG + ATAN components
b0/(x + 1) + sum( a/2 * (2*x + m.sum) / (x^2 + m.sum*x + 1) ) +
	- sum( D / ((x + m.shift)^2 + m.sq^2))

###########
### k0 = 2:
x/(x^n + 1) # ==
b0 + sum(a) - b0/(x + 1) +
	- sum( ((a*m.sum + b)*x + a) / (x^2 + m.sum*x + 1) )
# decomposed into: LOG + ATAN components
b0 + sum(a) - b0/(x + 1) +
	- sum( ((a*m.sum + b)/2*(2*x + m.sum)) / (x^2 + m.sum*x + 1) ) +
	+ sum( ((a*m.sum + b)*m.sum/2 - a) / (x^2 + m.sum*x + 1) )


####################
####################

### I(log(x^n + 1))

# sum(1/(n+1) - 1/(2*(2*n+1)) + 1/(3*(3*n+1)) - 1/(4*(4*n+1)) + ...)

### Integration by parts:
# x*log(x^n + 1) - n*x + n*I( 1 / (x^n + 1) )
# - for exact formula of the integral, see file:
#   Integrals.Fractions.Unity;


### Generalization:
### I(x^k * log(x^n + 1))

# sum(1/(n+1) - 1/(2*(2*n + k + 1)) + 1/(3*(3*n + k + 1)) - 1/(4*(4*n + k + 1)) + ...)

### Integration by parts:
# (k+1) * I(...) =>
# x^(k+1)*log(x^n + 1) - n*x^(k+1) / (k+1) + n*I( x^k / (x^n + 1) )
# - for exact formula of the integral, see file:
#   Integrals.Fractions.Unity;

intLog = function(n, k=0, lower=0, upper=1) {
	integrate(function(x) log(x^n + 1) * (if(k == 0) 1 else x^k), lower=lower, upper=upper);
}
intFr = function(n, k=0, lower=0, upper=1) {
	f = function(x) x^(k+1) * log(x^n + 1) - n*x^(k+1) / (k+1);
	r = f(upper) - f(lower) +
		+ n * integrate(function(x) (if(k==0) 1 else x^k) / (x^n + 1), lower=lower, upper=upper)$value;
	r = r / (k+1);
	return(r)
}
sumLogExp = function(n, k=0, x=1, iter=1000) {
	sign = rep(c(1,-1), iter %/% 2);
	if(iter %% 2 == 1) sign = c(sign, 1);
	xn = x^n;
	sign = if(k == 0) sign * x else sign * x^(k+1);
	sum( sign * xn^(seq(iter)) / (seq(iter)*(n*seq(iter) + k + 1)) ) 
}

###
n = 2
intLog(n=n)
intFr(n=n)
sumLogExp(n)

###
n = 5
intLog(n=n)
intFr(n=n)
sumLogExp(n)

###
n = 5; k = 2;
intLog(n=n, k=k)
intFr(n=n, k=k)
sumLogExp(n, k=k)

###
n = 5; k = 2; x0 = 0.75
intLog(n=n, k=k, upper=x0)
intFr(n=n, k=k, upper=x0)
sumLogExp(n, k=k, x=x0)


#################
#################


# library(pracma)

source("Polynomials.Helper.R")

### Cancellation using Roots of Unity:
### I[0, 1] + I[0, m] + I[0, m^2]
lineI = function(n, lower=0, upper=1) {
	line_integral(function(x)  log(x^n + 1), c(lower, upper));
}

sumLogExpM = function(n, m, k=0, x=1, iter=1000) {
	sign = rep(c(1,-1), iter %/% 2);
	if(iter %% 2 == 1) sign = c(sign, 1);
	xn = x^n;
	sign = if(k == 0) sign * x else sign * x^(k+1);
	sign[(seq(iter) %% m) != 1] = 0; # TODO: with k;
	sum( sign * xn^(seq(iter)) / (seq(iter)*(n*seq(iter) + k + 1)) ) 
}

### Examples:


n = 5
intLog(n=n)
lineI(n=n)

nm = 3 # != n;
m = unity(nm, all=F)
lineI(n=n) + lineI(n=n, upper=m) + lineI(n=n, upper=m^2)
sumLogExpM(n, m=nm) * 3

