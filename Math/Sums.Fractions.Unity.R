########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Roots of Unity
###
### draft v.0.1j


### Infinite Sums
### A.) Sum( (-1)^n / (k1*n + k0))
### B.) Based on I(x^k0 * log(x^n + 1))
### C.) Term-cancellation using Roots of unity

### Theory:
# - can be transformed into:
#   Integral( x^(k0 - 1) / (x^n + 1) );
# - exact formulas based on exact formula
#   for this integral;


### Integral( P(x) / (x^n + 1) ):
# - see file:
#   Integrals.Fractions.Unity;


#######################

### Helper Functions

# required:
source("Polynomials.Helper.R")


### Other

# Fraction: used for integration
sum.frn = function(x, n, k0=1) {
	x0 = if(k0 == 1) 1 else x^(k0-1);
	x0 / (x^n + 1)
}
# Sum: converges slowly:
sum.basicFr = function(n, iter=8000, k0=1) {
	sign.fr = rep(c(1,-1), iter %/% 2);
	if(iter %% 2 == 1) sign.fr = c(sign.fr, 1)
	1/k0 - sum(sign.fr / (n*seq(iter) + k0))
}
# Exact Integral: Fraction decomposition
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

### Sum( (-1)^j / (k1*j + k0) )

### Integral:
# I = I[from=0, to=1] x^(k0 - 1) * sum( (-1)^j * x^(k1*j) ) dx
# Expansion =>
# I = Sum( (-1)^j / (k1*j + k0) )
# j = from 1 to Inf;


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


##############
### Larger k0:
# Sum( (-1)^j / (k1*j + k0) )

# any odd integer
n = 9
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b;
m.sum = cf$m.sum; m.shift = cf$m.shift; m.sq = cf$m.sq;

# Formula differs based on k0

###########
### k0 == 2

k0 = 2 # fixed
cP0 = a;
cP1 = (a*m.sum + b) / 2;
cP2 = ((a*m.sum + b)*m.sum/2 - a);
#
integrate(sum.frn, lower=0, upper=1, n=n, k0=k0)
sum.basicFr(n, k0=k0)

### Exact formula:
b0 - b0*log(2) + sum(cP0) +
	- sum( cP1 * log(1 + 1 + m.sum) ) +
	+ sum( cP2 / m.sq * atan((1 + m.shift)/m.sq) ) +
	- sum( cP2 / m.sq * atan((0 + m.shift)/m.sq) )


###########
### k0 == 3

k0 = 3
cP0 = a/2 - b - a*m.sum;
cP1 = (- a + b*m.sum + a*m.sum^2) / 2;
cP2 = (2*b + 3*a*m.sum - b*m.sum^2 - a*m.sum^3) / 2;
#
integrate(sum.frn, lower=0, upper=1, n=n, k0=k0)
sum.basicFr(n, k0=k0)

### Exact formula:
-1/2 * b0 + b0*log(2) + sum(cP0) +
	+ sum(cP1 * log(1 + 1 + m.sum)) +
	+ sum( cP2 / m.sq * atan((1 + m.shift)/m.sq) ) +
	- sum( cP2 / m.sq * atan((0 + m.shift)/m.sq) )


###########
### k0 == 4

k0 = 4
cP0 = a/3 - (b + a*m.sum)/2 + (a*m.sum^2 - a + b*m.sum);
cP1 = (b + 2*a*m.sum - b*m.sum^2 - a*m.sum^3) / 2;
cP2 = (2*a - 3*b*m.sum - 4*a*m.sum^2 + b*m.sum^3 + a*m.sum^4) / 2;
#
integrate(sum.frn, lower=0, upper=1, n=n, k0=k0)
sum.basicFr(n, k0=k0)

### Exact formula:
5/6 * b0 - b0*log(2) + sum(cP0) +
	+ sum(cP1 * log(1 + 1 + m.sum)) +
	+ sum( cP2 / m.sq * atan((1 + m.shift)/m.sq) ) +
	- sum( cP2 / m.sq * atan((0 + m.shift)/m.sq) )


### Derivation:
k0 = -1 + 4
p1 = toPoly.pm("x^k0 * (a*x - b)")
p2 = toPoly.pm("x^2 + m.sum*x + 1")
pR = split.pm.fraction(p1, p2, by="x")
pR$D$coeff = 2* pR$D$coeff # Div by 2!
pR$Ct$coeff = 2* pR$Ct$coeff # Div by 2!

cP0 = x^2*a - (b + a*m.sum)*x + (a*m.sum^2 - a + b*m.sum);


##################

### Even n:
# n = 2*(2*m + 1);
n = 6; # 10, 14, 18, ...;
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b; ai = 1i * a;
m.sum = cf$m.sum; m.shift = 1i*m.sum/2;
D = b + a/2*m.sum; m.sq = sqrt(1 + m.shift^2);

###########
### k0 == 1
# [fixed]

integrate(sum.frn, lower=0, upper=1, n=n)
sum.basicFr(n)
# Exact formula:
2*b0*pi/4 + sum( ai/2 * log(-1i*m.sum) ) +
	+ sum( D / (2*m.sq) * log((1 + m.shift - m.sq)/(1 + m.shift + m.sq)) ) +
	+ sum( D / (2*m.sq) * log((0 + m.shift + m.sq)/(0 + m.shift - m.sq)) );


##################
##################

### Test
### Fraction decompositions:

### Odd n:

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


###########
### Even n:

# TODO:
# - resolve issues for n = 4*k;
# - the "i"-trick works only for n = 2*(2*k+1);
# TODO: full root of "-1" of order n;

n = 6; # 10, 14, 18
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b; ai = 1i * a;
m.sum = cf$m.sum; m.shift = 1i*m.sum/2;
D = b + a/2*m.sum; m.sq = sqrt(1 + m.shift^2);

# only a test value
x = 2.1

1/(x^n + 1) # ==
2*b0/(x^2 + 1) + sum( (ai*x + b) / (x^2 + 1i*m.sum*x - 1) )
# decomposed into: LOG + "ATAN" components
2*b0/(x^2 + 1) + sum( ai/2 * (2*x + 1i*m.sum) / (x^2 + 1i*m.sum*x - 1) ) +
	+ sum( D / ((x + m.shift)^2 - m.sq^2))


#######################
#######################

#######################
### I(log(x^n + 1)) ###
#######################

# sum(1/(n+1) - 1/(2*(2*n+1)) + 1/(3*(3*n+1)) - 1/(4*(4*n+1)) + ...)

### Integration by parts:
# x*log(x^n + 1) - n*x + n*I( 1 / (x^n + 1) )
# - for exact formula of the integral, see file:
#   Integrals.Fractions.Unity;


### Generalization:
### I(x^k * log(x^n + 1))

# both Harmonic & non-Harmonic:
# sum(1/(n+k+1) - 1/(2*(2*n + k + 1)) + 1/(3*(3*n + k + 1)) - 1/(4*(4*n + k + 1)) + ...)
# sum(1/(n+k+1) + 1/(2*(2*n + k + 1)) + 1/(3*(3*n + k + 1)) + 1/(4*(4*n + k + 1)) + ...)

### Integration by parts:
# (k+1) * I(...) =>
# x^(k+1)*log(x^n + 1) - n*x^(k+1) / (k+1) + n*I( x^k / (x^n + 1) )
# - for exact formula of the integral, see file:
#   Integrals.Fractions.Unity;

intLog = function(n, k=0, lower=0, upper=1, harmonic=TRUE) {
	f = if(harmonic) function(x) log(x^n + 1) * (if(k == 0) 1 else x^k)
		else function(x) - log(-x^n + 1) * (if(k == 0) 1 else x^k)
	integrate(f, lower=lower, upper=upper);
}
intFr = function(n, k=0, lower=0, upper=1, harmonic=TRUE) {
	sign = if(harmonic) 1 else -1;
	f = function(x) x^(k+1) * log(sign * x^n + 1) - n*x^(k+1) / (k+1);
	fdf = function(x) (if(k==0) 1 else x^k) / (sign * x^n + 1);
	r = f(upper) - f(lower) +
		+ n * integrate(fdf, lower=lower, upper=upper)$value;
	r = sign * r / (k+1);
	return(r)
}
sumLogExp = function(n, k=0, x=1, harmonic=TRUE, iter=1000) {
	sign = if( ! harmonic) 1 else rep(c(1,-1), iter %/% 2);
	if(harmonic && iter %% 2 == 1) sign = c(sign, 1);
	xn = x^n;
	sign = if(k == 0) sign * x else sign * x^(k+1);
	sum( sign * xn^(seq(iter)) / (seq(iter)*(n*seq(iter) + k + 1)) ) 
}
terms.sum = function(k1, k0=0, n=20, harmonic=TRUE) {
	sign = if( ! harmonic) "+" else rep(c("+","-"), n %/% 2);
	if(harmonic && n %% 2 == 1) sign = c(sign, "+");
	m1 = seq(n); m2 = k1*m1 + k0 + 1;
	m = data.frame(sign, m1, m2);
	return(m)
}
cat.sum = function(x, nr.split=6, sep.add=" + ") {
	isSplit = FALSE;
	if(nrow(x) > nr.split) {
		len = nrow(x); isSplit = TRUE;
		lenHalf = len %/% 2;
		lenH1 = lenHalf; x.tmp = x[seq(lenH1 + 1, len), ];
		if(len %% 2 == 1) {
			lenH1 = lenH1 + 1;
			x.tmp = rbind(x.tmp, c("", 0, 0));
		}
		x = cbind(x[seq(lenH1), ], x.tmp);
	}
	x = sapply(x, function(x) format(x));
	sep = c(" 1 / (", " * ");
	sep.add = paste0(")", sep.add, "\n")
	sep = if(isSplit) c(sep, ")   ", sep, sep.add) else c(sep, sep.add);
	cat(t(x), sep=sep); cat(")\n");
	invisible()
}

### Examples:

###
n = 2
intLog(n=n)
intFr(n=n)
sumLogExp(n)
cat.sum(terms.sum(n, k0=0))

###
n = 5
intLog(n=n)
intFr(n=n)
sumLogExp(n)
cat.sum(terms.sum(n, k0=0))

###
n = 5; k = 2;
intLog(n=n, k=k)
intFr(n=n, k=k)
sumLogExp(n, k=k)

### Non-Harmonic
n = 5; k = 2;
intLog(n=n, k=k, harmonic=FALSE)
intFr(n=n, k=k, upper = 1 - 1E-7, harmonic=FALSE)
sumLogExp(n, k=k, harmonic=FALSE)

### Test formulas:
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


#########################
#########################

