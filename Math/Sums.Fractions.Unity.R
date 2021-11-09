########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Roots of Unity
###
### draft v.0.1c-fix


### Infinite Sums
### Sum( (-1)^n / (k1*n + k0))

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
	m = unity(n, all=TRUE)
	m = m[-1] # all roots of unity (without 1)
	len = (n-1)/2; # only half the roots of unity
	m.conj = m[1:len];
	m.conj = cbind(m.conj, 1/m.conj);
	m.sum = (m.conj[,1] + m.conj[,2]);
	# Coefficients
	b0 = 1/n; b = -2*b0;
	a = b0 * m.sum;
	m.shift = m.sum/2; m.sq = sqrt(1 - m.shift^2);
	D = b + a*m.shift;
	coeffs = list(a=a, b0=b0, b=b, D=D, m.sum=m.sum, m.shift=m.shift, m.sq=m.sq);
	return(coeffs)
}


########################
########################


### Sum( (-1)^n / (5*n + k0) )

### Integral:
# k0 = 1 =>
# I[0, 1] sum( x^(5*n) ) dx


### Examples:

n = 5;
# Roots of unity
cf = coeffs.frn(n=n)
a = cf$a; b0 = cf$b0; b = cf$b; D = cf$D;
m.sum = cf$m.sum; m.shift = cf$m.shift; m.sq = cf$m.sq;

### Sum vs Integral:
integrate(sum.frn, lower=0, upper=1, n=n)
sum.basicFr(n)

### Exact formula:
b0*log(2) + sum( a/2 * log(1 + 1 + m.sum) ) +
	- sum(D / m.sq * atan((1 + m.shift)/m.sq)) +
	+ sum(D / m.sq * atan((0 + m.shift)/m.sq));

### Note:
# sum( a/2 * log(m.conj[,1] * m.conj[,2]) ) == 0 !


#################

### Sum( (-1)^n / (5*n + k0) )
n = 5
k0 = 2
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
integrate(sum.frn, lower=1E-8, upper=1, n=5, k0=k0)
sum.basicFr(5, k0=k0)

# TODO

##################

### Test
### Fraction decompositions:

# only a test value
x = 3;

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

