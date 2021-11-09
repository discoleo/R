########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Roots of Unity


### Infinite Sums
### Sum( 1 / (k1*n + k0))

### Integral( 1 / (x^n + 1) ):
# - see file:
#   Integrals.Fractions.Unity;


#######################

### Helper Functions

# required:
source("Polynomials.Helper.R")


sum.frn = function(x, n) {
	1 / (x^n + 1)
}
# converges slowly:
sum.basicFr = function(n, iter=8000, k0=1) {
	sign.fr = rep(c(1,-1), iter %/% 2);
	if(iter %% 2 == 1) sign.fr = c(sign.fr, 1)
	1/k0 - sum(sign.fr / (n*seq(iter) + k0))
}


########################
########################


### Sum( 1 / (5*n + 1)))

### Integral:
# I[0, 1] sum( x^(5*n) ) dx

integrate(sum.frn, lower=0, upper=1, n=5)
sum.basicFr(5)

### Formula:

n = 5;
# Roots of unity
m = unity(n, all=TRUE)
m = m[-1] # all roots of unity (without 1)
len = (n-1)/2; # only half the roots of unity
m.conj = m[1:len]
m.conj = cbind(m.conj, 1/m.conj)
# Coefficients
b0 = 1/n; b = -2*b0;
m.sum = (m.conj[,1] + m.conj[,2]);
a = b0 * m.sum; m.shift = m.sum/2;
D = b + a*m.shift; m.sq = sqrt(1 - m.shift^2);

### Exact formula:
b0*log(2) + sum( a/2 * log((1 + m.conj[,1]) * (1 + m.conj[,2])) ) +
	- sum(D / m.sq * atan((1 + m.shift)/m.sq)) +
	+ sum(D / m.sq * atan((0 + m.shift)/m.sq));

### Note:
# sum( a/2 * log(m.conj[,1] * m.conj[,2]) ) == 0 !
	

### Test decomposition:
x = 3; # only a test value
1/(x^n + 1) # ==
b0/(x + 1) + sum( (a*x - b) / ((x + m.conj[,1]) * (x + m.conj[,2])) )
# decomposed into: LOG + ATAN components
b0/(x + 1) + sum( a/2*(2*x + m.sum) / ((x + m.conj[,1]) * (x + m.conj[,2])) ) +
	- sum( D / ((x + m.shift)^2 + m.sq^2))


