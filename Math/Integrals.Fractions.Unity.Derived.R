
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### of Polynomial Fractions
### derived from Roots of Unity
###
### draft v.0.1e-tan



### Pre-Requisits

### Roots of Unity:
# - Base-Integrals:
#   x^k / (x^n - 1) dx;
#   x^k / (x^n + 1) dx;
#   see file: Integrals.Fractions.Unity.R;
# - Powers of Roots of Unity:
#   x^k / (x^n - 1)^p dx, x^k / (x^n + 1)^p dx
#   see file: Integrals.Fractions.Unity.Powers.R;


###############
### History ###

# draft v.0.1e-tan:
# - added some trigonometric derivatives
#   (based actually on Cardan-polynomials);
# draft v.0.1e:
# - various radicals:
#   sqrt(x + s) / (x^n - 1);
#   (x + s)^(1/p) / (x^n - 1);
# draft v.0.1d - v.0.1d-std:
# - started trigonometric derivatives:
#   cos(x)^(2*n-2) / (1 - cos(x)^(2*n));
# - more work on polynomial/radical derivatives:
#   x^(1/p) / ((x + s)^n - 1);
# - standardize shift notation; [v.0.1d-bis & v.0.1d-std]
# draft v.0.1c - v.0.1c-sqrt:
# - added/documented some derived polynomials;
# - fraction decomposition for:
#   1 / ((x^2 + s)^n - 1);
#   1 / ((x^p + s)^n - 1); [initial thoughts: v.0.1c-bis]
#   sqrt(x) / ((x + s)^n - 1); [initial thoughts: v.0.1c-sqrt]
# draft v.0.1b:
# - added more details to the examples with fractional powers;
# draft v.0.1a:
# - moved examples from Integrals.Fractions.Unity.R
#   to this file;
# - added examples with fractional powers;


################
### Introduction


### A.) Fractional Powers
### B.) Trigonometric Derivations
### C.) Derived Polynomials

### Terminology

### I[k, n]
# I[k, n] = Integral x^k / (x^n + 1) dx;
# or
# I[k, n] = Integral x^k / (x^n - 1) dx;
# (depending on context)


#########################
### A.) Fractional Powers

### I[k/m, j/n]
# Integral x^(k/m) / (x^(j/n) - 1) dx
# =>
# y = x^(1/(m*n))
# x = y^(m*n)
# dx = m*n*y^(m*n - 1) * dy
# =>
# I[k/m, j/n] = I  m*n * y^(k*n + m*n - 1) / (x^(m*j) - 1) dy
# = m*n * I[k*n + m*n - 1, m*j];
# Note: limits of integration also change!


#################################
### B.) Trigonometric Derivations
# - TODO: document + expand;
# - see also file:
#   Integrals.Exercises.R;


###########################
### C.) Derived Polynomials
# - TODO

### C.1.) Div by (x^n + (1-x)^n)

### x^(n-2) / (x^n + (1-x)^n)
### x^(n-2) / (x^n - (1-x)^n)
# I[2*p+1, 2*n] =>
# y = 1 - 1 / (x^2 + 1)
# x^2 = 1 / (1 - y) - 1
# 2*x*dx = (1-y)^(-2) * dy;
# =>
# I[2*p+1, 2*n] = 1/2 * Int y^k / (1-y)^(p+2) / (y^n / (1-y)^n + 1) dy
#  = 1/2 * Int y^p * (1-y)^(n - p - 2) / (y^n + (1-y)^n) dy

### Examples
n = 5
lim = c(2, 3)
# Variant: "+"
integrate(function(x) 2 * x^(2*n-3)/(x^(2*n) + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^(n - 2)/(x^n + (1-x)^n), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# Variant: "-"
integrate(function(x) 2 * x^(2*n-3)/(x^(2*n) - 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^(n - 2)/(x^n - (1-x)^n), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))


### C.2.) Int x^k / ((x^2 + k)^n - 1)

### Fraction Decomposition:
# F[0] = 1 / ((x^2 + k)^n - 1)
# = b0 / (x^2 + k - 1) + sum( (a*(x^2 + k) + b) / ((x^2 + k - m^j)*(x^2 + k - m^(-j))) )
# [n = odd powers]
# where b0, a, b = coefficients as per:
# Integrals.Fractions.Unity.R;


######################
######################

### helper function
I.f = function(f, lim) {
	integrate(f, lower=lim[1], upper=lim[2])$value
}
I.pf = function(b, n, lim) {
	if(length(b) > 1 && length(n) > 1 && length(b) != length(n)) {
		stop("Differing lengths!")
	}
	pow = n + 1
	coeff = b / pow
	sum(coeff * lim[2]^pow) - sum(coeff * lim[1]^pow)
}
###
roots.conj = function(n) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	n_1 = n - 1; n.half = floor(n_1/2)
	i = 1:(n.half)
	m.h = m^(1:n.half)
	m.sum = m.h + 1/m.h
	m.h = cbind(m.h, 1/m.h)
	#
	r = list(m=m, m.sum=m.sum, m.half = m.h)
	return(r)
}
decompose.fr = function(n, type="minus") {
	# TODO: n = even powers, type = "plus";
	m = roots.conj(n)
	b0 = 1/n
	a = b0 * m$m.sum
	b = -2*b0
	return(list(b0=b0, a=a, b=b, m=m$m, m.sum=m$m.sum, m.half=m$m.half))
}

#####################

################
### Examples ###
################

### A.) Fractional Powers

lim = c(1.2, 1.5)
### Example: x^(3/5)/(x^(5/4) + 1)
integrate(function(x) x^(3/5)/(x^(5/4) + 1), lower=lim[1], upper=lim[2]) # ==
integrate(function(x) 20 * x^31/(x^25 + 1), lower=lim[1]^(1/20), upper=lim[2]^(1/20)) # ==
integrate(function(x) 20*x^6 - 20 * x^6/(x^25 + 1), lower=lim[1]^(1/20), upper=lim[2]^(1/20)) # ==
I.pf(20, 6, lim^(1/20)) - 20 * I.f(function(x) x^6/(x^25 + 1), lim^(1/20))
# I[6] = Integral x^6/(x^25 + 1) can be computed based on I[0] = 1/(x^25 + 1);
#        or using direct polynomial division (within the I[0] decomposition);
# I[0] => [1/x] => I[23] => [half] => I[11] => [half] => I[5]
#      => [half] => I2 => [1/x] => I[21] => [half] => I[10]
#      => [1/x] => I[13] => [half] => I[6];
# I[6] = f(I[0]) using the reverse sequence as per:
#   Integrals.Fractions.Unity.R;
# I[0] can be computed using the fraction decomposition as per:
#   Integrals.Fractions.Unity.R;


lim = c(1.2, 1.5)
### Example: (x^(3/5) + x^(1/3)) / (x^(5/4) + 1)
integrate(function(x) (x^(3/5) + x^(1/3)) / (x^(5/4) + 1), lower=lim[1], upper=lim[2]) # ==
I.f(function(x) 20 * x^31/(x^25 + 1), lim=lim^(1/20)) +
	I.f(function(x) 12 * x^15/(x^15 + 1), lim=lim^(1/12)) # ==
I.f(function(x) 20*x^6 - 20 * x^6/(x^25 + 1), lim=lim^(1/20)) +
	I.f(function(x) 12 - 12/(x^15 + 1), lim=lim^(1/12)) # ==
I.pf(20, 6, lim=lim^(1/20)) - 20*I.f(function(x) x^6/(x^25 + 1), lim=lim^(1/20)) +
	I.pf(12, 0, lim=lim^(1/12)) - 12*I.f(function(x) 1/(x^15 + 1), lim=lim^(1/12))
# I[6, 25] & I[0, 15] can be computed exactly as per:
# Integrals.Fractions.Unity.R;


###########################

###########################
### C.) Derived Polynomials


### Simple Fraction
# [Test]
x = 3
n = 5
# Test
r = decompose.fr(n)
#
1 / (x^n - 1)
r$b0/(x - 1) + sum( (r$a*x + r$b) / (x^2 - r$m.sum * x + 1) )
#
1 / (x^n + 1)
r$b0/(x + 1) + sum( (r$a*x - r$b) / (x^2 + r$m.sum * x + 1) )


### Derived Polynomials

### Fraction: 1 / ( (x^2 + s)^n - 1 )
x = 2
s = 1.5
n = 5
# Test
r = decompose.fr(n)
#
1 / ((x^2 + s)^n - 1)
r$b0/(x^2 + s - 1) + sum( (r$a*(x^2+s) + r$b) / ((x^2 + s)^2 - r$m.sum * (x^2 + s) + 1) )
r$b0/(x^2 + s - 1) + sum( (r$a*(x^2+s) + r$b) * (1/(x^2 + s - r$m.half[,1]) - 1/(x^2 + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )
r$b0/(x^2 + s - 1) + sum( ((r$a*r$m.half[,1] + r$b)/(x^2 + s - r$m.half[,1]) -
	(r$a*r$m.half[,2] + r$b)/(x^2 + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )


### Fraction: 1 / ( (x^p + s)^n - 1 )
x = 1.3
s = -1.1
n = 5
p = 4
# Test
r = decompose.fr(n)
#
1 / ((x^p + s)^n - 1)
r$b0/(x^p + s - 1) + sum( (r$a*(x^p+s) + r$b) / ((x^p + s)^2 - r$m.sum * (x^p + s) + 1) )
r$b0/(x^p + s - 1) + sum( (r$a*(x^p+s) + r$b) * (1/(x^p + s - r$m.half[,1]) - 1/(x^p + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )
r$b0/(x^p + s - 1) + sum( ((r$a*r$m.half[,1] + r$b)/(x^p + s - r$m.half[,1]) -
	(r$a*r$m.half[,2] + r$b)/(x^p + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )
# TODO: re-scaling to 1/(x^p + 1);
x0.s = if(p %% 2 == 0 && s-1 < 0) x / (s - 1 + 0i)^(1/p) else x / (s - 1)^(1/p) # TODO: robust root.f()
x.s = x / (s - r$m.half)^(1/p) # for n = 5: 2 columns;
r$b0/(s - 1) / (x0.s^p + 1) + sum( 1 / (r$m.half[,1] - r$m.half[,2]) *
	((r$a*r$m.half[,1] + r$b)/(s - r$m.half[,1]) / (x.s[,1]^p + 1) -
	(r$a*r$m.half[,2] + r$b)/(s - r$m.half[,2]) / (x.s[,2]^p + 1)) )


### Higher Terms

### Method 1:
# - direct polynomial division;
### Method 2:
# - TODO: search for compact/simplified formulas;

### Fraction: x^2 / ( (x^2 + s)^n - 1 )
### Fraction: sqrt(x) / ( (x + s)^n - 1 )
### Fraction: sqrt(x + s) / (x^n - 1)
###
### Int: x^(2*k + 1) / ( (x^2 + s)^n - 1 ) [trivial]
### Int: x^(2*k) / ( (x^2 + s)^n - 1 ) [Methods 1 or 2]
x = 2
s = 1.5
n = 5
lim = c(1.1, 3)
# Test
r = decompose.fr(n)
#
x^2 / ((x^2 + s)^n - 1)
r$b0 * x^2/(x^2 + s - 1) + sum( (r$a*(x^2+s) + r$b) * x^2 / ((x^2 + s)^2 - r$m.sum * (x^2+s) + 1) )
r$b0 * x^2/(x^2 + s - 1) + sum( (r$a*(x^2+s) + r$b) * x^2 * (1/(x^2 + s - r$m.half[,1]) - 1/(x^2 + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )
r$b0 * x^2/(x^2 + s - 1) + sum( (r$a*(x^2+s) + r$b) / (r$m.half[,2] - r$m.half[,1]) *
	((s - r$m.half[,1])/(x^2 + s - r$m.half[,1]) - (s - r$m.half[,2])/(x^2 + s - r$m.half[,2])) )
r$b0 - r$b0*(s-1)/(x^2 + s - 1) + sum( 1 / (r$m.half[,2] - r$m.half[,1]) *
	((s - r$m.half[,1])*(r$a + (r$b + r$a * r$m.half[,1])/(x^2 + s - r$m.half[,1])) -
	(s - r$m.half[,2])*(r$a + (r$b + r$a * r$m.half[,2])/(x^2 + s - r$m.half[,2]))) )
### sqrt(x)
integrate(function(x) sqrt(x) / ((x+s)^n - 1), lower=lim[1], upper=lim[2])
integrate(function(x) 2 * x^2 / ((x^2+s)^n - 1), lower=sqrt(lim[1]), upper=sqrt(lim[2]))
#
integrate(function(x) sqrt(x) / ((x+s)^n + 1), lower=lim[1], upper=lim[2])
integrate(function(x) 2 * x^2 / ((x^2+s)^n + 1), lower=sqrt(lim[1]), upper=sqrt(lim[2]))
### sqrt(x + s)
integrate(function(x) sqrt(x+s) / (x^n - 1), lower=lim[1], upper=lim[2])
integrate(function(x) sqrt(x) / ((x-s)^n - 1), lower=lim[1]+s, upper=lim[2]+s)
integrate(function(x) 2 * x^2 / ((x^2-s)^n - 1), lower=sqrt(lim[1] + s), upper=sqrt(lim[2] + s))
#
integrate(function(x) sqrt(x+s) / (x^n + 1), lower=lim[1], upper=lim[2])
integrate(function(x) sqrt(x) / ((x-s)^n + 1), lower=lim[1]+s, upper=lim[2]+s)
integrate(function(x) 2 * x^2 / ((x^2-s)^n + 1), lower=sqrt(lim[1] + s), upper=sqrt(lim[2] + s))
#
k = 3
integrate(function(x) 2 * x^k / ((x^2+s)^n + 1), lower=lim[1], upper=lim[2])
integrate(function(x) (x-s)^((k-1)/2) / (x^n + 1), lower=(lim[1])^2 + s, upper=(lim[2])^2 + s)


### Generalisation: x^p
### Fraction: x^p / ( (x^p + s)^n - 1 )
### Fraction: x^(1/p) / ( (x + s)^n - 1 )
### Fraction: (x+s)^(1/p) / (x^n - 1)
x = 2
s = 1.5
n = 5
p = 3
lim = c(1.1, 3)
# Test
r = decompose.fr(n)
#
x^p / ((x^p + s)^n - 1)
r$b0 * x^p/(x^p + s - 1) + sum( (r$a*(x^p+s) + r$b) * x^p / ((x^p + s)^2 - r$m.sum * (x^p + s) + 1) )
r$b0 * x^p/(x^p + s - 1) + sum( (r$a*(x^p+s) + r$b) * x^p * (1/(x^p + s - r$m.half[,1]) - 1/(x^p + s - r$m.half[,2])) / (r$m.half[,1] - r$m.half[,2]) )
r$b0 * x^p/(x^p + s - 1) + sum( (r$a*(x^p+s) + r$b) / (r$m.half[,2] - r$m.half[,1]) *
	((s - r$m.half[,1])/(x^p + s - r$m.half[,1]) - (s - r$m.half[,2])/(x^p + s - r$m.half[,2])) )
### x^(1/p)
integrate(function(x) x^(1/p) / ((x+s)^n - 1), lower=lim[1], upper=lim[2])
integrate(function(x) p * x^p / ((x^p+s)^n - 1), lower=(lim[1])^(1/p), upper=(lim[2])^(1/p))
#
integrate(function(x) x^(1/p) / ((x+s)^n + 1), lower=lim[1], upper=lim[2])
integrate(function(x) p * x^p / ((x^p+s)^n + 1), lower=(lim[1])^(1/p), upper=(lim[2])^(1/p))
### (x + s)^(1/p)
integrate(function(x) (x + s)^(1/p) / (x^n - 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^(1/p) / ((x-s)^n - 1), lower=lim[1] + s, upper=lim[2] + s)
integrate(function(x) p * x^p / ((x^p-s)^n - 1), lower=(lim[1] + s)^(1/p), upper=(lim[2] + s)^(1/p))
#
integrate(function(x) (x + s)^(1/p) / (x^n + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^(1/p) / ((x-s)^n + 1), lower=lim[1] + s, upper=lim[2] + s)
integrate(function(x) p * x^p / ((x^p-s)^n + 1), lower=(lim[1] + s)^(1/p), upper=(lim[2] + s)^(1/p))


######################
######################

### Trigonometric

### cos(x)^(2*n-2) / (1 - cos(x)^(2*n))
### (1 + cos(x)^(2*n) * sin(x)^2) / (1 - cos(x)^(2*n))
#
# from Int 1 / ( (x^2 + 1)^n - 1 )
# x = tan(y)
# dx = 1/cos(y)^2 dy
# =>
# Int cos(x)^(2*n-2) / (1 - cos(x)^(2*n)) dx

n = 5
lim = c(1, 2)
#
integrate(function(x) 1 / ((x^2 + 1)^n - 1), lower=lim[1], upper=lim[2])
integrate(function(x) cos(x)^(2*n-2) / (1 - cos(x)^(2*n)), lower=atan(lim[1]), upper=atan(lim[2]))
# TODO: 1 / (...)
integrate(function(x) 1 / ((x^2 + 1)^n - 1), lower=lim[1], upper=lim[2])
I.pf(-1,0, lim=atan(lim)) + I.f(function(x) (1 + cos(x)^(2*n-2) * sin(x)^2) / (1 - cos(x)^(2*n)), lim=atan(lim))



######################
######################

### Other

### n = 5
# (and some n = 10 equivalent to n = 5)

lim = c(2, 3)
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 2 * x^7/(x^10 + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^3/(x^5 + 1), lower=lim[1]^2, upper=lim[2]^2)
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) 2 * sin(x)^7*cos(x)/(sin(x)^10 + cos(x)^10), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) x^3/(x^5 + (1-x)^5), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) x^3/(x^5 + (1-x)^5), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) (x - 1)^3/(x^5 + (1-x)^5), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) (x - 1)^3/(x^5 + (1-x)^5), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^3 / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) 2*(1+x)^3/((1+x)^5 + (1-x)^5), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))
integrate(function(x) (1+x)^3/(5*x^4 + 10*x^2 + 1), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))

# derived: (x^3 - (x-1)^3) / Q(x) = (3*x^2 - 3*x + 1) / Q(x)
# the 2nd Integral can be easily computed exactly;
I.f(function(x) (3*x^2 - 3*x + 1)/(x^5 + (1-x)^5), lim) # ==
I.f(function(x) x^3/(x^5 + 1), lim = c(1/(1-lim[1]) - 1, 1/(1-lim[2]) - 1)) -
	I.f(function(x) x^3/(x^5 + 1), lim = c(1/lim[1] - 1, 1/lim[2] - 1))



### n == 6 & 12
# paramater: can be changed;
n = 12
#
lim = c(2, 3)
#
n.half = n/2
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 1/2 * x^(n.half-2)/(x^n.half + 1), lower=lim[1]^2, upper=lim[2]^2)
integrate(function(x) x^(n-3)/(x^n + 1), lower=lim[1], upper=lim[2])
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) sin(x)^(n-3)*cos(x)/(sin(x)^n + cos(x)^n), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) -1/2 * (1 - x)^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) -1/2 * (1 - x)^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^(n.half - 2) / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) (1+x)^(n.half - 2)/((1+x)^n.half + (1-x)^n.half), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))



### n = 7
# (and some n = 14 equivalent to n = 7)
# see above for parametric version;

lim = c(2, 3)
# trivially computable: see Integrals.Fractions.Unity.R;
integrate(function(x) 2 * x^11/(x^14 + 1), lower=lim[1], upper=lim[2])
integrate(function(x) x^5/(x^7 + 1), lower=lim[1]^2, upper=lim[2]^2)
### derived Integrals: *ALL* are equal!
# tan(y) = x => dx = 1/cos(y)^2 dy;
integrate(function(x) 2 * sin(x)^11*cos(x)/(sin(x)^14 + cos(x)^14), lower=atan(lim[1]), upper=atan(lim[2]))
# y = sin(x)^2 => dy = 2*sin(x)*cos(x)*dx;
integrate(function(x) x^5/(x^7 + (1-x)^7), lower=sin(atan(lim[1]))^2, upper=sin(atan(lim[2]))^2)
integrate(function(x) x^5/(x^7 + (1-x)^7), lower=1-1/(lim[1]^2+1), upper=1-1/(lim[2]^2+1))
# y = cos(x)^2 => dy = -2*sin(x)*cos(x)*dx
integrate(function(x) (x - 1)^5/(x^7 + (1-x)^7), lower=cos(atan(lim[1]))^2, upper=cos(atan(lim[2]))^2)
integrate(function(x) (x - 1)^5/(x^7 + (1-x)^7), lower=1/(lim[1]^2+1), upper=1/(lim[2]^2+1))
# from x^5 / Q(x): x = (y + 1)/2 => dx = 1/2 * dy
integrate(function(x) 2*(1+x)^5/((1+x)^7 + (1-x)^7), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))
integrate(function(x) (1+x)^5/(7*x^6 + 35*x^4 + 21*x^2 + 1), lower=1-2/(lim[1]^2+1), upper=1-2/(lim[2]^2+1))

# derived: (x^5 - (x-1)^5) / Q(x) = (5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1) / Q(x)
# the 2nd Integral can be easily computed exactly;
I.f(function(x) (5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1)/(x^7 + (1-x)^7), lim) # ==
I.f(function(x) x^5/(x^7 + 1), lim = c(1/(1-lim[1]) - 1, 1/(1-lim[2]) - 1)) -
	I.f(function(x) x^5/(x^7 + 1), lim = c(1/lim[1] - 1, 1/lim[2] - 1))


###########

### [old]
# TODO: cleanup

n = 11/2
#
n.half = n/2
integrate(function(x) x^(n-3)/(x^n + 1), lower=lower, upper=upper)
integrate(function(x) sin(x)^(n-3)*cos(x)/(sin(x)^n + cos(x)^n), lower=atan(lower), upper=atan(upper))
integrate(function(x) 1/2 * x^(n.half - 2)/(x^n.half + (1-x)^n.half), lower=sin(atan(lower))^2, upper=sin(atan(upper))^2)

###########
###########

lower = 0.2 # 0.6
upper = 0.8
integrate(function(x) 1/2 * x^3/(x^5 + 1), lower=2/(1-lower)-1, upper=2/(1-upper)-1 )
integrate(function(x) (1+x)^3/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)

### x^3 + x
log( ((1+upper)^5 + (1-upper)^5) / ((1+lower)^5 + (1-lower)^5) ) * 1/40
integrate(function(x) (x^3 + x)/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)

### only x
(atan((upper^2+1)*sqrt(5)/2i) - atan((lower^2+1)*sqrt(5)/2i)) * sqrt(5)/2i / 20
integrate(function(x) x/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)
# => x^3, 3*x^2 + 1

### 5*x^4 + 10*x^2 + 1
(upper - lower)/2
integrate(function(x) (5*x^4 + 10*x^2 + 1)/((1+x)^5 + (1-x)^5), lower=lower, upper=upper)


############
############

n = 5
lower = 0.2
upper = 0.8
integrate(function(x) 1/(x^(n/2) - 1), lower=lower, upper=upper)
integrate(function(x) 2 * x/(x^n - 1), lower=sqrt(lower), upper=sqrt(upper))
integrate(function(x) (x^(n/2) + 1)/(x^n - 1), lower=lower, upper=upper)


integrate(function(x) 1/(x^(n/3) - 1), lower=lower, upper=upper)
integrate(function(x) 3 * x^2/(x^n - 1), lower=lower^(1/3), upper=upper^(1/3))


integrate(function(x) x^(1/2)/(x^n - 1), lower=lower, upper=upper)
integrate(function(x) 2 * x^2/(x^(2*n) - 1), lower=sqrt(lower), upper=sqrt(upper))


integrate(function(x) x^(1/2)/(x^(n/3) - 1), lower=lower, upper=upper)
integrate(function(x) 6 * x^8/(x^(2*n) - 1), lower=lower^(1/6), upper=upper^(1/6))


####################

### Derived from
### Cardan Polynomials


lim = c(2.5, 3)
#
I.f(function(x) -(cos(x) + 2) / (8 - 6*cos(x) - 12*cos(x)^2 - 5*sin(2*x)), lim=pi/2 - 2*atan(lim))
I.f(function(x)  (sin(x) + 2) / (8 - 6*sin(x) - 12*sin(x)^2 - 5*sin(2*x)), lim=2*atan(lim))
# Derivation
I.f(function(x) (sin(2*x) + 2) / (4 - 3*sin(2*x) - 6*sin(2*x)^2 - 5/2*sin(4*x)), lim=atan(lim))
I.f(function(x) 1/4 * (sin(2*x) + 2) / (1 + sin(x)*cos(x) - 6*sin(x)^2*cos(x)^2 - 5*sin(x)*cos(x)^3), lim=atan(lim))
1/2 * I.f(function(x) (x^3 - 1)/(x^5 - 5*x^3 + 5*x - 1), lim=lim)
### Variants
I.f(function(x)  sin(x) / (8 - 6*sin(x) - 12*sin(x)^2 - 5*sin(2*x)), lim=2*atan(lim))
1/2 * I.f(function(x) (x^2 - x)/(x^5 - 5*x^3 + 5*x - 1), lim=lim)
### Decomposition
I.f(function(x)  1 / (8 - 6*sin(x) - 12*sin(x)^2 - 5*sin(2*x)), lim=2*atan(lim))
1/4 * I.f(function(x) (x^3 - x^2 + x - 1)/(x^5 - 5*x^3 + 5*x - 1), lim=lim)

###
I.f(function(x) -(cos(x) - 2) / (8 + 6*cos(x) - 12*cos(x)^2 + 5*sin(2*x)), lim=pi/2 - 2*atan(lim))
I.f(function(x)  (sin(x) - 2) / (8 + 6*sin(x) - 12*sin(x)^2 + 5*sin(2*x)), lim=2*atan(lim))
# Derivation
I.f(function(x) (sin(2*x) - 2) / (4 + 3*sin(2*x) - 6*sin(2*x)^2 + 5/2*sin(4*x)), lim=atan(lim))
I.f(function(x) 1/4 * (sin(2*x) - 2) / (1 - sin(x)*cos(x) - 6*sin(x)^2*cos(x)^2 + 5*sin(x)*cos(x)^3), lim=atan(lim))
1/2 * I.f(function(x) -(x^3 + 1)/(x^5 - 5*x^3 + 5*x + 1), lim=lim)
### Variants
I.f(function(x)  sin(x) / (8 + 6*sin(x) - 12*sin(x)^2 + 5*sin(2*x)), lim=2*atan(lim))
1/2 * I.f(function(x) (x^2 + x)/(x^5 - 5*x^3 + 5*x + 1), lim=lim)
### Decomposition
I.f(function(x)  1 / (8 + 6*sin(x) - 12*sin(x)^2 + 5*sin(2*x)), lim=2*atan(lim))
1/4 * I.f(function(x) (x^3 + x^2 + x + 1)/(x^5 - 5*x^3 + 5*x + 1), lim=lim)

