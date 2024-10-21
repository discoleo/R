########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Fractions: Other
###
### draft v.0.2d


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


### Helper Constants
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


### History:

# [refactor] Integrals of type log(...) * log(...)
# moved to new file: Integrals.Log.Prod.R;


###################
###################


### I( |log( (1-x) / (1+x) )|^n ) on [0, 1]
# Maths 505: This integral is actually one of your favorite constants
# https://www.youtube.com/watch?v=83mUOaF7G9A

n = sqrt(3)
integrate(\(x) abs(log((1-x)/(1+x)))^n, 0, 1)
integrate(\(x) log((1+x)/(1-x))^n, 0, 1)
2*gamma(n+1)*(1 - 2^(1-n)) * pracma::zeta(n)


### I( x * log((1+x)/(1-x))^2 / (1-x^2) )
integrate(\(x) x * log((1+x)/(1-x))^2 / (1-x^2) - 1/2 * log((1-x)/2)^2/(1-x), 0, 1)
- 1/3 * gamma(4)*(1 - 1/4) * pracma::zeta(3) + log(2)^3 / 6

# library(Rmpfr)
x = integrate(\(x) {
	x = mpfr(x, 240);
	y = x * log((1+x)/(1-x))^(n-1) / (1-x^2) - 1/2*log((1-x)/2)^2/(1-x);
	as.numeric(y); }, 0, 1, rel.tol=1E-8)
x$value + 1/3 * gamma(4)*(1 - 1/4) * pracma::zeta(3) - log(2)^3 / 6


### I( log((x - 1) / (x + 1)) / x ) on [1, Inf]
# Maths 505: A RIDICULOUSLY AWESOME LOG INTEGRAL!!!
# https://www.youtube.com/watch?v=tTu6hedSlm0

# numerical issue:
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, Inf)
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, 100)
# alternative:
integrate(\(x) log((x - 1) / (x + 1)) / x, 1, Inf)
integrate(\(x) log((1 - x) / (x + 1)) / x, 0, 1)
- pi^2/4


### Simple:
integrate(\(x) log((1 - x) / (x + 1)), 0, 1)
- 2*log(2)

###
integrate(\(x) log(x) * log((1 - x) / (x + 1)), 0, 1)
2*log(2) - pi^2/12


###################
###################

### Types: log(x +/- 1) / (x^2+1)

# I( x * log(1+x) / (x^2 + 1) ) on [0, 1]
# I( x * log(1-x) / (x^2 + 1) ) on [0, 1]
# Cipher: Integrate xln(1-x)/(1 + x^2)dx from 0 to 1
# https://www.youtube.com/watch?v=CMRvFM2N7sw

###
integrate(\(x) log(1+x) / (x^2+1), 0, 1)
pi * log(2)/8

###
integrate(\(x) log(1-x) / (x^2+1), 0, 1)
pi*log(2)/8 - Catalan

###
integrate(\(x) x * log(1+x) / (x^2+1), 0, 1)
pi^2 / 96 + log(2)^2 / 8

###
integrate(\(x) x * log(1-x) / (x^2+1), 0, 1)
- 5/96 * pi^2 + log(2)^2 / 8


### I( log(1+x) / (x^2 + x + 1) )
integrate(\(x) log(1+x) / (x^2+x+1), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36

### I( log(1-x) / (x^2 - x + 1) )
integrate(\(x) log(1-x) / (x^2-x+1), 0, 1)
- (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 6

### I( x * log(1+x) / (x^2 + x + 1) )
integrate(\(x) x * log(1+x) / (x^2+x+1), 0, 1)
# TODO

### I( x * log(1-x) / (x^2 - x + 1) )
integrate(\(x) x * log(1-x) / (x^2-x+1), 0, 1)
pi^2/12 - pracma::psi(1, 1/3) / 6


### I( log(1-x) / (x^3 + 1) )
integrate(\(x) log(1-x) / (x^3 + 1), 0, 1)
- pi^2/18 + log(2)^2/6 +
	- pracma::psi(1, 1/3) / 18 + pracma::psi(1, 2/3) / 9;


### Li2(1/2)
integrate(\(x) log(1+x) / (1-x) - log(2)/(1-x), 0, 1)
- pi^2/12 + log(2)^2/2
#
integrate(\(x) log(1-x) / (x+1), 0, 1)
- pi^2/12 + log(2)^2/2

#
integrate(\(x) log(1+x) * (x^2 + 2) / (1-x^3) - log(2)/(1-x), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 36 +
	- pi^2/12 + log(2)^2/2;


### Pow = 4

###
integrate(\(x) log(1-x) / (x^4+1), 0, 1)

###
integrate(\(x) x * log(1-x) / (x^4+1), 0, 1)

###
integrate(\(x) x^2 * log(1-x) / (x^4+1), 0, 1)

### Helper
integrate(\(x) x * log(1-x^4) / (x^4+1), 0, 1)
pi*log(2)/8 - Catalan / 2

#
integrate(\(x) x * log(1-x^2) / (x^4+1), 0, 1)
pi*log(2)/16 - Catalan/2


##################
##################

### I( log((b - x)/(b + x)) / (x * sqrt(1 - x^2)) ) on [0, 1]
# Maths 505: DESTROYING a MONSTER integral using Feynman's technique
# https://www.youtube.com/watch?v=PWtzeKvOVEo

b = sqrt(3)
integrate(\(x) log((b + x)/(b - x)) / (x * sqrt(1 - x^2)), 0, 1)
pi^2 / 2 - pi*acos(1/b)

### Derived:
b = sqrt(3)
integrate(\(x) 1 / ((b^2 - x^2) * sqrt(1 - x^2)), 0, 1)
pi / sqrt(b^2 - 1) / (2*b)


################
################

# TODO: move to Trig

### I( atan(x) * log((1-x)/(1+x)) )
# Maths 505: Another INSANE integral!
# https://www.youtube.com/watch?v=KEDEzVqlAYU

integrate(\(x) atan(x) * log((1-x)/(1+x)), 0, 1)
pi^2/16 - pi/4*log(2) - Catalan;


### I( atan(x) * log(1-x) )
integrate(\(x) atan(x) * log(1-x), 0, 1)
integrate(\(x) ((1-x)*log(1-x) + x) / (x^2+1) - pi/4, 0, 1)
5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan;

### I( atan(x) * log(1+x) )
integrate(\(x) atan(x) * log(1+x), 0, 1)
-1/96*pi^2 + 3/8*pi*log(2) - log(2)^2 / 8 + log(2)/2 - pi/4;


# Gen 1: TODO
k = 2
integrate(function(x) atan(x) * log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


###
integrate(function(x) atan(x) / x, 0, 1)
Catalan

# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


##################

### I( log(x) * atan(x) / x )
# 1.) Maths 505: a superb integral sprinkled with some fourier analysis
# https://www.youtube.com/watch?v=aHCLS4l65VE
# 2.a) see also: Sums.Fractions.Higher.R;
# 2.b) Flammable Maths: An AMAZING Journey of Series Evaluation!
# Calculating Euler's Sum! [ Series pi^3/32 (-1)^k/(2k+1)^3 ]

integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32


### I( log(x) * atan(x) )
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2/48 - pi/4 + log(2)/2

# Derivation:
id = seq(0, 20000)
- sum((-1)^id * (1/(2*id+1) - 1/(2*id+2))) + pi^2/48
- sum((-1)^id / (2*id+1)) + pi^2/48 + log(2)/2
pi^2/48 - pi/4 + log(2)/2


### I( log(x) * atan(x^p) / x )
integrate(\(x) log(x) * atan(x^2) / x, 0, 1)
- pi^3 / 128

### Gen:
p = sqrt(5)
integrate(\(x) log(x) * atan(x^p) / x, 0, 1)
- pi^3 / (32 * p^2)

### Other:
integrate(\(x) log(x) * x / (x^2 + 1), 0, 1)
- pi^2/48


### Pow = 2
p = sqrt(5)
integrate(\(x) log(x)^2 * x^(p-1) / (x^(2*p) + 1), 0, 1)
pi^3 / (16 * p^3)

###
integrate(\(x) log(x)^2 / (x^2 + 1), 0, 1)
pi^3 / 16


###
integrate(\(x) log(x^2 + 1) / x, 0, 1)
pi^2/24


######################
######################

### I( log(sqrt(x^2 + 1) + x) / (x^2 + 1) )
# Maths 505: An AWESOME log integral with a surprising result
# https://www.youtube.com/watch?v=a2On9rAvqks

integrate(\(x) log(sqrt(x^2 + 1) + x) / (x^2 + 1), 0, Inf)
integrate(\(x) log(sqrt(x^2 + 1) + 1) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - x) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - 1) / (x^2 + 1), 0, Inf)
2*Catalan


### Derived / Sub-Integrals:
integrate(\(x) log((1 + cos(x)) / cos(x)), 0, pi/2)
2*Catalan

###
# - see also file: Integrals.Log.Trig.R
integrate(\(x) log(1 + cos(x)), 0, pi/2)
integrate(\(x) 4*log(cos(x)), 0, pi/4)$value + pi*log(2)/2
- pi*log(2)/2 + 2*Catalan


#######################
#######################

### Varia:

###
n = sqrt(5)
integrate(function(x) log(2 - x^n) / (2 - x^n)^(1+1/n), lower=0, upper=1)
integrate(function(x) - n*(n+1) * x^n * log(x) / (2 - x^n)^(2+1/n), lower=0, upper=1)
# TODO


#######################
#######################

### I( x * log(x) / cosh(k*x)^2 )
# Maths 505: Feynman's trick solves an AWESOME integral
# https://www.youtube.com/watch?v=qZywLtQaDGQ

###
k = sqrt(5)
integrate(\(x) x * log(x) / cosh(k*x)^2, 0, 100)
- (log(2)*log(k) + 3/2*log(2)^2 - log(2)) / k^2

### Helper:
k = sqrt(3)
integrate(\(x) x * exp(k*x) / (exp(k*x) + 1)^2, 0, 100)
log(2) / k^2


#######################
#######################

### 1 / LOG()

###
integrate(\(x) 1/log(x) + 1/(1-x), 0, 1)
Euler


### I( (x^n - 1) / log(x) )
# Michael Penn: My favorite way "around" Feynman's trick -- and a challenge for you!
# https://www.youtube.com/watch?v=EseFhG1QARM

n = sqrt(5)
integrate(\(x) (x^n - 1) / log(x), 0, 1)
log(n+1)

###
n = -1/2
integrate(\(x) x^n / log(x) + 1/(1-x), 0, 1, rel.tol=1E-8)
log(n+1) + Euler

###
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1) / log(x), 0, 1)
log(n+p+1) - log(p+1)


### I( 1 / ((x+1) * log(x)) )
integrate(\(x) 1 / ((x+1) * log(x)) + 1/(1-x^2), 0, 1, rel.tol=1E-8)
- log(gamma(1/2)) - digamma(1/2)/2

### I( (x-1) / ((x+1) * log(x)) )
integrate(\(x) (x-1) / ((x+1) * log(x)), 0, 1, rel.tol=1E-8)
log(pi/2)


### I( x^p * (x^n - 1)^2 / log(x)^2 )
# Maths 505: How Richard Feynman would solve this awesome golden integral
# https://www.youtube.com/watch?v=g2NPdw4ig5M

n = sqrt(5)
integrate(\(x) (x^n - 1)^2 / log(x)^2, 0, 1)
(2*n+1)*log(2*n+1) - (2*n+2)*log(n+1)


### I( 1 / log(x)^2 )
integrate(\(x) 1 / log(x)^2 - 1/(1-x)^2 + 1/(1-x), 0, 1, rel.tol=1E-6)
Euler - 1/2


### I( x^p / log(x)^2 )
p = sqrt(3)
integrate(\(x) x^p / log(x)^2 - 1/log(x)^2 - p/log(x), 0, 1)
(p+1)*log(p+1) - p;

### I( x^p * (x^n - 1)^2 / log(x)^2 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^2 / log(x)^2, 0, 1)
(2*n+p+1)*log(2*n+p+1) - 2*(n+p+1)*log(n+p+1) + (p+1)*log(p+1)


### Pow = 3

### I( (x^n - 1)^3 / log(x)^3 )
n = sqrt(5)
integrate(\(x) (x^n - 1)^3 / log(x)^3, 0, 1)
1/2*(3*n+1)^2*log(3*n+1) - 3/2*(2*n+1)^2*log(2*n+1) + 3/2*(n+1)^2*log(n+1)

### I( x^p * (x^n - 1)^3 / log(x)^3 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^3 / log(x)^3, 0, 1)
1/2*(3*n+p+1)^2*log(3*n+p+1) - 3/2*(2*n+p+1)^2*log(2*n+p+1) +
	+ 3/2*(n+p+1)^2*log(n+p+1) - 1/2*(p+1)^2*log(p+1)

# Note:
# shortcut: integrate (f(3*n) - 3*f(2*n) + 3*f(n)) dp;


### I( x^p * (x^n - 1)^4 / log(x)^4 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^4 / log(x)^4, 0, 1)
1/6*(4*n+p+1)^3*log(4*n+p+1) - 2/3*(3*n+p+1)^3*log(3*n+p+1) +
	+ (2*n+p+1)^3*log(2*n+p+1) - 2/3*(n+p+1)^3*log(n+p+1) +
	+ 1/6*(p+1)^3*log(p+1);


### I( x^p * (x^n - 1)^5 / log(x)^5 )
n = sqrt(5); p = sqrt(3)
integrate(\(x) x^p * (x^n - 1)^5 / log(x)^5, 0, 1)
((5*n+p+1)^4*log(5*n+p+1) - 5*(4*n+p+1)^4*log(4*n+p+1) +
	+ 10*(3*n+p+1)^4*log(3*n+p+1) - 10*(2*n+p+1)^4*log(2*n+p+1) +
	+ 5*(n+p+1)^4*log(n+p+1) - (p+1)^4*log(p+1)) / 24;


##################
##################

### I( 1 / (1 - x*log(x)) ) on [0, 1]
# Flammable Maths: I bet BPRP can not solve this Integral
# https://www.youtube.com/watch?v=p1uJkif5zE0

###
id = seq(0, 30)
integrate(\(x) 1 / (1 - x*log(x)), 0, 1)
sum( (-1)^id * gamma(id + 1) / (id + 1)^(id + 1) )

###
integrate(\(x) 1 / (1 + x*log(x)), 0, 1)
sum( gamma(id + 1) / (id + 1)^(id + 1) )

# TODO: any closed formulas?


##################
##################

### I( 1 / ((x+1) * (|log(x)|^n + b)) )
# Maths 505: A fun little integral exploration
# https://www.youtube.com/watch?v=Hl_ko1RSO3I

b = sqrt(3)
n = 6
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf)
pi/n * b^(1/n-1) / sin(pi/n)

#
b = sqrt(5)
n = sqrt(7)
# accuracy issues with small n;
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf, rel.tol = 1E-8, subdivisions = 1025)
pi/n * b^(1/n-1) / sin(pi/n)

#
b = 5
n = sqrt(7)
integrate(\(x) 1 / ((x+1) * (abs(log(x))^n + b)), 0, Inf, rel.tol = 1E-8, subdivisions = 1025)
pi/n * b^(1/n-1) / sin(pi/n)

