########################
###
### Leonard Mada
### [the one and only]
###
### Integral Tricks
###
### draft v.0.3r-var2

# TODO: refactor & split file;


### various Integral Tricks

### draft v.0.3d:
# - moved section on Moebius Transform
#   to separate file:
#   Integrals.Moebius.Transform.R;
### draft v.0.3c - v.0.3c-ex2:
# - definite integrals using Moebius Transform:
#   I( log(4*x-3) / (x*(3*x+1)) ) on c(1, 4);
#   I( log(5*x-3) / (x*(x+1)) ) on c(1, 3);
#   I( log(x-3) / (x*(3*x+1)) ) on c(5, 8);
### draft v.0.3a - v.0.3b:
# - Transformations:
#   f((a*x + b) / (c*x + d)) = - f(x);
# - added 2nd example;
# - acos() variant; [v.0.3b]
### draft v.0.2f:
# - started work on:
#   I( log(x-1) / x ) dx;
#   [see v.0.3a-d for related results]
### draft v.0.2d - v.0.2e:
# - more exponential tricks:
#   I( atan(exp(x)) ) dx;
#   I( atan(exp(x)) / (x^2 + 1) ) dx;
### draft v.0.2c:
# - partial results:
#   I( atan(x)^3 / x ) dx;
### draft v.0.2b:
# - just for completion:
#   I( atan(x) / x ) dx;
### draft v.0.2a:
# - sin(log(x)) type:
#   I( sin(log(x)) / (x^2+b*x+1) ) dx;
### draft v.0.1i - v.0.1k:
# - atan() trick:
#   I( atan(x) / (x^2+b*x+1) ) dx;
#   I( atan(sqrt(x)) / (x^2+b*x+1) ) dx;
#   I( atan(x) / sqrt(x^4+b*x^3+b*x+1) ) dx;
### draft v.0.1f - v.0.1h:
# - experiments with complex integrals;
#   I( log(x) / (x^2 - 1) ) dx over or using complex path;
### draft v.0.1e:
# - extension of the log() trick;


##################
##################

##################
### Logarithms ###
##################

### Feynman trick & Other tricks
# - usually (much) simpler than Contour integration;

### Basic / Helper
b = sqrt(3)
integrate(function(x) log(x) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b)/(2*b)
#
integrate(function(x) log(x) / (b^2*x^2 + 1), lower=0, upper=Inf)
- pi*log(b)/(2*b)

### TODO:
# - using Feynman's trick;
# - Note: Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)
Catalan = 0.915965594177219015054603514;
integrate(function(x) log(x + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)/4 + Catalan;

### TODO:
Catalan = 0.915965594177219015054603514;
integrate(function(x) log(x + 1) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b^2 + 1)/(4*b) + log(b)*atan(1/b)/b +
	- integrate(function(x) log(x) / (x^2 + b^2), 0, 1)$value;
pi*log(b^2 + 1)/(4*b) +
	- integrate(function(x) log(x) / (x^2 + 1), 0, 1/b)$value / b;

### TODO:
a = sqrt(5); b = sqrt(3)
integrate(function(x) log(x + a) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(a^2 + b^2)/(4*b) + log(b)*atan(a/b)/b +
	- integrate(function(x) log(x) / (x^2 + b^2), 0, a)$value;

#############
### log(P[2])
# - using Feynman's trick;
b = sqrt(3); a = sqrt(2)
integrate(function(x) log(x^2 + a^2) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(a + b)/b

### Helper
integrate(function(x) log(x^2 + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)

###
b = sqrt(3)
integrate(function(x) log(x^2 + 1) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b+1)/b

### [Special Case]
b = sqrt(3)
integrate(function(x) log(x^2 + b^2) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(2*b)/b

# "cyclic redundancy" / Extra Cases
b = sqrt(3)
integrate(function(x) log(x^2 + 1) / (x^2 + b^2), lower=0, upper=Inf)
integrate(function(x) log(x^2 + 1) / (1 + b^2*x^2), lower=0, upper=Inf)$value +
	+ pi*log(b)/b;
integrate(function(x) log(x^2 + b^2) / (x^2 + 1), lower=0, upper=Inf)$value / b;
pi*log(b+1)/b


### Generalization:
# - using Feynman's trick;
# - using formulas for Fractions with roots of unity, see:
#   Integrals.Fractions.Unity.R;
n = 3
a = sqrt(3); b = sqrt(5)
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)

###############
### Case: n = 3
n = 3
a = sqrt(3); b = sqrt(5)
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)
sqrt(3)*pi/3 * log((a^n - b^n)/(a - b)) / b^2 - pi^2/3/b^2 +
	+ 2*pi/3 * atan((2*a/b + 1)/sqrt(3)) / b^2;


# Helper:
integrate(function(x) log(x) / (x^3 + 1), lower=0, upper=Inf)
- 2*pi^2/27

#
b = sqrt(5)
integrate(function(x) log(x) / (x^3 + b^3), lower=0, upper=Inf)
- 2*pi^2/27/b^2 + 2*pi/(3*sqrt(3))*log(b)/b^2

# Derivation:
- 2*pi^2/27/b^2 + log(b)/b^2 * integrate(function(x) 1 / (x^3 + 1), lower=0, upper=Inf)$value
# - 2*pi^2/27/b^2 + log(b)/b^2/3 * (log(x+1) - 1/2*log(x^2-x+1) + 3/sqrt(3)*atan((2*x-1)/sqrt(3)))
- 2*pi^2/27/b^2 + log(b)/b^2*(pi/2 - atan(-1/sqrt(3)))/sqrt(3)

# dI: evaluated at Inf & at 0;
n = 3
integrate(function(x) n*a^(n-1)/(b^n-a^n)* (1/(x^n + a^n) -  1/(x^n + b^n)), lower=0, upper=Inf)
# [not run]
1/2/(b^n - a^n) * log((x+a)^2/(x^2-a*x+a^2)) +
	- 1/2*a^2/b^2/(b^n - a^n)*log((x+b)^2/(x^2-b*x+b^2)) +
	+ 3/sqrt(3)/(b^n - a^n)*(atan(2*x/(a*sqrt(3)) - 1/sqrt(3)) - a^2/b^2*atan(2*x/(b*sqrt(3)) - 1/sqrt(3)))
3/sqrt(3)/(b^n - a^n)*(atan(2*x/(a*sqrt(3)) - 1/sqrt(3)) - a^2/b^2*atan(2*x/(b*sqrt(3)) - 1/sqrt(3)))
# =>
sqrt(3)*(1 - a^2/b^2)/(b^n - a^n)*(pi/2 + atan(1/sqrt(3)))
2*sqrt(3)*pi/3 * (1 - a^2/b^2)/(b^n - a^n)

# back-Integration:
2*sqrt(3)*pi/3 * integrate(function(x) (x^2/b^2 - 1)/(x^n - b^n), lower=0, upper=a)$value +
	- 2*pi^2/9/b^2 + 2*pi/sqrt(3)*log(b)/b^2;
2*sqrt(3)*pi/9 * log((a^n - b^n)/(a - b)) / b^2 +
	- 2*pi^2/9/b^2 + 2*sqrt(3)*pi/9*log(b) / b^2 +
	+ 2*sqrt(3)*pi/9/b^2 * integrate(function(x) (x + 2*b)/(x^2+b*x+b^2), lower=0, upper=a)$value
sqrt(3)*pi/3 * log((a^n - b^n)/(a - b)) / b^2 - 2*pi^2/9/b^2 +
	+ 2*pi/3 *(atan((2*a/b + 1)/sqrt(3)) - atan(1/sqrt(3))) / b^2;


### Case: n = 5
n = 5
# dI:
integrate(function(x) n*a^(n-1)/(b^n-a^n)* (1/(x^n + a^n) -  1/(x^n + b^n)), lower=0, upper=Inf)
# [more complicated] x^2 - (m+m^4)*a*x + a^2 =>
# TODO:
sqrt(n)*(1 - (a/b)^(n-1))/(b^n - a^n) * (pi/2 + ...)

# Helper
integrate(function(x) log(x)/(x^5 + 1), 0, Inf)
- pi^2*cos(pi/5)/sin(pi/5)^2 / 25

#
b = 3^(1/4);
integrate(function(x) log(x)/(x^5 + b^5), 0, Inf)
- pi^2*cos(pi/5)/sin(pi/5)^2 / 5^2 / b^4 +
	+ pi*log(b)/sin(pi/5) / 5 / b^4;


# Derivation:
m = cos(pi/5) + 1i*sin(pi/5)
v = pi/sin(pi/5)/5
integrate(function(x) log(x)/(x^5 + 1), 0, Inf)
2i*pi*(log(m)/((m+1)*prod(m - m^c(3,7,9))) + (1/5)*v*exp(2i*pi/5)) / (1 - exp(2i*pi/5))
2i/25*pi^2*(1i/m^4 + exp(2i*pi/5)/sin(pi/5)) / (1 - exp(2i*pi/5))
2i/25*pi^2*(1i*exp(-4i*pi/5) + exp(2i*pi/5)/sin(pi/5)) / (1 - exp(2i*pi/5))
- pi^2*(exp(1i*pi/5)/sin(pi/5) - 1i) / sin(pi/5) / 25



### Generalization
n = 7
integrate(function(x) log(x)/(x^n + 1), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2


#
n = 7; b = 3^(1/4);
integrate(function(x) log(x)/(x^n + b^n), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2 / b^(n-1) +
	+ pi*log(b)/sin(pi/n) / n / b^(n-1);


##################
### Transforms ###
##################

###
Catalan = 0.915965594177219015054603514;

###
integrate(function(x) log(x+1) / (x^2 + 1), lower=-1, upper=1)
pi*log(2)/4 - Catalan;

# Base:
integrate(function(x) log(x + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)/4 + Catalan;
# x => (1-x)/(1+x)
integrate(function(x) (log(2) - log(x+1)) / (x^2 + 1), lower=-1, upper=1)
# =>
pi*log(2)/2 - integrate(function(x) log(x+1) / (x^2 + 1), lower=-1, upper=1)$value


### x => (a-x)/(a+x)
a = 7/5
integrate(function(x) log(x+a) / (x^2 + a^2), lower=-a, upper=a)
pi*log(a)/(2*a) + pi*log(2)/(4*a) - Catalan/a;


### x => (1-x)/(a+x)
a = 7/5
integrate(function(x) log(x+a) / (2*x^2 + 2*(a-1)*x + a^2 + 1), lower=-a, upper=1)
pi*log(a+1)/(2*(a+1)) - pi*log(2)/(4*(a+1)) - Catalan/(a+1);


##################
##################

### Example 1:
# - based on: "This trick is new to me!"
#   https://www.youtube.com/watch?v=JCuplIQ6JG4


### I (log(x) / (x^2 + b*x + 1)) dx
# lower = 1/a, upper = a, I = 0;
b = 3
lim = 2
integrate(function(x) log(x) / (x^2 + b*x + 1), lower=1/lim, upper=lim)

b = c(5, 2)
lim = 2
integrate(function(x) log(x) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)
integrate(function(x) log(b[1]) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)


### Higher Powers
a  = 3
b2 = 5
#
b  = c(1, b2, b2*a, a^3)
bf = c(1, b2/a, b2/a, 1)
# Note: polynomial fraction can be simplified;
integrate(function(x) sapply(x,
	function(x) (x+a)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)

integrate(function(x) sapply(x,
	function(x) log(a)/a * (x+1)/sum(bf*x^(3:0))), lower=0, upper=Inf)

###
integrate(function(x)
	sapply(x, function(x) sqrt(x)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)
integrate(function(x)
	sapply(x, function(x) log(a)*sqrt(a)/a^2 *sqrt(x)/ sum(bf*x^(3:0))), lower=0, upper=Inf)


#########
### Ex 2:
a = 3
b1n = 4
b1 = 2
# b2 can be added as well;

integrate(function(x) (x^2 + b1n*x + a^2)*log(x)/
	(x^4 + b1*x^3 + b1*a^2*x + a^4), lower=0, upper=Inf)
integrate(function(x) log(a)/a * (x^2 + b1n/a*x + 1)/
	(x^4 + b1/a*x^3 + b1/a*x + 1), lower=0, upper=Inf)

### only x:
integrate(function(x) x*log(x)/
	(x^4 + b1*x^3 + b1*a^2*x + a^4), lower=0, upper=Inf)
integrate(function(x) log(a)/a^2 * x /
	(x^4 + b1/a*x^3 + b1/a*x + 1), lower=0, upper=Inf)


#########
### Ex 3:
integrate(function(x) (x^3+8)*log(x)/(x^5 + 3*x^4 + 24*x + 32), lower=0, upper=Inf)
integrate(function(x) log(2)/2 * (x^3+1)/(x^5 + 3/2*x^4 + 3/2*x + 1), lower=0, upper=Inf)


###############

library(pracma)

# Finite Integrals

### I( log(x) / (x^2 - 1) ) dx
line_integral(function(x)  log(x) / (x^2+1), c((1-1i)/sqrt(2), (1+1i)/sqrt(2)))
line_integral(function(x) -log(-x*1i)*1i / (x^2-1), c((-1+1i)/sqrt(2), (1+1i)/sqrt(2)))
line_integral(function(x) (-1i*log(x) - pi/2) / (x^2-1), c((-1+1i)/sqrt(2), (1+1i)/sqrt(2)))


a = 2
line_integral(function(x)  log(x) / (x^2+1), c(1/a*(1-1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) -log(-x*1i)*1i / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) (-1i*log(x) - pi/2) / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))

a = 2
line_integral(function(x) log(x) / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) 1i*pi/2 / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
1i*pi/4 * log((a*(1+1i) - sqrt(2))*(1/a*(-1+1i) + sqrt(2)) /
	((1/a*(-1+1i) - sqrt(2))*(a*(1+1i) + sqrt(2))))

lim = c(1.1, 3) # lim[1] > 1; # strictly greater than 1!
line_integral(function(x) log(x) / (x^2-1), c(lim[1]*1i, lim[2]*1i))
1i*pi/4 * log((lim[2]*1i - 1)*(lim[1]*1i + 1) /
	((lim[2]*1i + 1)*(lim[1]*1i - 1))) +
	+ 1i * integrate(function(x) log(x) / (x^2+1), lower=1/lim[2], upper=lim[1])$value


lim = c(1.1, 3) # lim[1] > 1; # strictly greater than 1!
integrate(function(x) log(x) / (x^2-1), lower=lim[1], upper=lim[2])
1i * line_integral(function(x) log(-1i*x) / (x^2+1), c(lim[1]*1i, lim[2]*1i))
1i * line_integral(function(x) log(x) / (x^2+1), c(lim[1]*1i, lim[2]*1i)) +
	+ pi/2 * (atan(lim[2]*1i) - atan(lim[1]*1i))


### TODO:
# - find full path;


########################

# Michael Penn: I really like this integral
# https://www.youtube.com/watch?v=5sxtrPgHFw8

###
n = 4
a = sqrt(2)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


###
n = sqrt(5)
a = sqrt(2)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


###
n = 3
a = 2^(1/3)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


########################

# qncubed3: Complex Analysis: An Integral from @MichaelPennMath
# https://www.youtube.com/watch?v=LH4i9XJsz_I
# - Generalization of the 2nd (sub-) Integral;

p = 1/5
k = 3
integrate(function(x) log(x)/x^p/(x+k)^2, lower=0, upper=Inf)
pi*((p*log(k) - 1)*sin(pi*p) + pi*p*cos(pi*p)) / k^(p+1) / sin(pi*p)^2


# Residue: x = - k;
# x^(p-1)*(1 - p*log(x))
2i*pi*k^(-p-1)*exp(-1i*pi*(p+1))*(1 - p*log(k) - 1i*pi*p)

# Div:
(1 - exp(-2i*pi*p))

# Helper
integrate(function(x) x^(-p)/(x+k)^2, lower=0, upper=Inf)
pi*p/sin(pi*p)/k^(p+1)


########################

# qncubed3: Complex Analysis: Is this solution way too overkill?
# https://www.youtube.com/watch?v=vfa5wGE7fRI

###
integrate(function(x) log(x)^2/(x^2 + 1), lower=0, upper=Inf)
pi^3/8

### Gen 1:
b = sqrt(5)
integrate(function(x) log(x)^2/(x^2 + b^2), lower=0, upper=Inf)
pi^3/(2^3*b) + pi/2 * log(b)^2/b


### Gen: Power of log

###
b = sqrt(5)
integrate(function(x) log(x)^3/(x^2 + b^2), lower=0, upper=Inf)
pi*log(b)^3 / (2*b) + 3*pi^3*log(b)/(8*b)


###
integrate(function(x) log(x)^4/(x^2 + 1), lower=0, upper=Inf)
5*pi^5/32

#
b = sqrt(5)
integrate(function(x) log(x)^4/(x^2 + b^2), lower=0, upper=Inf)
5*pi^5/(2^5*b) + pi/2 * log(b)^4/b + 6*pi^3*log(b)^2/(8*b)


###
integrate(function(x) log(x)^6/(x^2 + 1), lower=0, upper=Inf)
61*pi^7/128

#
b = sqrt(5)
integrate(function(x) log(x)^6/(x^2 + b^2), lower=0, upper=Inf)
61*pi^7/(2^7*b) + pi/2 * log(b)^6/b + 75*pi^5*log(b)^2/(2^5*b) + 15*pi^3*log(b)^4/(8*b)


###
integrate(function(x) log(x)^8/(x^2 + 1), lower=0, upper=Inf)
1385*pi^9/512

#
b = sqrt(5)
integrate(function(x) log(x)^8/(x^2 + b^2), lower=0, upper=Inf)
1385*pi^9/(2^9*b) + pi/2 * log(b)^8/b + 28*61*pi^7*log(b)^2/(2^7*b) +
	+ 70*5*pi^5*log(b)^4/(2^5*b) + 28*pi^3*log(b)^6/(8*b)


###
integrate(function(x) log(x)^10/(x^2 + 1), lower=0, upper=Inf,
	abs.tol=1E-8, rel.tol=1E-8, subdivisions=4096)
50521*pi^11/2^11


### Gen 2: Both Powers

###
integrate(function(x) log(x)^3/(x^3 + 1), lower=0, upper=Inf)
pi^4*(25/3^4 - 1) / 12


# Aux. Polynomial:
x^4 - 4i*pi*x^3 - 4*pi^2*x^2

# Derivation:
- pi^4*(1 +
	- (1/3^4 - 4/3^3 + 4/9)*m +
	- (5^4/3^4 - 4*5^3/3^3 + 4*5^2/9)*m^5) / 12
- pi^4*(1 - 25*m/3^4 - 25*m^5/3^4) / 12
- pi^4*(1 - 25/3^4) / 12


###
integrate(function(x) log(x)^5/(x^5 + 1), lower=0, upper=Inf)
pi^6*(2*4779*cos(pi/5) - 4*31311*cos(pi/5)^2 + 2*31311 - 3*5^6) / (6*5*5^6)

# Aux. Polynomial:
x^6 - 6i*pi*x^5 - 10*pi^2*x^4 - 8*pi^4*x^2
b5 = - 3*k; b4 = 5/2*k^2; b3 = 0; b2 = -1/2*k^4;

# Derivation:
source("Polynomials.Helper.R")

n = 6
px = toPoly.pm(paste0("x^", seq(n-1), "*b", seq(n-1), collapse="+"));
px = sum.pm(px, toPoly.pm("x^n"));
pd = diff.pm(replace.pm(px, toPoly.pm("x+k"), "x"), px);
toCoeff(pd, "x")

px = toPoly.pm("x^6 - 6i*pi*x^5 - 10*pi^2*x^4 - 8*pi^4*x^2")
eval.pm(px, list(x=1i*pi/5, pi=pi)) / pi^6 * 5^6
# =>
- pi^6*(3*5^6/5 + 4779/(5*m^4) + 31311/(5*m^12) + 31311/(5*m^28) + 4779/(5*m^36)) / (6*5^6)
- pi^6*(3*5^6 - 4779*m - 31311*m^3 + 31311*m^2 + 4779*m^4) / (6*5*5^6)
- pi^6*(3*5^6 - 2*4779*cos(pi/5) + 2*31311*cos(2*pi/5)) / (6*5*5^6)
pi^6*(2*4779*cos(pi/5) - 4*31311*cos(pi/5)^2 + 2*31311 - 3*5^6) / (6*5*5^6)



########################
########################

#####################
### Trigonometric ###
#####################

### sin(log(x))
a = 3
b = 3
integrate(function(x) sin(log(x)) / (x^2+b*x+1), lower=1/a, upper=a)
# 0
integrate(function(x) cos(log(exp(pi/2)/x)) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) cos(log(exp(pi/2)*x)) / (x^2+b*x+1), lower=1/a, upper=a)
# 0

integrate(function(x) sin(log(x)) / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)
# 0


########################
########################

### ARCTAN

a = 3
integrate(function(x) atan(x) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 / x, lower=1/a, upper=a)
pi/2*log(a)


a = 3
b = 3
integrate(function(x) atan(x) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / (x^2+b*x+1), lower=1/a, upper=a)


a = 3
b = 3
integrate(function(x) x*atan(x) / (x^4+b*x^3+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 * x / (x^4+b*x^3+b*x+1), lower=1/a, upper=a)


### Sqrt
a = 3
integrate(function(x) atan(sqrt(x)) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 / x, lower=1/a, upper=a)

a = 3
b = 3
integrate(function(x) atan(sqrt(x)) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / (x^2+b*x+1), lower=1/a, upper=a)


a = 3
integrate(function(x) atan(x) * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(function(x) pi/2 * sqrt(x^4+1) / x^2, lower=1/sqrt(a), upper=sqrt(a))


a = 3
b = 3
integrate(function(x) atan(x) / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)


### Power 3
a = 3
integrate(function(x) atan(x)^3 / x, lower=1/a, upper=a)
integrate(function(x) (3*pi/4*atan(x)^2 - pi^3/32) / x, lower=1/a, upper=a)
integrate(function(x) 3*pi/4*atan(x)^2 / x, lower=1/a, upper=a)$value - pi^3/16*log(a)


### Exp
a = 3
integrate(function(x) (1 - exp(-pi/2))*exp(atan(x)) / x, lower=1/a, upper=a)
integrate(function(x) (exp(pi/2) - 1)*exp(-atan(x)) / x, lower=1/a, upper=a)

integrate(function(x) (exp(atan(x)) - exp(pi/2)*exp(-atan(x))) / x, lower=1/a, upper=a)
# 0

a = 3
integrate(function(x) exp(atan(x)) / (x * (exp(atan(x)) + exp(-atan(x)))), lower=1/a, upper=a)
# TODO


### atan(exp(x))
a = 3
integrate(function(x) atan(exp(x)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)), lower=-a, upper=a)
a * pi/2

integrate(function(x) atan(exp(x^3)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x^3)), lower=-a, upper=a)
a * pi/2


a = 3
integrate(function(x) atan(exp(x)) / (x^2 + 1), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)) / (x^2 + 1), lower=-a, upper=a)
pi/2 * atan(a)

integrate(function(x) atan(exp(x)) * x^2 / (x^2 + 1), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)) * x^2 / (x^2 + 1), lower=-a, upper=a)
pi/2 * (a - atan(a))


####################
####################

a = 2
integrate(function(x) log(x-1) / x, lower=a, upper=a+1)
integrate(function(x) log(x-1) / (x*(x-1)), lower=a/(a-1), upper=(a+1)/a)

i.f = function(x) 1/2*log(x-1)^2
i.f((a+1)/a) - i.f(a/(a-1))
integrate(function(x) log(x-1) / x, lower=a, upper=a+1)$value +
	+ integrate(function(x) log(x-1) / x, lower=a/(a-1), upper=(a+1)/a)$value
1/2*log(a)^2 - 1/2*log(a-1)^2

# TODO: both limits;


########################
########################
