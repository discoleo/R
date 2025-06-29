

### Trig Fractions

# I( sin(k1*x) / (exp(k2*x) - 1) )
# I( sin(k1*x) / (exp(k2*x) + 1) )
# I( cos(k1*x) / (exp(k2*x) - 1) )
# I( cos(k1*x) / (exp(k2*x) + 1) )
# I( sin(k1*x)^2 / (exp(k2*x) - 1) )


##############

### History

# - [refactor] Basic functions moved
#   to new file: Integrals.Exp.Trig.Basic.R;


####################

### Helper Constants

Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;


#####################
#####################

### I( sin(x) / (exp(k*x) - 1) )
# Maths 505: One of the coolest integrals on YouTube!
# https://www.youtube.com/watch?v=I7a-6VwKhgo

### I( sin(x) / (exp(x) - 1))
integrate(\(x) sin(x) / (exp(x) - 1), 0, Inf)
- (pracma::psi(0, 1+1i) - pracma::psi(0, 1-1i)) * 1i / 2
Im(pracma::psi(0, 1+1i))
#
pi/2 / tanh(pi) - 1/2
(digamma((1/2+1)/2) - digamma(1/4))/2 - 1/2
- (pracma::psi(0, 1i) - pracma::psi(0, - 1i)) * 1i / 2 - 1


### I( cos(x) / (exp(x) - 1) )
integrate(\(x) cos(x) / (exp(x) - 1) - exp(-x)/x, 0, Inf, rel.tol = 1E-9)
- Re(pracma::psi(1i))


### I( sin(k*x) / (exp(x) - 1) ) on [0, Inf]
k = sqrt(3)
integrate(\(x) sin(k*x) / (exp(x) - 1), 0, Inf)
- (pracma::psi(0, 1i*k + 1) - pracma::psi(0, - 1i*k + 1)) * 1i / 2;
- (pracma::psi(0, 1i*k) - pracma::psi(0, - 1i*k)) * 1i / 2 - 1/k;
Im(pracma::psi(0, 1 + 1i*k));


### I( sin(k1*x) / (exp(k2*x) - 1) ) on [0, Inf]
k1 = sqrt(3); k2 = sqrt(5);
integrate(\(x) sin(k1*x) / (exp(k2*x) - 1), 0, Inf)
Im(pracma::psi(0, 1 + 1i*k1/k2)) / k2;


### I( cos(k*x) / (exp(x) - 1) ) on [0, Inf]
k = sqrt(3)
integrate(\(x) cos(k*x) / (exp(x) - 1) - exp(-x)/x, 0, Inf)
- Re(pracma::psi(0, 1 + 1i*k));


### I( cos(k1*x) / (exp(k2*x) - 1) ) on [0, Inf]
k1 = sqrt(3); k2 = sqrt(5);
integrate(\(x) cos(k1*x) / (exp(k2*x) - 1) - 1/k2 * exp(-k2*x)/x, 0, Inf)
- Re(pracma::psi(0, 1 + 1i*k1/k2)) / k2;
#
integrate(\(x) cos(k1*x) / (exp(k2*x) - 1) - 1/k2 * exp(-x)/x, 0, Inf)
- Re(pracma::psi(0, 1 + 1i*k1/k2)) / k2 - log(k2)/k2;

# Note:
k = sqrt(7);
integrate(\(x) (exp(-x) - exp(-k*x)) / x, 0, Inf)
log(k);


### Pow: I( sin(k*x)^2 / (exp(x) - 1) )
k = sqrt(3)
integrate(\(x) sin(k*x)^2 / (exp(x) - 1), 0, Inf)
integrate(\(x) (1 - cos(k*x)^2) / (exp(x) - 1), 0, Inf)
(Re(pracma::psi(0, 1 + 2i*k)) + Euler) / 2;


### Gen:

### I( sin(x) / (exp(k*x) - 1) ) on [0, Inf]
integrate(\(x) sin(x) / (exp(k*x) - 1), 0, Inf)
- (pracma::psi(0, 1i/k) - pracma::psi(0, - 1i/k)) * 1i / (2*k) - 1;

### I( sin(x) / (exp(k*x) + 1) ) on [0, Inf]
k = sqrt(3)
integrate(\(x) sin(x) / (exp(k*x) + 1), 0, Inf)
(pracma::psi(0, 1i/(2*k)) - pracma::psi(0, - 1i/(2*k))) * 1i / (2*k) +
	- (pracma::psi(0, 1i/k) - pracma::psi(0, - 1i/k)) * 1i / (2*k) + 1;


### I( sin(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(2); k = sqrt(3)
integrate(\(x) sin(q*x) / (exp(k*x) + 1), 0, Inf)
(pracma::psi(0, q*1i/(2*k)) - pracma::psi(0, - q*1i/(2*k)) +
	- pracma::psi(0, q*1i/k) + pracma::psi(0, - q*1i/k)) * 1i / (2*k) + 1/q;

### I( cos(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(2); k = sqrt(3)
integrate(\(x) cos(q*x) / (exp(k*x) + 1), 0, Inf)
- (pracma::psi(0, (q*1i/k + 1)/2) + pracma::psi(0, (- q*1i/k + 1)/2) +
	- pracma::psi(0, q*1i/k/2) - pracma::psi(0, - q*1i/k/2)) / (4*k)

# Derivation
integrate(\(x) Re(x^(q*1i)) / x / (x^k + 1), 1, Inf)
integrate(\(x) Re(x^(q*1i - 1)) - Re(x^(q*1i - 1)) / (x^k + 1), 0, 1)
- Re(pracma::psi(0, (q*1i/k + 1)/2) - pracma::psi(0, q*1i/k/2)) / (2*k);
# direct:
- (pracma::psi(0, (q*1i/k + 1)/2) + pracma::psi(0, (- q*1i/k + 1)/2) +
	- pracma::psi(0, q*1i/k/2) - pracma::psi(0, - q*1i/k/2)) / (4*k)


### I( x * Trig(q*x) / ... )

### I( x * sin(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(2); k = sqrt(3)
integrate(\(x) x * sin(q*x) / (exp(k*x) + 1), 0, Inf)
(pracma::psi(1, (q*1i/k + 1)/2) - pracma::psi(1, (- q*1i/k + 1)/2) +
	- pracma::psi(1, q*1i/k/2) + pracma::psi(1, - q*1i/k/2)) * 1i / (8*k^2)

### I( x * cos(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(2); k = sqrt(3)
integrate(\(x) x * cos(q*x) / (exp(k*x) + 1), 0, Inf)
- (pracma::psi(1, q*1i/(2*k)) + pracma::psi(1, - q*1i/(2*k)) +
	- 2*pracma::psi(1, q*1i/k) - 2*pracma::psi(1, - q*1i/k)) / (4*k^2) - 1/q^2;


### Variant: x^2
# see Sections below for the Generalization x^p (p = integer);

### I( x^2 * sin(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(3); k = sqrt(2)
integrate(\(x) x^2 * sin(q*x) / (exp(k*x) + 1), 0, Inf)
(pracma::psi(2, q*1i/(2*k)) - pracma::psi(2, - q*1i/(2*k)) +
	- (pracma::psi(2, q*1i/k) - pracma::psi(2, - q*1i/k)) * 2^p ) * 1i / (2*k)^3 +
	- 2/q^3;

### I( x^2 * cos(q*x) / (exp(k*x) + 1) ) on [0, Inf]
q = sqrt(2); k = sqrt(3)
integrate(\(x) x^2 * cos(q*x) / (exp(k*x) + 1), 0, Inf)
- (pracma::psi(2, (q*1i/k + 1)/2) + pracma::psi(2, (- q*1i/k + 1)/2) +
	- pracma::psi(2, q*1i/k/2) - pracma::psi(2, - q*1i/k/2)) / (16*k^3);


##########
### Other:
integrate(\(x) sin(x) / (exp(x) * (exp(x) - 1)), 0, Inf)
pi/2 / tanh(pi) - 1
#
integrate(\(x) cos(x) / (exp(x) * (exp(x) - 1)) - exp(-x) / x, 0, Inf)
- Re(pracma::psi(1i)) - Re(1/(1i+1))

### x^p
integrate(\(x) x * sin(x) / (exp(x) - 1), 0, Inf)
id = seq(10000)
2 * sum(id/(id^2 + 1)^2)
# TODO: see file: Sums.Fractions.Higher.R;
(pracma::psi(1, 1i) - pracma::psi(1, - 1i)) * 1i / 2


### I( x * sin(k*x) / (exp(x) - 1) )
integrate(\(x) x * sin(2*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, 2i) - pracma::psi(1, - 2i)) * 1i / 2
#
k = sqrt(2)
integrate(\(x) x * sin(k*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i / 2


### Generalization:
n = sqrt(3); # n >= 0 !!!
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
# [TODO] formula for sum();
# [done using polygamma function]
id = seq(10000)
gamma(n+1) * sum(sin((n+1)*pi/2 - (n+1)*atan(id)) / (id^2 + 1)^((n + 1)/2))

# Note:
# pracma::psi(n, ...): n must be integer;

### n = ODD
n = 3;
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
(pracma::psi(n, 1i) - pracma::psi(n, - 1i)) * 1i / 2
#
n = 5; k = sqrt(3)
integrate(\(x) x^n * sin(k*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(n, 1i*k) - pracma::psi(n, - 1i*k)) * 1i / 2
#
n = 5; q = sqrt(5); k = sqrt(3)
integrate(\(x) x^n * sin(q*x) / (exp(k*x) - 1), 0, Inf)
(pracma::psi(n, 1i*q/k) - pracma::psi(n, - 1i*q/k)) * 1i / (2*k^(n+1))


### n = EVEN
n = 2;
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
- (pracma::psi(n, 1i + 1) - pracma::psi(n, - 1i + 1)) * 1i / 2
#
n = 4; k = sqrt(3)
integrate(\(x) x^n * sin(k*x) / (exp(x) - 1), 0, Inf)
- (pracma::psi(n, 1i*k + 1) - pracma::psi(n, - 1i*k + 1)) * 1i / 2
#
n = 4; q = sqrt(5); k = sqrt(3)
integrate(\(x) x^n * sin(q*x) / (exp(k*x) - 1), 0, Inf)
- (pracma::psi(n, 1i*q/k + 1) - pracma::psi(n, - 1i*q/k + 1)) * 1i / (2*k^(n+1))


### Pow of sin

### n = ODD
n = 5; k = sqrt(3)
integrate(\(x) x^n * sin(k*x)^3 / (exp(x) - 1), 0, Inf)
(pracma::psi(n, 1i*k) - pracma::psi(n, - 1i*k)) * 3i / 8 +
	- (pracma::psi(n, 3i*k) - pracma::psi(n, - 3i*k)) * 1i / 8

###
n = 5; k = sqrt(3)
integrate(\(x) x^n * sin(x)^3 / (exp(k*x) - 1), 0, Inf)
(pracma::psi(n, 1i/k) - pracma::psi(n, - 1i/k)) * 3i / (8*k^(n+1)) +
	- (pracma::psi(n, 3i/k) - pracma::psi(n, - 3i/k)) * 1i / (8*k^(n+1))
#
integrate(\(x) x^n * sin(x)^3 / (exp(k*x) + 1), 0, Inf)
(pracma::psi(n, 1i/k) - pracma::psi(n, - 1i/k)) * 3i / (8*k^(n+1)) +
	- (pracma::psi(n, 3i/k) - pracma::psi(n, - 3i/k)) * 1i / (8*k^(n+1)) +
	- (pracma::psi(n, 1i/(2*k)) - pracma::psi(n, - 1i/(2*k))) * 3i / (4*(2*k)^(n+1)) +
	+ (pracma::psi(n, 3i/(2*k)) - pracma::psi(n, - 3i/(2*k))) * 1i / (4*(2*k)^(n+1))


###############

### Derivation:
# Decomposition into sum:
integrate(\(x) sin(x) / exp(x), 0, Inf)
1/2
#
integrate(\(x) x * sin(x) / exp(x), 0, Inf)
1/2
#
integrate(\(x) x^2 * sin(x) / exp(x), 0, Inf)
1/2
#
sapply(seq(9), \(p) round(integrate(\(x) x^p * sin(x) / exp(x), 0, Inf)$value, 3))
sapply(seq(9), \(p) Im(gamma(p + 1)/(1-1i)^(p + 1)))
#
sapply(seq(9), \(n) round(integrate(\(x) x * sin(x) / exp(n*x), 0, Inf)$value, 5))
sapply(seq(9), \(n) round(Im(gamma(2)/(n - 1i)^(1 + 1)), 5))
sapply(seq(9), \(n) round(2*n*gamma(2)/(n^2 + 1)^(1 + 1), 5))

###
k = sqrt(5)
integrate(\(x) x * sin(x) / exp(k*x), 0, Inf)
2*k*gamma(2)/(k^2 + 1)^(1 + 1)

###
k = sqrt(5); n = sqrt(3)
integrate(\(x) x^n * sin(x) / exp(k*x), 0, Inf)
Im(gamma(n+1)/(k - 1i)^(n + 1))
gamma(n+1) / (k^2 + 1)^((n + 1)/2) * sin((n+1)*pi/2 - (n+1)*atan(k))


#########################
#########################

### I( cos(m*x) / (cosh(k*x) + cos(phi)) ) on [0, Inf]
# Blagouchine IV. Rediscovery of Malmsten’s integrals, their evaluation
# by contour integration methods and some related results.
# Ramanujan J (2014) 35:21–110; DOI: I 10.1007/s11139-013-9528-5

# - see intermediate I for Exercise 2 (p 50 of article);
k = sqrt(5); m = sqrt(3); phi = 1/3;
integrate(\(x) cos(m*x) / (cosh(k*x) + cos(phi)), 0, Inf)
pi / (k*sin(phi)) * sinh(phi*m/k) / sinh(pi*m/k);

### Special Case: m = 0
k = sqrt(5); phi = 1/3;
integrate(\(x) 1 / (cosh(k*x) + cos(phi)), 0, Inf)
phi / (k*sin(phi))


############

# for Base-Integrals: see Integrals.Exp.R;

### I( x^p * (cos(x) - exp(-k*x)) / (cosh(k*x) - cos(x)) )
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p * (cos(x) - exp(-k*x)) / (cosh(k*x) - cos(x)), 0, Inf)
gamma(p + 1) * pracma::zeta(p + 1) * (1/(k+1i)^(p+1) + 1/(k-1i)^(p+1))

###
m = 3;
integrate(\(x) x^p * (cos(m*x) - exp(-k*x)) / (cosh(k*x) - cos(m*x)), 0, Inf)
gamma(p + 1) * pracma::zeta(p + 1) * (1/(k+1i*m)^(p+1) + 1/(k-1i*m)^(p+1))

### I( x^p * sin(m*x) / (cosh(k*x) - cos(m*x)) )
m = exp(1);
integrate(\(x) x^p * sin(m*x) / (cosh(k*x) - cos(m*x)), 0, Inf)
-2 * gamma(p + 1) * pracma::zeta(p + 1) * Im(1/(k+1i*m)^(p+1))


#########################
#########################

### Trig - Exp

### I( cos(x) / (exp(1/x) + 1) )
# 1. Michael Penn: is this integration trick TOO POWERFUL?
#    https://www.youtube.com/watch?v=PX2QXILRgsc
# 2. Maths 505: This OP trick solves IMPOSSIBLE integrals!
#    https://www.youtube.com/watch?v=25dQUgpeX74
# Note: only trick;

###
lim = 1
integrate(\(x) cos(x) / (exp(1/x) + 1), -lim, lim)
sin(lim)

###
lim = 1; k = 1/5;
integrate(\(x) cos(k*x) / (exp(1/x) + 1), -lim, lim)
sin(k*lim) / k;

### Various Variants:
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) + 1), -lim, lim)
lim

###
lim = 1
integrate(\(x) cos(x)^2 / (exp(1/x) + 1), -lim, lim)
sin(2*lim) / 4 + lim/2

###
lim = sqrt(pi)
integrate(\(x) sin(x)^2 / (exp(1/x) + 1), -lim, lim)
lim/2 - sin(2*lim) / 4


### Variants: exp() - 1

###
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) - 1), -lim, lim)
- lim

###
lim = sqrt(pi)
integrate(\(x) cos(x) / (exp(1/x) - 1), -lim, lim)
- sin(lim)


#####################
#####################

### Clausen-Function:
# Integrals for Higher Order

# Ref:
# Tričković, S. B., & Stanković, M. S. (2023).
# On the closed form of Clausen functions.
# Integral Transforms and Special Functions, 34(6), 469–477.
# https://doi.org/10.1080/10652469.2022.2149961


### I( x^3 * ... )
n = 7; k = 3;
# Upper = Inf;
integrate(\(x) x^3 * exp(x) / (exp(2*x) - 2*cos(k/n*pi)*exp(x) + 1), 0, 40)
integrate(\(x) log(x)^3 / (x^2 - 2*cos(k/n*pi)*x + 1), 1, Inf)
id = seq(n); sn = sin(k*id*pi/n); idn = id / (2*n);
sum(sn * (pracma::psi(3, idn) + (-1)^k * pracma::psi(3, 1/2 + idn))) / (2*n)^4 / sin(k/n*pi);


# Note:
# - can be reduced to fraction,
#   but rescaling the fraction produces a non-classic interval;
#   (which is interesting in itself)
# - can be obtained in principle from
#   a polynomial fraction of order (2*n): x^p * P(x) / Q[2*n](x)
#   & 3x differentiation w. respect to p;

