##################
##
## Integrals: Trig
## Polynomial Fractions
##
## Leonard Mada


## Note:
# - [refactor] moved Basic Fractions to:
#   Integrals.Trig.Fractions.Simple.R;
# - [refactor] moved Inverse Fractions to:
#   Integrals.Trig.Fractions.Inverse.R;


####################

### Helper Functions


Catalan = 0.915965594177219015054603514;


int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/(2*n))) / (2*n);
}
int.FrDU01 = function(n, p=0) {
	# TODO: p != 0
	digamma(1/n) + Euler + log(n);
}

zeta = function(x) {
	pracma::zeta(x);
}


####################
####################

### I( sin(x) / (x + 1) )
# Michael Penn: Check out this crazy integral trick!
# https://www.youtube.com/watch?v=nBGcO3VL8Kw

# Upper = Inf: numeric instability;
pracma::integral(\(x) sin(x) / (x + 1), 0, 1000000)
integrate(\(x) 1 / (log(x)^2 + 1), 0, 1)
integrate(\(x) exp(-x) / (x^2 + 1), 0, Inf)
# TODO: ???


####################
####################

### I( sin(x) / (x^2 + 1) )
# Various variants:

# 1) qncubed3:  TWO WAYS to destroy this INSANE integral feat. @maths_505
#    https://www.youtube.com/watch?v=Ry5q4NsZDx0
# 2) Maths 505: INSANE integral solved using 2 different methods
#    (feat. Feynman's technique and Lobachevsky)
#    https://www.youtube.com/watch?v=GMcCiR0y3wg


###
integrate(\(x) sin(x)^2 / (x^2*(x^2+1)), 0, Inf)
pi/4*(1 + exp(-2))


### Helper:
k = sqrt(3)
integrate(\(x) cos(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2 * exp(-k)


### I( x * sin(k*x) / (x^2 + 1) )
k = sqrt(3)
# Upper = Inf: numerical instability;
pracma::integral(\(x) x * sin(k*x) / (x^2 + 1), 0, 200000)
pi/2 * exp(-k)


### Gen 1:
k = sqrt(3)
b = sqrt(5)
integrate(function(x) cos(k*x) / (x^2 + b^2), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2 * exp(-k*b)/b


### TODO
k = 1
integrate(\(x) sin(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-5)
integrate(\(x) sapply(x, \(x) {
	1/(4*pi) * Im(sin(k*x) * (pracma::psi((x + 1i)/(2*pi)) - pracma::psi((x - 1i)/(2*pi))))
	}), 0, 2*pi)
(exp(2)*pracma::expint(1) - pracma::expint(-1) - 1i*pi) / (2*exp(1))
# TODO

###
k = sqrt(5)
integrate(\(x) sin(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-5)
integrate(\(x) sin(k*x) / (x^2 + 1), 0, 2000, subdivisions=4096, rel.tol=1E-6)
(exp(2*k)*pracma::expint(k) - pracma::expint(-k) - 1i*pi) / (2*exp(k))


### I( x * cos(x) / (x^2 + 1) )
# upper = Inf; # numerical instability!
integrate(\(x) x * cos(x) / (x^2 + 1), 0, 1E+5, subdivisions = 8000)
- Re(exp(2)*pracma::expint_Ei(-1) + pracma::expint_Ei(1)) / (2*exp(1))
(exp(2)*pracma::expint_E1(1) - pracma::expint_Ei(1)) / (2*exp(1))
(exp(2)*pracma::expint(1) + pracma::expint(-1) + 1i*pi) / (2*exp(1))


### Gen 1:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2 + 1), 0, Inf, subdivisions=4096, rel.tol=1E-6)
pi/4*(1 - exp(-2*k))


### Gen 2:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2*(x^2 + 1)), 0, Inf, subdivisions=1024, rel.tol=1E-8)
pi/4*(2*k - 1 + exp(-2*k))

# f = (1 - exp(2i*k*z)) / (z^2 * (z^2 + 1));
# Res at 1i = pi*(1 - exp(-2*k));


### I( sin(k*x) / (x*(x^2 + 1)) )
# D(k) =>
k = sqrt(3)
# Upper = Inf: Numeric instability;
pracma::integral(function(x) sin(k*x) / (x*(x^2 + 1)), 0, 100000)
pi/2*(1 - exp(-k))

### I( sin(tan(x)) / x )
# takes very long: upper = Inf
pracma::integral(\(x) sin(tan(x)) / x, 0, 2000)
pracma::integral(\(x) sin(x) / (x*(x^2 + 1)), 0, 10000)
pi/2*(1 - exp(-1))


### I( sin(k*x)^3 / (x*(x^2 + 1)) )
pracma::integral(function(x) sin(k*x)^3 / (x*(x^2 + 1)), 0, 100000)
pi/8*(2 - 3*exp(-k) + exp(-3*k))


### Helper:

### I( cos(k*x)/x - exp(-x)/x )
k = 3
# upper = Inf; numerical issues;
pracma::integral(\(x) cos(k*x)/x - exp(-x)/x, 0, 100000, no_=1024)
- log(k)

###
k = 3
# upper = Inf; numerical issues;
pracma::integral(\(x) cos(k*x)/x^2 - 1/x^2, 0, 200000)
- k*pi/2


########################

### Numerical Stability:

k = sqrt(3)
start = 0; steps = 1000
r = sapply(seq(start, steps), function(id) {
	startI = 2*pi*id/k;
	endI = 2*pi*(id+1)/k;
	integrate(function(x) cos(k*x) / (x^2 + 1), startI, endI, subdivisions=4096*2, rel.tol=1E-8)$value
})
r = sum(r)
# FAILS miserably!
r + integrate(function(x) cos(k*x) / (x^2 + 1), 2*pi*(steps+1)/k, Inf, subdivisions=4096*2, rel.tol=1E-8)$value
# residual: should be < 1E-8;
print(r, 12)
print(pi/2*exp(-k), 12)


#######################
#######################

### I( x / (1 - cos(x)) )
integrate(\(x) x / (1 - cos(x)) - 2/x, 0, pi)
2*log(2) + 2 - 2*log(pi);

### on [0, pi/2]
integrate(\(x) x / (1 - cos(x)) - 2/x, 0, pi/2)
3*log(2) - pi/2 + 2 - 2*log(pi)


# Base:
integrate(\(x) x / (1 - cos(x)) - 2/x + 2*log(pi)/pi, 0, pi)
2*log(2) + 2

integrate(\(x) x / (1 - cos(x)) - 2/x + 4*log(pi/2)/pi, 0, pi/2)
log(2) - pi/2 + 2


#######################
#######################

### ATAN

### I( atan(x^p) / x )
integrate(\(x) atan(x) / x, 0, 1)
Catalan
#
integrate(\(x) atan(x^2) / x, 0, 1)
Catalan / 2
#
p = sqrt(5)
integrate(\(x) atan(x^p) / x, 0, 1)
Catalan / p

### I( x^p * atan(x^n) )
n = sqrt(5)
integrate(\(x) atan(x^n), 0, 1)
pi/4 - n * int.FrU01(2*n, n)

###
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * atan(x^n), 0, 1)
pi/(4*(p+1)) +
	- (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1));
pi/(4*(p+1)) - n / (p+1) * int.FrU01(2*n, n + p);


###
# see file: Integrals.Trig.Tan.R;
integrate(\(x) x * atan(x) / (x^2 + 1), 0, 1)
Catalan / 2 - 1/8 * pi*log(2)

###
integrate(\(x) atan(x) / (x^2 + 1), 0, 1)
pi^2 / 32


### I( x^p * atan(x^n) * log(x) )
n = sqrt(3); p = sqrt(5);
# works only for n > 0;
# see also file: Integrals.Trig.Tan.R; # TODO: check which file;
integrate(\(x) x^p * atan(x^n) * log(x), 0, 1)
- pi/(4*(p+1)^2) +
	+ (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1)^2) +
	- (pracma::psi(1, ((p+1)/n + 3)/4) - pracma::psi(1, ((p+1)/n + 1)/4)) / (16*(p+1)*n);

