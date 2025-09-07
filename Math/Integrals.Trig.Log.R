


### I( sin(k*x) * log(x) / x^(p+1) )

###
p = sqrt(3) - sqrt(2)
# Up = Inf; numerical issues;
integrate(\(x) sin(x) * log(x) / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-p) * digamma(-p) * sin(pi*p/2) + pi/2 * gamma(-p) * cos(pi*p/2)


###
p = sqrt(3) - sqrt(2)
k = sqrt(3)
# Up = Inf; numerical issues;
integrate(\(x) sin(k*x) * log(x) / x^(p+1), 0, 100000, subdivisions=40000)
- gamma(-p) * digamma(-p) * sin(pi*p/2) * k^p +
	+ pi/2 * gamma(-p) * cos(pi*p/2) * k^p + gamma(-p) * sin(pi*p/2) * log(k) * k^p;


### Base:
# p in (-1, 1)
# see file: Integrals.Trig.Fractions.Poly.R;
p = 2 - sqrt(2)
k = sqrt(3)
integrate(\(x) sin(k*x) / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-p) * sin(pi*p/2) * k^p;


###################
###################

### I( cos(x) * log(x^2 + k^2) ) on [0, Inf]
# Maths 505: A savage log trig integral!
# https://www.youtube.com/watch?v=Jt1UVLnd1bA
# Note: Feynman trick on log(x^2 + b^2);
# - misses a factor of 2 => Const = 0;

# upper = Inf;
integrate(\(x) cos(x) * log(x^2+1), 0, 1024*pi, subdivisions = 1024, rel.tol=1E-8)
integrate(\(x) cos(x) * log(x^2+1), 1024*pi, 2*1024*pi, subdivisions = 1024, rel.tol=1E-8)
- pi * exp(-1);

### I( cos(x) * log(x^2+k^2) )
k = sqrt(3)
# upper = Inf;
integrate(\(x) cos(x) * log(x^2+k^2), 0, 1024*pi, subdivisions = 1024, rel.tol=1E-8)
integrate(\(x) cos(x) * log(x^2+k^2), 1024*pi, 2*1024*pi, subdivisions = 1024, rel.tol=1E-8)
- pi * exp(-k);

# Helper
k = sqrt(5)
# upper = Inf;
integrate(\(x) cos(k*x) / (x^2+1), 0, 2^16, subdivisions = 16024, rel.tol=1E-10)
pi/2 * exp(-k);

#
k = sqrt(5)
integrate(\(x) exp(-k*x) * log(x), 0, Inf, subdivisions = 16024, rel.tol=1E-10)
- (log(k) + Euler) / k;

### Bug:
k = sqrt(5)
# This must be a bug:
integrate(\(x) exp(-k*x) * log(x), 0, 2^16, subdivisions = 16024, rel.tol=1E-10)
# 6.527367e-67 with absolute error < 1.3e-66

integrate(\(x) exp(-k*x) * log(x), 0, Inf, subdivisions = 16024, rel.tol=1E-10)
# -0.61802 with absolute error < 7.8e-11
- (log(k) - digamma(1)) / k;
# -0.61802

integrate(\(x) exp(-k*x) * log(x), 0, 1, subdivisions = 16024, rel.tol=1E-10)
# -0.6338608 with absolute error < 3.9e-13
integrate(\(x) exp(-k*x) * log(x), 0, 2^12, subdivisions = 16024, rel.tol=1E-10)
# -0.61802 with absolute error < 6e-11


#################
#################

### I( cos(p * log(x)) / ... )

# Maths 505: A nice integral result
# https://www.youtube.com/watch?v=b-5f66DSD7g

# Generalization:

p = sqrt(3); n = 5*sqrt(11); k = 1/sqrt(2);
integrate(\(x) cos(p*log(x)) / (x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((1i*p+1)/n) * pracma::gammaz(k - (1i*p+1)/n) +
	pracma::gammaz((-1i*p+1)/n) * pracma::gammaz(k - (-1i*p+1)/n) ) / (2*n*gamma(k));

### I( x^p * cos(q * log(x)) / ... )
p = sqrt(2); q = sqrt(3);
n = sqrt(5); k = sqrt(7);
integrate(\(x) x^p * cos(q*log(x)) / (x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((p+1 + 1i*q)/n) * pracma::gammaz(k - (p+1 + 1i*q)/n) +
	+ pracma::gammaz((p+1 - 1i*q)/n) * pracma::gammaz(k - (p+1 - 1i*q)/n) ) /
	(2*n*gamma(k));

### I( x^p * sin(q * log(x)) / ... )
p = sqrt(2); q = sqrt(3);
n = sqrt(5); k = sqrt(7);
integrate(function(x) x^p * sin(q*log(x)) /(x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((p+1 + 1i*q)/n) * pracma::gammaz(k - (p+1 + 1i*q)/n) +
	- pracma::gammaz((p+1 - 1i*q)/n) * pracma::gammaz(k - (p+1 - 1i*q)/n) ) /
	(2i*n*gamma(k));


### I( x^p * cos(q * log(x^n + 1)) / (x^n + 1)^k )
p = 1/sqrt(3); q = sqrt(pi); n = sqrt(5); k = sqrt(2);
integrate(function(x) x^p * cos(q*log(x^n + 1)) / (x^n + 1)^k, lower=0, upper=Inf)
Re(pracma::gammaz(k - 1i*q - (p+1)/n) / pracma::gammaz(k - 1i*q)) *
	gamma((p+1)/n) / n;


integrate(function(x) x^p * sin(q * log(x^n + 1)) / (x^n + 1)^k, lower=0, upper=Inf)
Im(pracma::gammaz(k - 1i*q - (p+1)/n) / pracma::gammaz(k - 1i*q)) *
	gamma((p+1)/n) / n;


### I( x^p * cos(q * log(x / (x^n + 1))) / (x^n + 1)^k )
p = 1/sqrt(3); q = sqrt(pi); n = sqrt(5); k = sqrt(2);
integrate(function(x) x^p * cos(q*log(x / (x^n + 1))) / (x^n + 1)^k,
	lower=0, upper=Inf)
Re(pracma::gammaz(k - 1i*q - (p+1 - 1i*q)/n) / pracma::gammaz(k - 1i*q) *
	pracma::gammaz((p+1 - 1i*q)/n)) / n;


p = 1/sqrt(3); q = sqrt(pi); n = sqrt(5); k = sqrt(2);
integrate(function(x) x^p * sin(q*log(x / (x^n + 1))) / (x^n + 1)^k,
	lower=0, upper=Inf)
- Im(pracma::gammaz(k - 1i*q - (p+1 - 1i*q)/n) / pracma::gammaz(k - 1i*q) *
	pracma::gammaz((p+1 - 1i*q)/n)) / n;


### Power 2: I( Trig(p * log(x))^2 / ... )

### cos-Variant
p = sqrt(3); n = 5*sqrt(11); k = 1/sqrt(2);
integrate(function(x) cos(p*log(x))^2 / (x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((2i*p+1)/n) * pracma::gammaz(k - (2i*p+1)/n) +
	+ pracma::gammaz((-2i*p+1)/n) * pracma::gammaz(k - (-2i*p+1)/n) +
	+ 2*gamma(1/n) * gamma(k - 1/n) ) / (4*n*gamma(k));

### sin-Variant
p = sqrt(3); n = 5*sqrt(11); k = 1/sqrt(2);
integrate(function(x) sin(p*log(x))^2 / (x^n + 1)^k, lower=0, upper=Inf)
- (pracma::gammaz((2i*p+1)/n) * pracma::gammaz(k - (2i*p+1)/n) +
	+ pracma::gammaz((-2i*p+1)/n) * pracma::gammaz(k - (-2i*p+1)/n) +
	- 2*gamma(1/n) * gamma(k - 1/n) ) / (4*n*gamma(k));

