#########################
##
## Leonard Mada
##
## Integrals: Exp
## Various Combinations


####################

### Helper Functions

Euler   = 0.57721566490153286060651209008240243079;


####################
####################

### Zeta Function

### I( zeta(1 + exp(1i*x)) )
# Maths 505: This is one of the coolest integrals ever solved
# https://www.youtube.com/watch?v=YYhckN554v8
# Residue Theorem: Pole of Order 2 (1/z^2);

integrate(\(x) Re(pracma::zeta(1 + exp(1i*x))), 0, 2*pi)
2*pi*Euler

integrate(\(x) Im(pracma::zeta(1 + exp(1i*x))), 0, 2*pi)
0;


### I( zeta(2 + 2*exp(1i*x)) )
integrate(\(x) Re(pracma::zeta(2 + 2*exp(1i*x))), 0, 2*pi, rel.tol=1E-9)
2*pi*(pi^2/6 - 1)

integrate(\(x) Im(pracma::zeta(2 + 2*exp(1i*x))), 0, 2*pi)
0;


### I( zeta(2 + exp(1i*x)) )
integrate(\(x) Re(pracma::zeta(2 + exp(1i*x))), 0, 2*pi)
2*pi*(pi^2/6 - 1) + pi;
pi^3/3 - pi;

#
up = 2*pi - 1E-9; # Numerical instability;
integrate(\(x) Im(pracma::zeta(2 + exp(1i*x))), 0, up)
0;


### I( zeta(3 + n*exp(1i*x)) )

###
integrate(\(x) Re(pracma::zeta(3 + 1*exp(1i*x))), 0, 2*pi, rel.tol=1E-9)
2*pi*pracma::zeta(3)

###
integrate(\(x) Re(pracma::zeta(3 + 2*exp(1i*x))), 0, 2*pi, rel.tol=1E-9)
2*pi*pracma::zeta(3) - pi/2

###
integrate(\(x) Re(pracma::zeta(3 + 3*exp(1i*x))), 0, 2*pi, rel.tol=1E-9)
2*pi*pracma::zeta(3) - pi


####################
####################

### I( Ei(-x^2) * exp(-x^2) )
# Hmath: integral with integral exponential function
# https://www.youtube.com/watch?v=0pP94DgWdds
# [in Russian]
# Ei => Integral form => Double I() => polar coordinates;

### I( Ei(-x^2) * exp(-x^2) )
integrate(\(x) Re(pracma::expint(x^2)) * exp(-x^2), 0, Inf, rel.tol=1E-13)
log(1+sqrt(2)) * sqrt(pi);

### I( Ei(-x) * exp(-x) )
integrate(\(x) Re(pracma::expint(x)) * exp(-x), 0, Inf, rel.tol=1E-13)
log(2);

### I( Ei(-x^3) * exp(-x^3) )
integrate(\(x) Re(pracma::expint(x^3)) * exp(-x^3), 0, Inf, rel.tol=1E-13)
gamma(1/3) * integrate(\(x) 1 / (x^3 + 1)^(1/3), 0, 1, rel.tol=1E-13)$value;
# see file: Integrals.Fractions.Unity.Radicals.R;
gamma(1/3) * (
	- (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan(1/sqrt(3) * (2^(1/3) - 1) / (2^(1/3) + 1)));

### I( Ei(-x^4) * exp(-x^4) )
integrate(\(x) Re(pracma::expint(x^4)) * exp(-x^4), 0, Inf, rel.tol=1E-13)
gamma(1/4) * integrate(\(x) 1 / (x^4 + 1)^(1/4), 0, 1, rel.tol=1E-13)$value;
# TODO

