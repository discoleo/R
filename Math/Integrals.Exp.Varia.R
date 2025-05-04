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

