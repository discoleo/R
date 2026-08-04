

#################
### ATAN: Exp ###

# Note:
# - more atan() / EXP variants are in file:
#   Integrals.Exp.Advanced.R;


### I( atan(exp(-x))/x )
integrate(\(x) atan(exp(-x))/x - pi/4*exp(-x)/x, 0, Inf)
pi/4 * log(4*pi^3) - pi*log(gamma(1/4)) + pi*Euler/4


#########################

### I( x^2 * cosh(atan(x)) / sqrt(x^2 + 1) )
# Maths 505: A RIDICULOUSLY AWESOME INTEGRAL:
# int 0 to 1 (x^2cosh(arctan(x)))/sqrt(1+x^2)
# https://www.youtube.com/watch?v=u69PYnJzQvU
# Note: even function;

integrate(\(x) x^2 * cosh(atan(x)) / sqrt(x^2 + 1), 0, 1)
integrate(\(x) 1/2 * x^2 * exp(atan(x)) / sqrt(x^2 + 1), -1, 1)
integrate(\(x) 1/2 * sin(x)^2 / cos(x)^3 * exp(x), -pi/4, pi/4)
exp(-pi/4) * sqrt(2)/2;

### I( x^2 * sinh(atan(x)) / sqrt(x^2 + 1) )
integrate(\(x) x^2 * sinh(atan(x)) / sqrt(x^2 + 1), 0, 1)
1/2 - exp(-pi/4) * sqrt(2)/2;


### I( cosh(atan(x)) / sqrt(x^2 + 1) )
integrate(\(x) cosh(atan(x)) / sqrt(x^2 + 1), 0, 1)
integrate(\(x) 1/2 / cos(x) * exp(x), -pi/4, pi/4)
integrate(\(x) cosh(x) / cos(x), 0, pi/4)
# TODO


### I( x * sinh(atan(x)) / sqrt(x^2 + 1) )
integrate(\(x) x * sinh(atan(x)) / sqrt(x^2 + 1), 0, 1)
integrate(\(x) - cosh(x) / cos(x), 0, pi/4)$value + sinh(pi/4) * sqrt(2);
# TODO


#####################

### I( atan(atan(x) / atan(up-x)) ) on [0,up]
lim = exp(-1);
integrate(\(x) atan(atan(x) / atan(lim-x)), 0, lim, rel.tol=1E-13)
pi/4 * lim;

### I( atan(atan(x) / atan(1-x)) )
integrate(\(x) atan(atan(x) / atan(1-x)), 0, 1, rel.tol=1E-13)
pi/4;

### I( atan(atan(x) / atan(2-x)) ) on [0,2]
integrate(\(x) atan(atan(x) / atan(2-x)), 0, 2, rel.tol=1E-13)
pi/2;


### I( atan(log(x) / log(1-x)) )
up = exp(-1); # up <= 1;
integrate(\(x) atan(log(x) / log(up-x)), 0, up, rel.tol=1E-13)
pi/4 * up;

#
integrate(\(x) atan(log(x) / log(1-x)), 0, 1, rel.tol=1E-13)
pi/4;

#
integrate(\(x) atan(log(x) / log(1/2-x)), 0, 1/2, rel.tol=1E-13)
pi/8;

#
integrate(\(x) atan(log(x) / log(1/3-x)), 0, 1/3, rel.tol=1E-13)
pi/12;

