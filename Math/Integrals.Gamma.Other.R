


Catalan = 0.915965594177219015054603514;
dzeta2  = -0.937548254316;
### Glaisher–Kinkelin Constant:
# https://en.wikipedia.org/wiki/Glaisher%E2%80%93Kinkelin_constant
A = exp((log(2*pi) + Euler - 6*dzeta2/pi^2)/12);



### History:

# [refactor] moved integrals:
#   I( x^p / (exp(x) + 1) )
#   I( x^p / (exp(x) - 1) )
#   to file: Integrals.Exp.R;


####################
####################

### I( exp(-x^p) / (x^p + 1) ) on [0, Inf]
# Dr. Peyam: But I AM joking, Mr. Feynman!
# https://www.youtube.com/watch?v=sdWpmwutSHE

###
integrate(\(x) exp(-x^2) / (x^2+1), 0, Inf)
erf1 = integrate(\(x) exp(-x^2), 0, 1)$value;
gamma(1/2) * (gamma(1/2)/2 - erf1) * exp(1)


###
integrate(\(x) exp(-x^3) / (x^3+1), 0, Inf)
erf1 = integrate(\(x) x*exp(-x^3), 0, 1)$value;
gamma(1/3) * (gamma(2/3)/3 - erf1) * exp(1)


### Generalization:
p = sqrt(5)
integrate(\(x) exp(-x^p) / (x^p+1), 0, Inf)
erf1 = integrate(\(x) x^(p-2)*exp(-x^p), 0, 1)$value;
gamma(1/p) * (gamma(1-1/p)/p - erf1) * exp(1)


########################
########################

### I( log(gamma(x + k)) )
# Flammable Maths:
# WOW! The Most AMAZING Way of Solving an Integral Ever! Deriving Raabe's Integral Formula!
# https://www.youtube.com/watch?v=APRr6wmLyuk

# see also:
# 1. Junesang Choi, H.M. Srivastava.
#    A family of log-gamma integrals and associated results.
#    J. Math. Anal. Appl. 303 (2005)


###
k = sqrt(5)
integrate(\(x) log(gamma(x)), k, k+1)
integrate(\(x) log(gamma(x + k)), 0, 1)
log(2*pi)/2 + k*log(k) - k


### I( log(gamma(2*x + k)) )
k = sqrt(5)
integrate(\(x) log(gamma(2*x + k)), 0, 1)
((k+1)*log(k+1) + k*log(k)) / 2 - k + log(2*pi)/2 - 1/2

### I( log(gamma(3*x + k)) )
k = sqrt(5)
integrate(\(x) log(gamma(3*x + k)), 0, 1)
(k*log(k) + (k+1)*log(k+1) + (k+2)*log(k+2)) / 3 - k + log(2*pi)/2 - 1


### I( log(gamma(x/k1 + k2)) )
k = sqrt(5)
integrate(\(x) log(gamma(x/2 + k)), 0, 1)
# TODO


### I( log(gamma(x/2 + 1)) )
# see article;
integrate(\(x) 1/2 * log(gamma(x/2 + 1)), 0, 1)
integrate(\(x) log(gamma(x + 1)), 0, 1/2)
1/4 * log(pi) + 3/2 * log(A) - 7/24 * log(2) - 1/2
3/8 * log(pi) - (6*dzeta2/pi^2 - Euler) / 8 - 1/6 * log(2) - 1/2
#
integrate(\(x) 1/2 * log(gamma(x/2 + 3/2)), 0, 1)
integrate(\(x) log(gamma(x + 1)), 1/2, 1)
log(pi)/8 + (6*dzeta2/pi^2 - Euler) / 8 + 2/3 * log(2) - 1/2

###
integrate(\(x) 1/4 * log(gamma(x/4 + 1)), 0, 1)
integrate(\(x) log(gamma(x + 1)), 0, 1/4)
log(pi)/8 + 9/8 * log(A) - 3/8 * log(2) + Catalan / (4*pi) - 1/4

### on [0, z]
z = 1/sqrt(7);
integrate(\(x) log(gamma(x + 1 - z/2)), 0, z)
integrate(\(x) - 1/pi * log(sin(x)), 0, pi/2 * z)$value +
	+ log(pi/2 * z) * z/2 - z/2;


###############

### I( x * log(gamma(x)) )
integrate(\(x) x * log(gamma(x)), 0, 1)
log(2*pi)/4 - log(A);


### TODO:
integrate(\(x) log(1-x) * log(gamma(x)), 0, 1)
integrate(\(x) log(x) * log(gamma(x)), 0, 1)

integrate(\(x) log(gamma(x)) / (x + 1), 0, 1)
integrate(\(x) log(gamma(x)) / (x^2 + 1), 0, 1)
integrate(\(x) log(gamma(x)) / log(x), 0, 1)
integrate(\(x) log(gamma(x)) / gamma(x), 0, 1)
integrate(\(x) log(gamma(x)) * log(gamma(1-x)), 0, 1)

integrate(\(x) 1/gamma(x), 0, 1)
integrate(\(x) x / gamma(x), 0, 1)
integrate(\(x) log(x) / gamma(x), 0, 1)
integrate(\(x) sin(x) / gamma(x), 0, 1)

integrate(\(x) gamma(x) * x, 0, 1)
integrate(\(x) gamma(x) * sin(x), 0, 1)
integrate(\(x) gamma(x) * sin(pi*x), 0, 1)
integrate(\(x) gamma(x) * (exp(x) - 1), 0, 1)
integrate(\(x) gamma(x) * log(x + 1), 0, 1)
integrate(\(x) gamma(x) * log(x*exp(x) + 1), 0, 1)
integrate(\(x) gamma(x) * atan(x), 0, 1)
integrate(\(x) atan(gamma(x)), 0, 1)
integrate(\(x) gamma(x) / gamma(1/x), 0,1)
integrate(\(x) gamma(x) / pracma::zeta(x+1), 0, 1)
integrate(\(x) gamma(x) / pracma::zeta(1-x), 0, 1)


###
integrate(\(x) x * log(gamma(x)*gamma(1-x)), 0, 1)
log(2*pi)/2

integrate(\(x) x * log(sin(x)), 0, pi)
- pi^2 * log(2)/2
integrate(\(x) x * log(sin(pi*x)), 0, 1)
- log(2)/2


