


Euler   = 0.577215664901532860606512090;
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
integrate(\(x) exp(-x^2) / (x^2 + 1), 0, Inf)
erf1 = integrate(\(x) exp(-x^2), 0, 1)$value;
gamma(1/2) * (gamma(1/2)/2 - erf1) * exp(1)


###
integrate(\(x) exp(-x^3) / (x^3 + 1), 0, Inf)
erf1 = integrate(\(x) x*exp(-x^3), 0, 1)$value;
gamma(1/3) * (gamma(2/3)/3 - erf1) * exp(1)


### Generalization:
p = sqrt(5)
integrate(\(x) exp(-x^p) / (x^p + 1), 0, Inf)
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

### Lim: k -> 0
integrate(\(x) log(gamma(2*x)), 0, 1)
log(2*pi)/2 - 1/2
###
integrate(\(x) log(gamma(3*x)), 0, 1)
log(2*pi)/2 + 2/3*log(2) - 1


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
3/8 * log(pi) - (6*dzeta2/pi^2 - Euler) / 8 - 1/6 * log(2) - 1/2;
#
integrate(\(x) 1/2 * log(gamma(x/2 + 3/2)), 0, 1)
integrate(\(x) log(gamma(x + 1)), 1/2, 1)
log(pi)/8 + (6*dzeta2/pi^2 - Euler) / 8 + 2/3 * log(2) - 1/2;

### on [0, 1/2]
integrate(\(x) log(gamma(x)), 0, 1/2)
3/8 * log(pi) - (6*dzeta2/pi^2 - Euler) / 8 + log(2) / 3;


###
integrate(\(x) 1/4 * log(gamma(x/4 + 1)), 0, 1)
integrate(\(x) log(gamma(x + 1)), 0, 1/4)
log(pi)/8 + 9/8 * log(A) - 3/8 * log(2) + Catalan / (4*pi) - 1/4

### on [0, z]
z = 1/sqrt(7);
integrate(\(x) log(gamma(x + 1 - z/2)), 0, z)
integrate(\(x) z * log(gamma(z*x + 1 - z/2)), 0, 1)
integrate(\(x) - 1/pi * log(sin(x)), 0, pi/2 * z)$value +
	+ log(pi/2 * z) * z/2 - z/2;

# Note:
# I( log(sin(x)) ) is solvable for z = rational fraction,
# see file: Integrals.Log.Trig.R;

### on [0, 1/2]
integrate(\(x) log(gamma(x + 3/4)), 0, 1/2)
Catalan / (2*pi) + log(pi/2) /4 - 1/4


### on [0, 1/4] & [1/4, 1/2]
integrate(\(x) log(gamma(x/2 + 3/4)), 0, 1/2)
integrate(\(x) 2 * log(gamma(x + 3/4)), 0, 1/4)
integrate(\(x) - 2/pi * log(sin(x)), 0, pi/4)$value +
	+ log(pi) / 4 - log(2)/4 - 9/4 * log(A) - Catalan / (2*pi);
#
integrate(\(x) log(gamma(x + 3/4)), 1/4, 1/2)
log(pi) / 8 - 3/8 * log(2) + 9/8 * log(A) + Catalan / (4*pi) - 1/4;


# Derivation:
integrate(\(x) log(gamma(x/2 + 3/4)), 0, 1)
integrate(\(x) - 2/pi * log(sin(x)), 0, pi/4)$value +
	+ log(pi/4) /2 - 1/2;
#
integrate(\(x) log(gamma(x/2 + 3/4)), 1/2, 1)
log(pi)/4 + 9/4 * log(A) - 3/4 * log(2) + Catalan / (2*pi) - 1/2;


###############

### I( x * log(gamma(x)) )
integrate(\(x) x * log(gamma(x)), 0, 1)
log(2*pi)/4 - log(A);


### I( cos(x) * log(gamma(x/(2*pi))) )
integrate(\(x) cos(x) * log(gamma(x/(2*pi))), 0, 2*pi)
pi / 2

### I( sin(x) * log(gamma(x/(2*pi))) )
integrate(\(x) sin(x) * log(gamma(x/(2*pi))), 0, 2*pi)
log(2*pi) + Euler

###
integrate(\(x) cos(x) * log(gamma(x/pi)), 0, 2*pi)
integrate(\(x) - cos(x) * log(x/pi), 0, pi)
integrate(\(x) sin(pi*x) / x, 0, 1)
# TODO


### I( sin(x) * digamma(x/(2*pi)) )
integrate(\(x) sin(x) * digamma(x/(2*pi)), 0, 2*pi)
- pi^2

### I( cos(x) * digamma(x/(2*pi)) )
integrate(\(x) cos(x) * digamma(x/(2*pi)) + 2*pi/x, 0, 2*pi)
2*pi*Euler + 2*pi*log(2*pi);


### I( digamma(x/(2*pi)) )
integrate(\(x) digamma(x/(2*pi)) + 2*pi/x, 0, 1)
2*pi * (log(gamma(1/(2*pi))) - log(2*pi))

### I( digamma(k*x) )
k = 1 / sqrt(pi)
integrate(\(x) digamma(k*x) + 1 / (k*x), 0, 1)
(log(gamma(k)) + log(k)) / k;


### I( digamma(1 + 1i*x) )
integrate(\(x) sapply(x, \(x) Im(pracma::psi(0, 1 + 1i*x))), 0, 1)
log(sinh(pi)/pi) / 2;
#
integrate(\(x) sapply(x, \(x) Re(pracma::psi(0, 1 + 1i*x))), 0, 1)
- atan(Re(pracma::gammaz(1i)) / Im(pracma::gammaz(1i)))


### I( x * digamma(x) )
integrate(\(x) x * digamma(x), 0, 1)
- log(pi*2)/2

###
k = 3
integrate(\(x) x * digamma(k*x), 0, 1)
log(gamma(k))/k - integrate(\(x) 1/k * log(gamma(k*x)), 0, 1)$value
# TODO


### I( x * psi(1, k*x) )
k = 1/sqrt(pi)
integrate(\(x) x * pracma::psi(1, k*x) - 1 / (k^2*x), 0, 1)
(k*digamma(k) - log(gamma(k)) - log(k) + 1) / k^2


###############

### I( cos(pi*x)^2 * log(gamma(x)) )
# Maths 505: This integral is INSANE beyond measure!
# https://www.youtube.com/watch?v=Em5R4ckyqk0
# Note: uses Reflection of Gamma;

### I( cos(pi*x)^2 * log(gamma(x)) )
integrate(\(x) cos(pi*x)^2 * log(gamma(x)), 0, 1)
(log(2*pi) + 1/2)/4


### I( sin(pi*x)^2 * log(gamma(x)) )
integrate(\(x) sin(pi*x)^2 * log(gamma(x)), 0, 1)
(log(2*pi) - 1/2)/4


### I( sin(pi*x) * log(gamma(x)) )
integrate(\(x) sin(pi*x) * log(gamma(x)), 0, 1)
log(pi/2) / pi + 1/pi;


### I( cos(pi*x) * log(gamma(x)) )
integrate(\(x) cos(pi*x) * log(gamma(x)), 0, 1)
# TODO


### I( sin(2*pi*x) * log(gamma(x)) )
integrate(\(x) sin(2*pi*x) * log(gamma(x)), 0, 1)
(log(2*pi) + Euler) / (2*pi);

### I( cos(2*pi*x) * log(gamma(x)) )
integrate(\(x) cos(2*pi*x) * log(gamma(x)), 0, 1)
1/4;


### I( log(gamma(x)) * tan(pi/2*x) )
# Note: Reflection does NOT work;
integrate(\(x) log(gamma(x)) * tan(pi/2*x), 0, 1)
# TODO

### I( log(gamma(x)) * tan(pi*x) )
# TODO

# library(Rmpfr)
FUN = \(x) {
	x = mpfr(x, 240);
	y = log(gamma(x)) * tan(pi*x);
	as.numeric(y);
}
eps = 1E-7;
integrate(FUN, 0, 1/2 - eps, rel.tol=1E-9)$value +
integrate(FUN, 1/2 + eps, 1, rel.tol=1E-9)$value;


#################

### I( atan(gamma(x) / gamma(1-x)) )
integrate(\(x) atan(gamma(x) / gamma(1-x)), 0, 1)
pi/4;


##################

### I( sin(pi*x) * gamma(x) )
integrate(\(x) sin(pi*x) * gamma(x), 0, 1)
integrate(\(x) pi / gamma(x), 0, 1)
# TODO


##################

### I( |Gamma(1/2 + x*1i)|^2 )
# Hmath: Integral with the gamma function of a complex argument
# [in Russian]
# https://www.youtube.com/watch?v=AnjZQjTgHG8
# Note:
# abs(gamma(1/2 + 1i*x)) = sqrt(gamma(1/2 + ...) * gamma(1/2 - ...))
#  = sqrt(pi / sin(pi*(1/2 + 1i*x))) = sqrt(pi / cosh(pi*1i*x));

integrate(\(x) abs(pracma::gammaz(1/2 + x*1i))^2, 0, Inf)
pi/2;


### I( |Gamma(1/2 + x*1i)| )
integrate(\(x) abs(pracma::gammaz(1/2 + x*1i)), 0, Inf, rel.tol=1E-13)
beta(1/4, 1/2) / sqrt(4*pi);
gamma(1/4) / gamma(3/4) / 2;

### I( Gamma(1/2 + x*1i) )
integrate(\(x) Re(pracma::gammaz(1/2 + x*1i)), 0, Inf, rel.tol=1E-13)
pi / exp(1);

#
integrate(\(x) Im(pracma::gammaz(1/2 + x*1i)), 0, Inf, rel.tol=1E-13)
# TODO


##################

### I( log(zeta(x) / zeta(1-x)) )
# Math 505: A mathematical treat
# https://www.youtube.com/watch?v=N2J0kmQ0xSc
# Reflection of zeta => ...;

integrate(\(x) log(pracma::zeta(x) / pracma::zeta(1-x)), 0, 1/2)
integrate(\(x) x*log(2) + (x-1)*log(pi) + log(sin(pi*x/2) * gamma(1-x)), 0, 1/2)
- (log(pi)/8 + Catalan/pi + 3/2*log(A) + log(2)/12);


##############
##############

### TODO:
integrate(\(x) log(1-x) * log(gamma(x)), 0, 1)
integrate(\(x) log(x) * log(gamma(x)), 0, 1)

integrate(\(x) log(gamma(x)) * sin(pi*x), 0, 1) # OK
integrate(\(x) log(gamma(x)) / (x + 1), 0, 1)
integrate(\(x) log(gamma(x)) / (1 - x), 0, 1)
integrate(\(x) log(gamma(x)) / (x^2 + 1), 0, 1)
integrate(\(x) log(gamma(x)) / log(x), 0, 1)
integrate(\(x) log(gamma(x)) / gamma(x), 0, 1)
integrate(\(x) log(gamma(x)) * log(gamma(1-x)), 0, 1)

integrate(\(x) 1 / gamma(x), 0, 1)
integrate(\(x) x / gamma(x), 0, 1)
integrate(\(x) log(x) / gamma(x), 0, 1)
integrate(\(x) sin(x) / gamma(x), 0, 1)
# Gamma(x) * ...
integrate(\(x) gamma(x) - 1/x, 0, 1)
integrate(\(x) gamma(x) / x - 1/x^2 + Euler/x, 0, 1)
integrate(\(x) gamma(x) * x, 0, 1)
integrate(\(x) gamma(x) * sin(x), 0, 1)
integrate(\(x) gamma(x) * sin(pi*x), 0, 1)
integrate(\(x) gamma(x) * (exp(x) - 1), 0, 1)
integrate(\(x) gamma(x) * log(x + 1), 0, 1)
integrate(\(x) gamma(x) * log(x*exp(x) + 1), 0, 1)
integrate(\(x) gamma(x) * atan(x), 0, 1)
integrate(\(x) atan(gamma(x)), 0, 1)
integrate(\(x) gamma(x) / gamma(1/x), 0,1)
integrate(\(x) gamma(x) / gamma(1 - x) - 1/x, 0, 1)
integrate(\(x) gamma(log(x+1)) - 1/x, 0, 1)
integrate(\(x) gamma(log(x+1)) - 1/x, 0, exp(1)-1)
integrate(\(x) gamma(x) / pracma::zeta(x+1), 0, 1)
integrate(\(x) gamma(x) / pracma::zeta(1-x), 0, 1)
integrate(\(x) gamma(gamma(x) - 1/x) + 1 / ((1-x)*(1-Euler)), 0, 1)
integrate(\(x) (gamma(x) - 1)*(gamma(1-x) - 1), 0, 1)

integrate(\(x) - gamma(x) / (1-x) + 1/(x*(1-x)), 0, 1)
# =>
integrate(\(x) gamma(x) - 1/x + 1/(x+1), -1, 0)

# Varia:
integrate(\(x) gamma(x)/log(gamma(x)) + 1/(x*log(x)) - (1-Euler)/Euler/(1-x), 0, 1, rel.tol=1E-5)


###
integrate(\(x) x * log(gamma(x)*gamma(1-x)), 0, 1)
log(2*pi)/2

###
integrate(\(x) gamma(x) * gamma(1-x) - 1/x - 1/(1-x), 0, 1)
- 2 * log(pi/2)

integrate(\(x) x * log(sin(x)), 0, pi)
- pi^2 * log(2)/2
integrate(\(x) x * log(sin(pi*x)), 0, 1)
- log(2)/2


