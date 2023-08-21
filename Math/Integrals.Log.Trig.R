

### Integrals: Log(Trig)

# - this file covers primarily
#   integrals of type: Log( Trig );


### Helper

Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)
# Catalan = (psi(1, 1/4) - psi(1, 3/4)) / 4^2


###################
###################

### I( log(sin(x)) ) on various intervals

# see e.g.:
# Maths 505: see 2nd & 3rd intermediate integrals;
# A MONSTER INTEGRAL!!! int from 0 to infty arctan(x)/(x(x+1)(x^2+1))
# https://www.youtube.com/watch?v=u-FKjn_83l8

# TODO: move to Trig;
integrate(\(x) x / (tan(x) + 1), 0, pi/2)
pi^2/16 + pi*log(2)/8 - Catalan/2


###
integrate(function(x) log(sin(x)), 0, pi/2)
integrate(function(x) log(cos(x)), 0, pi/2)
- pi/2*log(2)
#
integrate(\(x) log(1 + cos(x)), 0, pi/2)
integrate(\(x) log(1 + sin(x)), 0, pi/2)
- pi/2*log(2) + 2*Catalan

###
integrate(function(x) log(sin(x)), 0, pi/4)
- pi/4*log(2) - Catalan/2

###
integrate(function(x) log(cos(x)), 0, pi/4)
- pi/4*log(2) + Catalan/2

###
integrate(function(x) log(cos(x) + sin(x)), 0, pi/2)
- pi/4*log(2) + Catalan
#
integrate(function(x) log(cos(x) - sin(x)), 0, pi/4)
- pi/8*log(2) - Catalan/2
#
integrate(function(x) log(cos(x) + sin(x)), 0, pi/4)
- pi/8*log(2) + Catalan/2


### I( log(sin(x)) ) on [0, pi/8]
# TODO: C8;
C8 = integrate(\(x) log(x) / (x^2 + 1), 0, tan(pi/8))$value
#
integrate(function(x) log(sin(x)), 0, pi/8)
- pi/8*log(2) - Catalan/8 + C8/2
#
integrate(function(x) log(cos(x)), 0, pi/8)
- pi/8*log(2) - Catalan/8 - C8/2


### on [0, pi/3]

###
integrate(\(x) log(cos(x)), 0, pi/3)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 1/3) - pracma::psi(1, 2/3))

###
integrate(\(x) log(sin(x)), 0, pi/6)
- pi*log(2)/6 - sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 1/3) - pracma::psi(1, 2/3))

### I( log(sin(x)) )
integrate(\(x) log(sin(x)), 0, pi/3)
- pi*log(2)/3 + sqrt(3)/(8*6^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1, 7/12) +
	- 5 * pracma::psi(1, 1/6) + 5 * pracma::psi(1, 5/6) +
	- 3 * pracma::psi(1, 2/6) + 3 * pracma::psi(1, 4/6));


# Derivation:
# for technique see from 12:40 in the Ref. Maths 505;
id = seq(0, 40000)
- pi*log(2)/3 + sqrt(3)/4 * sum(1/(6*id+1)^2, -1/(6*id+5)^2, 1/(6*id+2)^2, - 1/(6*id+4)^2)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 2/6) - pracma::psi(1, 4/6))


### on [0, pi/5]
integrate(\(x) log(cos(x)), 0, pi/5)
sn = sin(2*pi/5 * c(1,2));
- pi*log(2)/5 +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100);

# Note: also alternating signs;
- pi*log(2)/5 + 1/2 * sum(
	sn[1]/(10*id+1)^2, - sn[1]/(10*id+9)^2, - sn[2]/(10*id+2)^2, sn[2]/(10*id+8)^2,
	- sn[2]/(10*id+3)^2, sn[2]/(10*id+7)^2, sn[1]/(10*id+4)^2, - sn[1]/(10*id+6)^2)


### on [0, pi * 2/5]
integrate(\(x) log(cos(x)), 0, 2*pi/5)
sn = sin(2*pi/5 * c(2,4));
- pi*log(2) * 2/5 +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100);

# Derivation:
- pi*log(2) * 2/5 + 1/2 * sum(
	sn[1]/(10*id+1)^2, - sn[1]/(10*id+9)^2, - sn[2]/(10*id+2)^2, sn[2]/(10*id+8)^2,
	- sn[2]/(10*id+3)^2, sn[2]/(10*id+7)^2, sn[1]/(10*id+4)^2, - sn[1]/(10*id+6)^2)


### on [0, pi/6]
integrate(\(x) log(cos(x)), 0, pi/6)
- pi*log(2)/6 + sqrt(3)/(4*12^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) +
	- pracma::psi(1, 1/6) + pracma::psi(1, 5/6) +
	+ pracma::psi(1, 2/6) - pracma::psi(1, 4/6));

#
id = seq(0, 40000)
sn = sin(2*pi/6 * c(1,2,3));
- pi*log(2)/6 + 1/2 * sum(
	sn[1]/(12*id+1)^2, - sn[1]/(12*id+11)^2, - sn[2]/(12*id+2)^2, sn[2]/(12*id+10)^2,
	sn[1]/(12*id+4)^2, - sn[1]/(12*id+8)^2, - sn[2]/(12*id+5)^2, + sn[2]/(12*id+7)^2)


### Ap / Cp Constants
# - useful for: I( log(gamma(x+k)) )
Ap = function(p, iter=10000) {
	n = iter; id = seq(n) + 1/p;
	sum(id * log(id)) - (n^2 + (2/p + 1)*n + 1/p^2 + 1/p + 1/6)*log(n + 1/p)/2 +
		+ n^2/4 + n/(2*p);
}
#
Ap(6) - Ap(-6)
log(6)/6 - sqrt(3)/(4*36*pi) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 1/3) - pracma::psi(1, 2/3))


### Varia: Clausen Function
# see e.g.
# 1. Michael Penn: When "normal" trig functions aren't enough -- the Clausen function.
#    https://www.youtube.com/watch?v=5kN4oH8W1r8


############

### I( x^p * log(sin(x)) )

###
integrate(function(x) x * log(sin(x)), 0, pi/2)
- 1/8*pi^2*log(2) + 7/16 * pracma::zeta(3)

###
integrate(function(x) x^2 * log(sin(x)), 0, pi/2)
- 1/(3*8)*pi^3*log(2) + 3/16 * pi * pracma::zeta(3)


### Helper
x  = pi/7
id = seq(10000);
log(sin(x))
- sum( cos(2*id*x) / id ) - log(2);


# Varia:
integrate(\(x) - log(x) / (x^2 + x + 1), 0, 1)
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9
# TODO: psi(1, x) - psi(1, 1-x) ???


#################

### I( log(tan(x)) )
# Michael Penn: Can you guess the trick for this integral?
# https://www.youtube.com/watch?v=8R0MiRYmjbk

# Intermediary:
integrate(function(x) log(tan(x)), 0, pi/4)
# ==
- Catalan


### I( log(cos x)^2 )
# Michael Penn: Integral of (ln(cos x))^2
# https://www.youtube.com/watch?v=ikyVHEHmgP8

integrate(\(x) log(cos(x))^2, 0, pi/2)
pi^3/24 + pi/2 * log(2)^2


#########################

### I( log(cos(x)) / tan(x) )
# Maths 505: A cool log trig integral
# https://www.youtube.com/watch?v=iPicVw4lOg0


### I( log(cos(x)) / tan(x) )
integrate(\(x) log(cos(x)) / tan(x), 0, pi/2)
integrate(\(x) 1/2 * log(1 - x^2) / x, 0, 1)
- pi^2/24

### I( log(cos(x)) / sin(x) )
integrate(\(x) log(cos(x)) / sin(x), 0, pi/2)
integrate(\(x) log(x) / (1 - x^2), 0, 1)
- pi^2/8

### I( log(cos(x)) / cos(x) )
integrate(\(x) log(cos(x)) / cos(x) - log(pi/2 - x)/(pi/2 - x), 0, pi/2)
integrate(\(x) pi/2*log(cos(pi/2*x)) / cos(pi/2*x) - (log(pi/2) + log(1-x))/(1 - x), 0, 1)
# TODO: ???

# equivalent:
integrate(\(x) 1/2 * log(1 - x^2) / (1 - x^2) - 1/4*(log(1-x) + log(2)) / (1-x), 0, 1)
3*log(2)^2/8 - pi^2/24

# Lim: x -> 1
x = 1 - 1E-6;
log(1 - x^2) / (1 + x) - 1/2*(log(1-x) + log(2))
log1p(- x^2) / (1 + x) - 1/2*(log1p(-x) + log(2))
0;


#########################
#########################

### Derived Integrals

### I( log(x) / sqrt(2*b*x - x^2) )
# Maths 505: A nice integral from 1886
# https://www.youtube.com/watch?v=CbAcyU8gFPw

###
integrate(\(x) log(x) / sqrt(4*x - x^2), 0, 2)
- 2*Catalan

###
integrate(\(x) log(x) / sqrt(6*x - x^2), 0, 3)
pi/2 *  log(3/2) - 2*Catalan

###
b = sqrt(5)
integrate(\(x) log(x) / sqrt(2*b*x - x^2), 0, b)
pi/2 *  log(b/2) - 2*Catalan

# TODO: log(sin)-interval [0, pi/3];

