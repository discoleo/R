

### Integrals: Log(Trig)


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
integrate(function(x) log(sin(x)), 0, pi/8)
- pi/8*log(2) - Catalan/8 + C8/2
#
integrate(function(x) log(cos(x)), 0, pi/8)
- pi/8*log(2) - Catalan/8 - C8/2


### on [0, pi/3]

integrate(\(x) log(cos(x)), 0, pi/3)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 1/3) - pracma::psi(1, 2/3))

# Derivation:
id = seq(0, 40000)
- pi*log(2)/3 + sqrt(3)/4 * sum(1/(6*id+1)^2, -1/(6*id+5)^2, 1/(6*id+2)^2, - 1/(6*id+4)^2)
- pi*log(2)/3 + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) + pracma::psi(1, 2/6) - pracma::psi(1, 4/6))

# Varia:
(pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / 9
integrate(\(x) - log(x) / (x^2 + x + 1), 0, 1)
# TODO


#################

### Michael Penn: Can you guess the trick for this integral?
# https://www.youtube.com/watch?v=8R0MiRYmjbk

# Intermediary:
integrate(function(x) log(tan(x)), 0, pi/4)
# ==
- Catalan
