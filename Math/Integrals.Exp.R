

### Helper

Euler = 0.57721566490153286060651209008240243079;
constEuler = Euler;


### Examples:
# I( x^p / (exp(k*x) + 1) ) on [0, Inf]
# I( x^p / (exp(k*x) - 1) ) on [0, Inf]


### History

# - [refactor] moved:
#    I( atan(x) / (exp(x) - 1) ) and
#    I( atan(x) / (exp(x) + 1) )
#    to new file: Integrals.Exp.Advanced.R;
# - [refactor] moved:
#    I( sin(x) / (exp(x) - 1) )
#    I( x * sin(k*x) / (exp(x) - 1) )
#    I( cos(x) / (exp(1/x) - 1) )
#    I( cos(x) / (exp(1/x) + 1) )
#    to file: Integrals.Exp.Trig.R;


####################
####################

### Basics

### on [0, Inf]

### I( exp(-x^n) / x )
integrate(\(x) exp(- x^2)/x - exp(-x)/x, 0, Inf)
Euler / 2

###
integrate(\(x) exp(- x^3)/x - exp(-x)/x, 0, Inf)
Euler * 2/3

###
n = sqrt(5)
integrate(\(x) exp(- x^n)/x - exp(-x)/x, 0, Inf)
Euler * (1 - 1/n)


### Log:
integrate(\(x) log(x) * exp(- x), 0, Inf)
- Euler


#############

### on [0, 1]

### I( x^p * exp(x) ) on [0, 1]
# using an Erf1-generalization:
p = sqrt(3)
integrate(\(x) x^p * exp(x), 0, 1, rel.tol = 1E-8)
exp(1) - integrate(\(x) exp(x^(1/p)), 0, 1, rel.tol = 1E-8)$value


###
gamma.p.part = function(up, n=2, rel.tol=1E-8) {
	integrate(\(x) exp(x^n), 0, up, rel.tol = rel.tol)$value;
}
gamma.part = function(up, n=2, rel.tol=1E-8) {
	integrate(\(x) exp(-x^n), 0, up, rel.tol = rel.tol)$value;
}

#
integrate(\(x) x^(1/2) * exp(x^2), 0, 1, rel.tol = 1E-8)
integrate(\(x) 2*x^2 * exp(x^4), 0, 1, rel.tol = 1E-8)
# direct:
integrate(\(x) 2/3 * exp(x^(2/(1/2 + 1))), 0, 1, rel.tol = 1E-8)

# Complicated way:
# Note: up = x^2;
gamma.p.part(1) - integrate(\(x) sapply(x^2, gamma.p.part), 0, 1)$value;
# Note:
gamma.p.part(1, 2) + 2/3 * gamma.p.part(1, 2/3) # == exp(1)

###
p = sqrt(3)
integrate(\(x) x^p * exp(x^2), 0, 1, rel.tol = 1E-8)
integrate(\(x) 1/(p+1) * exp(x^(2/(p+1))), 0, 1, rel.tol = 1E-8)
gamma.p.part(1, n = 2/(p+1)) / (p+1)

#
p = sqrt(3); n = sqrt(5);
integrate(\(x) x^p * exp(x^n), 0, 1, rel.tol = 1E-8)
gamma.p.part(1, n = n/(p+1)) / (p+1)


####################
####################

###
integrate(\(x) exp(-x)/x^2 - 1/x^2 + exp(-x)/x, 0, 1)
- exp(-1)

###
integrate(\(x) exp(-x)/x^2 - 1/x^2 + 1/x, 0, 1)
pracma::expint(1) - exp(-1) + Euler

###
integrate(\(x) exp(-x)/x - 1/x, 0, 1)
- (pracma::expint(1) + Euler)


####################
####################

### I( (1 - exp(-x) - exp(-1/x)) / x )
# Michael Penn: if you gamma be my constant, you gamma solve this integral.
# https://www.youtube.com/watch?v=heIrFaP05Lk

integrate(\(x) (1 - exp(-x) - exp(-1/x)) / x, 0, 1)
Euler


###
integrate(\(x) exp(-1/x), 0, 1)
exp(-1) - pracma::expint(1)


####################
####################

### I( x^p * (exp(-a*x) - exp(-b*x)) )
# Dr Peyam: A multivariable integral
# https://www.youtube.com/watch?v=NLNDdfbQBWU

### p = -1
a = c(3, 2)
integrate(\(x) (exp(-a[1]*x) - exp(-a[2]*x)) / x, 0, Inf)
log(a[2] / a[1])


### Generalization:
# I( I( x^(p+1) * exp(-x*y) ) )
# x on [0, Inf] & y on [a, b];

# p > -1;
a = sqrt(c(3, 2))
integrate(\(x) x^(-1/2) * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(1/2) * (a[1]^(-1/2) - a[2]^(-1/2))

### I( x^p * (exp(-a1*x) - exp(-a2*x)) )
p = -1/3; # p > -2 (strictly) && p != - 1 !!!
a = sqrt(c(5, 2))
integrate(\(x) x^p * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(p + 1) * (a[1]^(-p-1) - a[2]^(-p-1))

###
p = sqrt(2) - sqrt(3); # p > -2 (strictly) !!!
a = sqrt(c(5, 2))
integrate(\(x) x^p * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(p + 1) * (a[1]^(-p-1) - a[2]^(-p-1))

###
p = sqrt(7); # p > -2 (strictly) !!!
a = sqrt(c(5, 2))
integrate(\(x) x^p * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(p + 1) * (a[1]^(-p-1) - a[2]^(-p-1))

# TODO: Lim: p -> -1


#########################
#########################

### I( exp(- x^2 - a^2/x^2) )
# Maths 505: A classic integral solved using Feynman's technique
# https://www.youtube.com/watch?v=dvRarIKUMSE


a = 1/3
integrate(\(x) exp(- x^2 - a^2/x^2), 0, Inf)
gamma(1/2) * exp(-2*a) / 2


### [hidden x^2]
a = 1/3
integrate(\(x) exp(- x^3 - a^3/x^3) * x^(1/2), 0, Inf)
gamma(1/2) * exp(-2*a^(3/2)) * a


#########################
#########################

### I( x^p / (exp(k*x) + 1) )
# Michael Penn: WHY are we finding pi HERE?
# https://www.youtube.com/watch?v=amL2uBK2buU

p = 1/3
integrate(\(x) x^p / (exp(x) + 1), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p
gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p)

###
p = sqrt(5)
integrate(\(x) x^p / (exp(x) + 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p)

###
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (exp(k*x) + 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) / k^(p + 1)
#
integrate(\(x) x^p / (exp(k*x) - 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1) / k^(p + 1)

###
# - bringing pi back into the equation:
p = sqrt(5)
integrate(\(x) x^p / (1 + pi^x), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) / log(pi)^(p + 1)


### I( x^p / (exp(k*x) - 1) )
# Dr Peyam: A gamma zeta integral
# https://www.youtube.com/watch?v=hqo_-AX4IXg

p = sqrt(5)
integrate(\(x) x^p / (exp(x) - 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1)

#
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (exp(k*x) - 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1) / k^(p+1)


#################

### Higher Powers
### I( x^p / (exp(k*x) + 1)^r )

###
k = sqrt(3)
integrate(\(x) 1 / (exp(k*x) + 1), 0, 100)
log(2)/k

###
k = sqrt(3)
integrate(\(x) 1 / (exp(k*x) + 1)^2, 0, 100)
(log(2) - 1/2)/k

###
k = sqrt(3)
integrate(\(x) x * exp(k*x) / (exp(k*x) + 1)^2, 0, 100)
log(2) / k^2


### I( x^p / (exp(k*x) + 1)^2 )
k = sqrt(3); p = sqrt(5);
integrate(\(x) x^p / (exp(k*x) + 1)^2, 0, 100)
gamma(p + 1) *
	( pracma::zeta(p + 1) * (2^p - 1) +
	- pracma::zeta(p) * (2^p - 2) ) / (2^p * k^(p+1));


### I( x^p * exp(k*x) / (exp(k*x) + 1)^3 )
k = 1/sqrt(3); p = sqrt(5);
# Up = Inf; numerical issues;
integrate(\(x) x^p * exp(k*x) / (exp(k*x) + 1)^3, 0, 100)
gamma(p + 1) *
	( pracma::zeta(p) * (2^(p-1) - 1) +
	- pracma::zeta(p - 1) * (2^(p-1) - 2) ) / (2^p * k^(p+1));


### I( x^p / (exp(k*x) + 1)^3 )
# =>
k = 1/sqrt(3); p = sqrt(5);
# Up = Inf; numerical issues;
integrate(\(x) x^p / (exp(k*x) + 1)^3, 0, 100)
gamma(p + 1) *
	( pracma::zeta(p + 1) * (2^p - 1) - 3 * pracma::zeta(p) * (2^(p-1) - 1) +
	+ pracma::zeta(p - 1) * (2^(p-1) - 2) ) / (2^p * k^(p+1));


#########################

### I( 1 / (exp(x^n) + 1) )
# Maths 505: A RIDICULOUSLY AWESOME integral from the quantum realm
# https://www.youtube.com/watch?v=oZWsDyN5ssI

n = sqrt(5)
integrate(\(x) 1 / (exp(x^n) + 1), 0, Inf)
gamma(1/n) * (1 - 2^(1-1/n)) * pracma::zeta(1/n) / n;

### Generalizations:

### I( x^p / (exp(x^n) + 1) )
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p / (exp(x^n) + 1), 0, Inf)
gamma((p+1)/n) * (1 - 2^(1-(p+1)/n)) * pracma::zeta((p+1)/n) / n;

### I( x^p / (exp(x^n) - 1) )
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p / (exp(x^n) - 1), 0, Inf)
gamma((p+1)/n) * pracma::zeta((p+1)/n) / n;


### Transforms

### I( x^p * log((exp(x) - 1) / (exp(x) + 1)) )
p = sqrt(5);
# numerical issues: upper -> Inf;
integrate(\(x) x^p * log((exp(x) - 1) / (exp(x) + 1)), 0, 200, rel.tol=1E-8)
- gamma(p + 1) * pracma::zeta(p + 2) * (2 - 1/2^(p+1))


### I( x^p * log((exp(x^n) - 1) / (exp(x^n) + 1)) )
p = sqrt(5); n = sqrt(3);
# numerical issues: upper -> Inf;
integrate(\(x) x^p * log((exp(x^n) - 1) / (exp(x^n) + 1)), 0, 40, rel.tol=1E-8)
- gamma((p + 1)/n) * pracma::zeta((p + 1)/n + 1) * (2 - 1/2^((p + 1)/n)) / n

