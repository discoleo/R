

constEuler = 0.57721566490153286060651209008240243079;

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

### I( x^p * exp(x) ) on [0, 1]
# using an Erf1-generalization:
p = sqrt(3)
integrate(\(x) x^p * exp(x), 0, 1, rel.tol = 1E-8)
exp(1) - integrate(\(x) exp(x^(1/p)), 0, 1, rel.tol = 1E-8)$value


###
gamma.part = function(up, n=2, rel.tol=1E-8) {
	integrate(\(x) exp(x^n), 0, up, rel.tol = rel.tol)$value;
}

#
integrate(\(x) x^(1/2) * exp(x^2), 0, 1, rel.tol = 1E-8)
integrate(\(x) 2*x^2 * exp(x^4), 0, 1, rel.tol = 1E-8)
# direct:
integrate(\(x) 2/3 * exp(x^(2/(1/2 + 1))), 0, 1, rel.tol = 1E-8)

# Complicated:
# Note: up = x^2;
gamma.part(1) - integrate(\(x) sapply(x^2, gamma.part), 0, 1)$value;
# Note:
gamma.part(1, 2) + 2/3 * gamma.part(1, 2/3) # == exp(1)

###
p = sqrt(3)
integrate(\(x) x^p * exp(x^2), 0, 1, rel.tol = 1E-8)
integrate(\(x) 1/(p+1) * exp(x^(2/(p+1))), 0, 1, rel.tol = 1E-8)
gamma.part(1, n = 2/(p+1)) / (p+1)

#
p = sqrt(3); n = sqrt(5);
integrate(\(x) x^p * exp(x^n), 0, 1, rel.tol = 1E-8)
gamma.part(1, n = n/(p+1)) / (p+1)


####################
####################

### I( (1 - exp(-x) - exp(-1/x)) / x )
# Michael Penn: if you gamma be my constant, you gamma solve this integral.
# https://www.youtube.com/watch?v=heIrFaP05Lk

integrate(\(x) (1 - exp(-x) - exp(-1/x)) / x, 0, 1)
constEuler


#########################
#########################

### I( x^p / (exp(x) + 1) )
# Michael Penn: WHY are we finding pi HERE?
# https://www.youtube.com/watch?v=amL2uBK2buU

p = 1/3
integrate(\(x) x^p / (1 + exp(x)), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p
gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p)

###
p = sqrt(5)
integrate(\(x) x^p / (1 + exp(x)), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p

###
# - bringing pi back into the equation:
p = sqrt(5)
integrate(\(x) x^p / (1 + pi^x), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p / log(pi)^(p + 1)


### I( x^p / (exp(x) - 1) )
# Dr Peyam: A gamma zeta integral
# https://www.youtube.com/watch?v=hqo_-AX4IXg

p = sqrt(5)
integrate(\(x) x^p / (exp(x) - 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * pracma::zeta(p + 1)


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


#####################
#####################

### I( exp(-x^n) / (x^n + 1) ) on [0, Inf]
# Maths 505: A crazy yet perfect integral
# https://www.youtube.com/watch?v=regnh-HOR0c
# [the sub-integral]

###
n = sqrt(3)
integrate(\(x) exp(-x^n) / (x^n + 1), 0, Inf)
int = integrate(\(x) x^(-1/n) * exp(-x), 0, 1)$value
exp(1) * (pi/sin(pi/n) - int*gamma(1/n))/n

### [stable]
n = sqrt(3) - sqrt(2)
integrate(\(x) exp(-x^n) / (x^n + 1), 0, Inf)
int = integrate(\(x) n*x^(n-2) * exp(-x^n), 1, Inf)$value
exp(1) * int * gamma(1/n) / n

# only: n > 1
# int = integrate(\(x) x^(-1/n) * exp(-x), 0, 1)$value
int = integrate(\(x) n*x^(n-2) * exp(-x^n), 0, 1)$value
exp(1) * (pi/sin(pi/n) - int*gamma(1/n)) / n

