

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


### I( x^p * log((exp(x) - 1) / (exp(x) + 1)) )
p = sqrt(5);
# numerical issues: upper -> Inf;
integrate(\(x) x^p * log((exp(x) - 1) / (exp(x) + 1)), 0, 200, rel.tol=1E-8)
- gamma(p + 1) * pracma::zeta(p + 2) * (2 - 1/2^(p+1))


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


