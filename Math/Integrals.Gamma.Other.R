

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

