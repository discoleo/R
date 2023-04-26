

### I( x^p / (1 + exp(x)) )
# Michael Penn: WHY are we finding pi HERE?
# https://www.youtube.com/watch?v=amL2uBK2buU

p = 1/3
integrate(\(x) x^p / (1 + exp(x)), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p

###
p = sqrt(5)
integrate(\(x) x^p / (1 + exp(x)), 0, Inf, rel.tol=1E-8)
gamma(p) * pracma::zeta(p + 1) * (1 - 1/2^p) * p
