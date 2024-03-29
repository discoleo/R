



### I( x^2 / cosh(x^2)^2 )
# Maths 505: This EPIC integral is the best thing you'll see today
# https://www.youtube.com/watch?v=aQEJevhQjTQ

integrate(\(x) x^2 / cosh(x^2)^2, 0, Inf)
sqrt(pi/8)*(1 - sqrt(2)) * pracma::zeta(1/2)


###
p = sqrt(5)
integrate(\(x) x^p / cosh(x^2)^2, 0, Inf)
gamma(p/2 + 1/2) * pracma::zeta(p/2 - 1/2) * (1 - 2^(3/2 - p/2)) / 2^(p/2 - 1/2)


###
p = sqrt(5)
n = sqrt(3)
integrate(\(x) x^p / cosh(x^n)^2, 0, Inf)
gamma((p + 1)/n) * pracma::zeta((p + 1)/n - 1) * (1 - 2^(2 - (p + 1)/n)) / 2^((p + 1)/n - 1) * 2/n

# TODO: even more generalizations


######################
######################

### I( 1 / cosh(x)^s )
# Maths 505: An epic exponential function integral
# https://www.youtube.com/watch?v=rKkTGwIpODs

###
integrate(\(x) 1 / (2*cosh(x))^(1/2), 0, Inf)
gamma(1/4)^2 / gamma(1/2) /4


###
integrate(\(x) 1 / (2*cosh(x))^(1/3), 0, Inf)
gamma(1/6)^2 / gamma(1/3) / 4


### Gen:
s = sqrt(3)/sqrt(5);
integrate(\(x) 1 / (2*cosh(x))^s, 0, Inf)
gamma(s/2)^2 / gamma(s) / 4


### Full Gen:
# see file: Integrals.ComplexAnalysis.R;
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
integrate(function(x) cosh(p*x) / cosh(n*x)^k, lower=0, upper=100) # BUG: upper = Inf
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) / gamma(k) * 2^(k - 2) / n

### I( x * sinh(p*x) / cosh(n*x)^(1/k) )
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
# BUG: upper = Inf
integrate(function(x) x * sinh(p*x) / cosh(n*x)^k, lower=0, upper=100)
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) *
	(digamma(k/2 + p/(2*n)) - digamma(k/2 - p/(2*n))) / gamma(k) * 2^(k - 2) / (2*n^2)


### I( cosh(p*x) * log(cosh(n*x)) / cosh(n*x)^k ) on [0, Inf]
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
# BUG: upper = Inf
integrate(function(x) cosh(p*x) * log(cosh(n*x)) / cosh(n*x)^k, lower=0, upper=100)
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) / gamma(k) * 2^(k - 3) / n *
	(2*digamma(k) - digamma(k/2 - p/(2*n)) - digamma(k/2 + p/(2*n)) - 2*log(2));


### I( log(x) / cosh(x) )
# see Vardi Integral
pracma::integral(\(x) log(x) / cosh(x), 0, Inf)
- log(gamma(1/4) / (gamma(3/4) * gamma(1/2) * sqrt(2))) * pi
