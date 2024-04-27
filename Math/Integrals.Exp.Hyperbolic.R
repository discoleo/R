



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


### Gen: I( 1 / cosh(x)^s )
s = sqrt(3)/sqrt(5);
integrate(\(x) 1 / cosh(x)^s, 0, Inf)
gamma(s/2)^2 / gamma(s) * 2^(s-2)


### Full Gen: I( cosh(p*x) / cosh(n*x)^k )
# see file: Integrals.ComplexAnalysis.R;
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
integrate(\(x) cosh(p*x) / cosh(n*x)^k, lower=0, upper=100) # BUG: upper = Inf
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) / gamma(k) * 2^(k - 2) / n

### SINH-Variants

### I( x * sinh(p*x) / cosh(n*x)^(1/k) )
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
# BUG: upper = Inf
integrate(function(x) x * sinh(p*x) / cosh(n*x)^k, lower=0, upper=100)
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) *
	(digamma(k/2 + p/(2*n)) - digamma(k/2 - p/(2*n))) / gamma(k) * 2^(k - 2) / (2*n^2)


### Simple:

### I( 1 / sinh(x) )
# BUG: upper = Inf;
integrate(\(x) 1 / sinh(x) - exp(-x)/x, 0, 20)
log(2) + Euler

###
b = 3
# BUG: upper = Inf;
integrate(\(x) 1 / sinh(b*x) - exp(-b*x)/(b*x), 0, 10)
(log(2) + Euler) / b


### I( sinh(b*x) * cos(w*x) / sinh(x) )
b = sqrt(3); w = sqrt(2);
# BUG: upper = Inf;
integrate(\(x) sinh(b*x) * cos(w*x) / sinh(pi*x), 0, 40)
1/2 * sin(b) / (cosh(w) + cos(b))

###
b = sqrt(3) / pi; w = sqrt(2);
# BUG: upper = Inf;
integrate(\(x) sinh(b*x) * cos(w*x) / sinh(x), 0, 40)
pi/2 * sin(pi*b) / (cosh(pi*w) + cos(pi*b))


### x-Variants

### I( x * sinh(b*x) * sin(w*x) / sinh(x) )
b = 1 / sqrt(3); w = sqrt(2);
# BUG: upper = Inf;
integrate(\(x) x * sinh(b*x) * sin(w*x) / sinh(x), 0, 40)
pi^2/2 * sinh(pi*w) * sin(pi*b) / (cosh(pi*w) + cos(pi*b))^2


####################
### Log-Variants ###

### I( cosh(p*x) * log(cosh(n*x)) / cosh(n*x)^k ) on [0, Inf]
n = sqrt(19); k = 1/sqrt(3); p = sqrt(2);
# BUG: upper = Inf
integrate(function(x) cosh(p*x) * log(cosh(n*x)) / cosh(n*x)^k, lower=0, upper=100)
gamma(k/2 - p/(2*n)) * gamma(k/2 + p/(2*n)) / gamma(k) * 2^(k - 3) / n *
	(2*digamma(k) - digamma(k/2 - p/(2*n)) - digamma(k/2 + p/(2*n)) - 2*log(2));


### I( log(x) / cosh(x) )
# see Vardi Integral
integrate(\(x) log(x) / cosh(x), 0, Inf, rel.tol=1E-8)
- log(gamma(1/4) / (gamma(3/4) * sqrt(2*pi))) * pi


### I( log(x^2 + b^2) / cosh(x) )
# Generalization of Vardi Integral
# see additional details in Integrals.Exp.Log.R;
b = sqrt(5);
integrate(\(x) log(x^2 + b^2) / cosh(x), 0, Inf)
2*pi * log(gamma(b/(2*pi) + 3/4) / gamma(b/(2*pi) + 1/4) * sqrt(2*pi))

