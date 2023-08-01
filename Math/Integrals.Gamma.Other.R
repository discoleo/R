


dzeta2   = -0.937548254316;
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


########################
########################

### I( log(gamma(x + k)) )
# Flammable Maths:
# WOW! The Most AMAZING Way of Solving an Integral Ever! Deriving Raabe's Integral Formula!
# https://www.youtube.com/watch?v=APRr6wmLyuk

###
k = sqrt(5)
integrate(\(x) log(gamma(x)), k, k+1)
integrate(\(x) log(gamma(x + k)), 0, 1)
log(2*pi)/2 + k*log(k) - k


### I( log(gamma(2*x + k)) )
k = sqrt(5)
integrate(\(x) log(gamma(2*x + k)), 0, 1)
((k+1)*log(k+1) + k*log(k)) / 2 - k + log(2*pi)/2 - 1/2

### I( log(gamma(3*x + k)) )
k = sqrt(5)
integrate(\(x) log(gamma(3*x + k)), 0, 1)
(k*log(k) + (k+1)*log(k+1) + (k+2)*log(k+2)) / 3 - k + log(2*pi)/2 - 1


###############

### I( x * log(gamma(x)) )
integrate(\(x) x * log(gamma(x)), 0, 1)
log(2*pi)/4 - log(A);


### TODO:
integrate(\(x) log(1-x) * log(gamma(x)), 0, 1)
integrate(\(x) log(x) * log(gamma(x)), 0, 1)

integrate(\(x) log(gamma(x)) / (x + 1), 0, 1)
integrate(\(x) log(gamma(x)) / (x^2 + 1), 0, 1)
integrate(\(x) log(gamma(x)) / log(x), 0, 1)
integrate(\(x) log(gamma(x)) / gamma(x), 0, 1)
integrate(\(x) log(gamma(x)) * log(gamma(1-x)), 0, 1)

integrate(\(x) 1/gamma(x), 0, 1)
integrate(\(x) x / gamma(x), 0, 1)
integrate(\(x) log(x) / gamma(x), 0, 1)
integrate(\(x) sin(x) / gamma(x), 0, 1)

integrate(\(x) gamma(x) * x, 0, 1)
integrate(\(x) gamma(x) * sin(x), 0, 1)
integrate(\(x) gamma(x) * sin(pi*x), 0, 1)
integrate(\(x) gamma(x) * (exp(x) - 1), 0, 1)
integrate(\(x) gamma(x) * log(x + 1), 0, 1)
integrate(\(x) gamma(x) * log(x*exp(x) + 1), 0, 1)
integrate(\(x) gamma(x) * atan(x), 0, 1)
integrate(\(x) atan(gamma(x)), 0, 1)
integrate(\(x) gamma(x) / gamma(1/x), 0,1)
integrate(\(x) gamma(x) / pracma::zeta(x+1), 0, 1)
integrate(\(x) gamma(x) / pracma::zeta(1-x), 0, 1)


###
integrate(\(x) x * log(gamma(x)*gamma(1-x)), 0, 1)
log(2*pi)/2

integrate(\(x) x * log(sin(x)), 0, pi)
- pi^2 * log(2)/2
integrate(\(x) x * log(sin(pi*x)), 0, 1)
- log(2)/2


