


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


