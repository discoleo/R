#################
##
## Integrals: Tan
## Powers


# For Polynomial Fractions:
# see file Integrals.Fractions.Unity.R;

##################

up = pi/5
integrate(\(x) tan(x), 0, up)
- log(sin(up))

###
up = pi/5
integrate(\(x) tan(x)^2, 0, up)
tan(up) - up


###
up = pi/5
integrate(\(x) sqrt(tan(x)), 0, up)
integrate(\(x) 2*x^2 / (x^4 + 1), 0, sqrt(tan(up)))


###
up = pi/5
integrate(\(x) tan(x)^(1/3), 0, up)
integrate(\(x) 3*x^3 / (x^6 + 1), 0, tan(up)^(1/3))


### on [0, pi/4]

n = 1/sqrt(3)
integrate(\(x) tan(x)^n, 0, pi/4)
integrate(\(x) 1/n * x^(1/n) / (x^(2/n) + 1), 0, 1)
(digamma(((n+1)/2 + 1)/2) - digamma((n+1)/4)) / 4
(digamma((n+3)/4) - digamma((n+1)/4)) / 4


### Derived Variants


n = 1/sqrt(3)
integrate(\(x) tan(x)^n * log(tan(x)), 0, pi/4)
(pracma::psi(1, (n+3)/4) - pracma::psi(1, (n+1)/4)) / 16

###
n = 1/sqrt(3)
integrate(\(x) tan(x)^n * log(tan(x))^2, 0, pi/4)
(pracma::psi(2, (n+3)/4) - pracma::psi(2, (n+1)/4)) / 64

