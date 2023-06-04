

constEuler = 0.57721566490153286060651209008240243079;

### I( (1 - exp(-x) - exp(-1/x)) / x )
# Michael Penn: if you gamma be my constant, you gamma solve this integral.
# https://www.youtube.com/watch?v=heIrFaP05Lk

integrate(\(x) (1 - exp(-x) - exp(-1/x)) / x, 0, 1)
constEuler


#########################
#########################

### I( cos(x) / (exp(1/x) + 1) )
# Michael Penn: is this integration trick TOO POWERFUL?
# https://www.youtube.com/watch?v=PX2QXILRgsc

###
lim = 1
integrate(\(x) cos(x) / (exp(1/x) + 1), -lim, lim)
sin(lim)

### Various Variants:
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) + 1), -lim, lim)
lim

###
lim = 1
integrate(\(x) cos(x)^2 / (exp(1/x) + 1), -lim, lim)
sin(2*lim) / 4 + lim/2

###
lim = sqrt(pi)
integrate(\(x) sin(x)^2 / (exp(1/x) + 1), -lim, lim)
lim/2 - sin(2*lim) / 4


### Variants: exp() - 1

###
lim = sqrt(pi)
integrate(\(x) 1 / (exp(1/x) - 1), -lim, lim)
- lim

###
lim = sqrt(pi)
integrate(\(x) cos(x) / (exp(1/x) - 1), -lim, lim)
- sin(lim)


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


#####################
#####################

### I( sin(x) / (exp(x) - 1))
# Maths 505: One of the coolest integrals on YouTube!
# https://www.youtube.com/watch?v=I7a-6VwKhgo

###
integrate(\(x) sin(x) / (exp(x) - 1), 0, Inf)
pi/2 / tanh(pi) - 1/2
(digamma((1/2+1)/2) - digamma(1/4))/2 - 1/2

###
integrate(\(x) sin(x) / (exp(x) * (exp(x) - 1)), 0, Inf)
pi/2 / tanh(pi) - 1

###
integrate(\(x) x * sin(x) / (exp(x) - 1), 0, Inf)
id = seq(10000)
2 * sum(id/(id^2 + 1)^2)
# TODO: see file: Sums.Fractions.Higher.R;
(pracma::psi(1, 1i) - pracma::psi(1, - 1i)) * 1i / 2


###
integrate(\(x) x * sin(2*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, 2i) - pracma::psi(1, - 2i)) * 1i / 2
#
k = sqrt(2)
integrate(\(x) x * sin(k*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i / 2


###
n = sqrt(3); # n >= 0 !!!
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
# TODO
id = seq(10000)
gamma(n+1) * sum(sin((n+1)*pi/2 - (n+1)*atan(id)) / (id^2 + 1)^((n + 1)/2))


### Derivation:
integrate(\(x) sin(x) / exp(x), 0, Inf)
1/2
#
integrate(\(x) x*sin(x) / exp(x), 0, Inf)
1/2
#
integrate(\(x) x^2*sin(x) / exp(x), 0, Inf)
1/2
#
sapply(seq(9), \(n) round(integrate(\(x) x^n * sin(x) / exp(x), 0, Inf)$value, 3))
sapply(seq(9), \(n) Im(gamma(n + 1)/(1-1i)^(n + 1)))
#
sapply(seq(9), \(n) round(integrate(\(x) x * sin(x) / exp(n*x), 0, Inf)$value, 5))
sapply(seq(9), \(n) round(Im(gamma(2)/(n - 1i)^(1 + 1)), 5))
sapply(seq(9), \(n) round(2*n*gamma(2)/(n^2 + 1)^(1 + 1), 5))

###
k = sqrt(5)
integrate(\(x) x * sin(x) / exp(k*x), 0, Inf)
2*k*gamma(2)/(k^2 + 1)^(1 + 1)

###
k = sqrt(5); n = sqrt(3)
integrate(\(x) x^n * sin(x) / exp(k*x), 0, Inf)
Im(gamma(n+1)/(k - 1i)^(n + 1))
gamma(n+1) / (k^2 + 1)^((n + 1)/2) * sin((n+1)*pi/2 - (n+1)*atan(k))


#######################
#######################

### I( atan(x) / (exp(x) - 1) )
# Flammable Maths: A ""-relaxing Integral Experience
# https://www.youtube.com/watch?v=QVgfL8Le0I0
# Flammable Maths:: One Spicy Class of Integrals.
# https://www.youtube.com/watch?v=HDaQAHhwtLo


###
integrate(\(x) atan(x) / (exp(2*pi*x) - 1), 0, Inf)
1/2 - log(2*pi)/4

###
integrate(\(x) atan(x) / (exp(x) - 1), 0, Inf)
pi*log(gamma(1/(2*pi))) - (pi - 1/2)*log(2*pi) + 1/2

###
n = sqrt(3)
integrate(\(x) atan(x) / (exp(n*x) - 1), 0, Inf)
log(gamma(n/(2*pi))) * pi/n + (1 - pi/n)/2*log(2*pi/n) - log(2*pi) * pi/(2*n) + 1/2


### Transformations
integrate(\(x) - atan(x) / (exp(x) - 1), 0, Inf)
#
integrate(\(x) log(1 - exp(-x)) / (x^2 + 1), 0, Inf)
integrate(\(x) log(1 - exp(-1/x)) / (x^2 + 1), 0, Inf)

### Perverse Transformations
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - x/(x^2+2), 0, Inf)
log(2)/2 - pi*log(gamma(1/(2*pi))) + (pi - 1/2)*log(2*pi) - 1/2

### including also Euler's constant:
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - exp(-1/x)/x, 0, Inf)
- pi*log(gamma(1/(2*pi))) + (pi - 1/2)*log(2*pi) - 1/2 + constEuler
# numerically problematic:
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - exp(-x)/x, 0, Inf)

