


#####################

### Simple

### I( sin(x)^n / x^n * exp(-x) )
# Maths 505: Feynman's technique is incredibly overpowered!
# https://www.youtube.com/watch?v=_d1PJiIm-lE

integrate(\(x) sin(x)^2/x^2 / exp(x), 0, Inf)
atan(2) - log(5)/4

###
integrate(\(x) sin(x)^3/x^3 / exp(x), 0, Inf)
atan(3) - 3*log(5)/8


# Derivation:
t = sqrt(2)
integrate(\(x) sin(x)^3 / exp(t*x), 0, Inf)
integrate(\(x) (3*sin(x) - sin(3*x))/4 / exp(t*x), 0, Inf)
Im(1/(t + 3i) - 3/(t + 1i)) / 4
# =>
3/4*(1/(t^2 + 1) - 1/(t^2 + 9))

#
integrate(\(x) sin(x)^3/x / exp(x), 0, Inf, rel.tol=1E-10)
# Constant of Integration = pi/4;
pi/4 - (atan(1) - 1/3*atan(1/3))*3/4
(atan(1/3) + pi/4)/4


###
t = sqrt(2)
integrate(\(x) sin(x)^5 / exp(t*x), 0, Inf)
integrate(\(x) (sin(5*x) + 20*sin(x)^3 - 5*sin(x))/16 / exp(t*x), 0, Inf)
integrate(\(x) (sin(5*x) - 5*sin(3*x) + 10*sin(x))/16 / exp(t*x), 0, Inf)
- Im(1/(t + 5i) - 5/(t + 3i) + 10/(t + 1i)) / 16
# =>
5/16*(1/(t^2 + 25) - 3/(t^2 + 9) + 2/(t^2 + 1))


###
t = sqrt(2)
integrate(\(x) sin(x)^7 / exp(t*x), 0, Inf)
Im(1/(t + 7i) - 7/(t + 5i) + 21/(t + 3i) - 35/(t + 1i)) / 4^3
# =>
-7/4^3 * (1/(t^2 + 49) - 5/(t^2 + 25) + 9/(t^2 + 9) - 5/(t^2 + 1))


###
t = sqrt(2)
integrate(\(x) sin(x)^9 / exp(t*x), 0, Inf)
- Im(1/(t + 9i) - 9/(t + 7i) + 36/(t + 5i) - 84/(t + 3i) + 126/(t + 1i)) / 4^4
9/4^4 * (1/(t^2 + 81) - 7/(t^2 + 49) + 20/(t^2 + 25) - 28/(t^2 + 9) + 14/(t^2 + 1))


### I[n]( atan(...) )
# x/a * atan(x/a) - 1/2*log(x^2 + a^2)
a = sqrt(3)
integrate(\(x) 1/a * atan(x/a), 0, 1)
1/a * atan(1/a) - log(a^2 + 1)/2 + log(a)

### I[2]
a = sqrt(3)
# x^2/(2*a) * atan(x/a) - a/2 * atan(x/a) - 1/2*x*log(x^2 + a^2) + x/2;
integrate(\(x) x/a * atan(x/a) - 1/2*log(x^2 + a^2), 0, 1)
1/(2*a) * atan(1/a) - a/2 * atan(1/a) - 1/2*log(a^2 + 1) + 1/2


#####################
#####################

### I( sin(x) / (exp(x) - 1))
# Maths 505: One of the coolest integrals on YouTube!
# https://www.youtube.com/watch?v=I7a-6VwKhgo

###
integrate(\(x) sin(x) / (exp(x) - 1), 0, Inf)
pi/2 / tanh(pi) - 1/2
(digamma((1/2+1)/2) - digamma(1/4))/2 - 1/2
- (pracma::psi(0, 1i) - pracma::psi(0, - 1i)) * 1i / 2 - 1

### I( sin(k*x) / (exp(x) - 1) ) on [0, Inf]
k = sqrt(3)
integrate(\(x) sin(k*x) / (exp(x) - 1), 0, Inf)
- (pracma::psi(0, 1i*k) - pracma::psi(0, - 1i*k)) * 1i / 2 - 1/k;


###
integrate(\(x) sin(x) / (exp(x) * (exp(x) - 1)), 0, Inf)
pi/2 / tanh(pi) - 1

###
integrate(\(x) x * sin(x) / (exp(x) - 1), 0, Inf)
id = seq(10000)
2 * sum(id/(id^2 + 1)^2)
# TODO: see file: Sums.Fractions.Higher.R;
(pracma::psi(1, 1i) - pracma::psi(1, - 1i)) * 1i / 2


### I( x * sin(k*x) / (exp(x) - 1) )
integrate(\(x) x * sin(2*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, 2i) - pracma::psi(1, - 2i)) * 1i / 2
#
k = sqrt(2)
integrate(\(x) x * sin(k*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i / 2


### Generalization:
n = sqrt(3); # n >= 0 !!!
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
# TODO: formula for sum();
id = seq(10000)
gamma(n+1) * sum(sin((n+1)*pi/2 - (n+1)*atan(id)) / (id^2 + 1)^((n + 1)/2))


### n = ODD
n = 3;
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
(pracma::psi(n, 1i) - pracma::psi(n, - 1i)) * 1i / 2
#
n = 5; k = sqrt(3)
integrate(\(x) x^n * sin(k*x) / (exp(x) - 1), 0, Inf)
(pracma::psi(n, 1i*k) - pracma::psi(n, - 1i*k)) * 1i / 2


### n = EVEN

# TODO: ???
n = 2;
integrate(\(x) x^n * sin(x) / (exp(x) - 1), 0, Inf)
#
n = 4; k = sqrt(3)
integrate(\(x) x^n * sin(k*x) / (exp(x) - 1), 0, Inf)


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


#########################
#########################

### Trig - Exp

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

