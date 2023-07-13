########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Fractions: Other
###
### draft v.0.1a


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


### Helper Constants
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


################
################

### I( log(x) * log(1 - x^n) / x )
# Maths 505: A beautiful log-trig integral featuring an important constant
# https://www.youtube.com/watch?v=2bwvuWsQSEY
# [the intermediate integral]

###
integrate(\(x) log(x) * log(1 - x^2) / x, 0, 1)
pracma::zeta(3) / 4

###
n = sqrt(5)
integrate(\(x) log(x) * log(1 - x^n) / x, 0, 1)
pracma::zeta(3) / n^2

###
integrate(\(x) log(x) * log(1 + x) / x, 0, 1)
- pracma::zeta(3) * 3/4

###
n = sqrt(5)
integrate(\(x) log(x) * log(1 + x^n) / x, 0, 1)
- pracma::zeta(3) * 3/4 / n^2

### Derived:
integrate(\(x) log(x) * log(x^2 + x + 1) / x, 0, 1)
- pracma::zeta(3) * 8/9
#
integrate(\(x) log(x) * log(x^2 - x + 1) / x, 0, 1)
pracma::zeta(3) * 2/3


### Series: Pow = 2
integrate(\(x) log(x) * log(1 - x^2) / x^2, 0, 1)
pi^2 / 4 - 2*log(2)

###
n = sqrt(3)
integrate(\(x) log(x) * log(1 - x^(2*n)) / x^(n+1), 0, 1)
(pi^2 / 4 - 2*log(2)) / n^2

# Partial Variants:
# =>
integrate(\(x) log(x) * log(1 + x) / x^2 - log(x)/x + log(x)/2, 0, 1)
pi^2/12 - 2*log(2) + 0.5
integrate(\(x) log(x) * log(1 - x) / x^2 + log(x)/x + log(x)/2, 0, 1)
pi^2/6 - 1.5


### Other:
integrate(function(x) log(x) * log(1 - x), 0, 1)
2 - pi^2/6
#
integrate(function(x) log(x) * log(1 - x), 0, 1/2)
integrate(function(x) log(x) * log(1 - x), 1/2, 1)
integrate(function(x) log(1/2 - x) * log(1/2 + x), 0, 1/2)
1 - pi^2/12


################
################

### I( log((x - 1) / (x + 1)) / x ) on [1, Inf]
# Maths 505: A RIDICULOUSLY AWESOME LOG INTEGRAL!!!
# https://www.youtube.com/watch?v=tTu6hedSlm0

# numerical issue:
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, Inf)
integrate(\(x) log((exp(x) - 1) / (exp(x) + 1)), 0, 100)
# alternative:
integrate(\(x) log((x - 1) / (x + 1)) / x, 1, Inf)
integrate(\(x) log((1 - x) / (x + 1)) / x, 0, 1)
- pi^2/4

###
integrate(\(x) log(x) * log((1 - x) / (1 + x)) / x, 0, 1)
7/4 * pracma::zeta(3)

###
integrate(\(x) log(x)^2 * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol = 1E-8)
- 15/4 * pracma::zeta(4)

###
integrate(\(x) log(x)^3 * log((1 - x) / (1 + x)) / x, 0, 1, rel.tol = 1E-8)
93/8 * pracma::zeta(5)

### =>
integrate(\(x) log(x) * log(1 - x) / x, 0, 1)
pracma::zeta(3)
#
integrate(\(x) log(x) * log(1 + x) / x, 0, 1)
-3/4 * pracma::zeta(3)

### =>
integrate(\(x) log(x)^3 * log(1 - x) / x, 0, 1)
6 * pracma::zeta(5)
#
integrate(\(x) log(x)^3 * log(1 + x) / x, 0, 1)
-45/8 * pracma::zeta(5)

### Generalized:
p = sqrt(5)
integrate(\(x) abs(log(x))^p * log(1 - x) / x, 0, 1)
- gamma(p+1) * pracma::zeta(p+2)


### I( |log(x)|^p * log(1 - x^n) / x )
p = sqrt(3); n = sqrt(7);
integrate(\(x) abs(log(x))^p * log(1 - x^n) / x, 0, 1, rel.tol=1E-8)
- gamma(p+1) * pracma::zeta(p+2) / n^(p+1)


### I( |log(x)|^p * log(1 + x) / x )
# workout of log(1 - x^2) =>
p = sqrt(7)
integrate(\(x) abs(log(x))^p * log(1 + x) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1))
#
n = sqrt(5)
integrate(\(x) abs(log(x))^p * log(1 + x^n) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1)) / n^(p+1)

#
integrate(\(x) abs(log(x))^p * log(1 + x + x^2) / x, 0, 1, rel.tol=1E-8)
gamma(p+1) * pracma::zeta(p+2) * (1 - 1/3^(p+1))
#
integrate(\(x) abs(log(x))^p * log(1 - x + x^2) / x, 0, 1, rel.tol=1E-8)
- gamma(p+1) * pracma::zeta(p+2) * (1 - 1/2^(p+1)) * (1 - 1/3^(p+1))


### Simple:
integrate(\(x) log((1 - x) / (x + 1)), 0, 1)
- 2*log(2)

###
integrate(\(x) log(x) * log((1 - x) / (x + 1)), 0, 1)
2*log(2) - pi^2/12


################
################

# TODO: move to Trig

# Maths 505: Another INSANE integral!
# https://www.youtube.com/watch?v=KEDEzVqlAYU

integrate(function(x) atan(x) * log((1-x)/(1+x)), 0, 1)
pi^2/16 - pi/4*log(2) - Catalan;


# Gen 1: TODO
k = 2
integrate(function(x) atan(x)*log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


###
integrate(function(x) atan(x) / x, 0, 1)
Catalan

# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


##################

### I( log(x) * atan(x) / x )
# 1.) Maths 505: a superb integral sprinkled with some fourier analysis
# https://www.youtube.com/watch?v=aHCLS4l65VE
# 2.a) see also: Sums.Fractions.Higher.R;
# 2.b) Flammable Maths: An AMAZING Journey of Series Evaluation!
# Calculating Euler's Sum! [ Series pi^3/32 (-1)^k/(2k+1)^3 ]

integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32


### I( log(x) * atan(x) )
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2/48 - pi/4 + log(2)/2

# Derivation:
id = seq(0, 20000)
- sum((-1)^id * (1/(2*id+1) - 1/(2*id+2))) + pi^2/48
- sum((-1)^id / (2*id+1)) + pi^2/48 + log(2)/2
pi^2/48 - pi/4 + log(2)/2


### I( log(x) * atan(x^p) / x )
integrate(\(x) log(x) * atan(x^2) / x, 0, 1)
- pi^3 / 128

### Gen:
p = sqrt(5)
integrate(\(x) log(x) * atan(x^p) / x, 0, 1)
- pi^3 / (32 * p^2)

### Other:
integrate(\(x) log(x) * x / (x^2 + 1), 0, 1)
- pi^2/48

###
integrate(\(x) log(x^2 + 1) / x, 0, 1)
pi^2/24


####################
####################

### I( log(x) * log(x+1) )
# Maths 505: A surprisingly wonderful log integral
# https://www.youtube.com/watch?v=v6iNE-heksU

integrate(\(x) log(x) * log(x+1), 0, 1)
2 - 2*log(2) - pi^2/12

###
a = sqrt(5)
integrate(\(x) log(x) * log(x+a), 0, a)
integrate(\(x) a*log(a*x) * log(a*x+a), 0, 1)
a*(2 + 2*log(2)*log(a) - 2*log(2) - 2*log(a) + log(a)^2 - pi^2/12)


###
a = 1/sqrt(3)
integrate(\(x) log(x) * log(x+1), 0, a)
# TODO:
id = seq(10000)
log(a) * sum((-a)^(id+1) / (id*(id+1))) +
	sum((-1)^id * a^(id+1) / (id*(id+1)^2))

###
integrate(\(x) log(x) * log(x+1), 0, 1/2)
# TODO: find closed form?
Li2 = - integrate(\(x) x / (2*exp(x) + 1), 0, Inf)$value;
(1/2 - 3/2*log(3/2))*log(2) - 3/2*log(3/2) + Li2 + 1
# alternative for Li2:
id = seq(10000)
Li2 = sum((-1/2)^id / id^2);


###
integrate(\(x) log(x) * log(1-x) / x^2, 1/2, 1)
pi^2/12 + 2*log(2)^2 - 2*log(2)

#
integrate(\(x) log(x)^2 / x^2, 1/2, 1)
2*log(2)^2 - 4*log(2) + 2


######################
######################

### I( log(sqrt(x^2 + 1) + x) / (x^2 + 1) )
# Maths 505: An AWESOME log integral with a surprising result
# https://www.youtube.com/watch?v=a2On9rAvqks

integrate(\(x) log(sqrt(x^2 + 1) + x) / (x^2 + 1), 0, Inf)
integrate(\(x) log(sqrt(x^2 + 1) + 1) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - x) / (x^2 + 1), 0, Inf)
integrate(\(x) - log(sqrt(x^2 + 1) - 1) / (x^2 + 1), 0, Inf)
2*Catalan


### Derived / Sub-Integrals:
integrate(\(x) log((1 + cos(x)) / cos(x)), 0, pi/2)
2*Catalan

###
# - see also file: Integrals.Log.Trig.R
integrate(\(x) log(1 + cos(x)), 0, pi/2)
integrate(\(x) 4*log(cos(x)), 0, pi/4)$value + pi*log(2)/2
- pi*log(2)/2 + 2*Catalan


#######################
#######################

### I( x * log(x) / cosh(k*x)^2 )
# Maths 505: Feynman's trick solves an AWESOME integral
# https://www.youtube.com/watch?v=qZywLtQaDGQ

###
k = sqrt(5)
integrate(\(x) x * log(x) / cosh(k*x)^2, 0, 100)
- (log(2)*log(k) + 3/2*log(2)^2 - log(2)) / k^2

### Helper:
k = sqrt(3)
integrate(\(x) x * exp(k*x) / (exp(k*x) + 1)^2, 0, 100)
log(2) / k^2


#######################
#######################

### I( (x^n - 1) / log(x) )
# Michael Penn: My favorite way "around" Feynman's trick -- and a challenge for you!
# https://www.youtube.com/watch?v=EseFhG1QARM

n = sqrt(5)
integrate(\(x) (x^n - 1) / log(x), 0, 1)
log(n+1)

