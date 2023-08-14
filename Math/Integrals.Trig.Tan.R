


### Trig: TAN


### Helper Constants
Catalan = 0.915965594177219015054603514;


##################
##################

###
integrate(\(x) x * tan(x), 0, pi/4)
Catalan / 2 - 1/8 * pi*log(2)

### [by Parts]
integrate(\(x) x^2 * tan(x)^2, 0, pi/4)
- 1/3*(pi/4)^3 + pi^2 / 16 + 1/4 * pi*log(2) - Catalan


###
integrate(\(x) x^2 * tan(x), 0, pi/4)
- 21/64 * zeta(3) - 1/32 * pi^2 * log(2) + 1/4 * pi * Catalan

### [by Parts]
integrate(\(x) x^3 * tan(x)^2, 0, pi/4)
(pi/4)^3 - 1/4*(pi/4)^4 + 63/64 * zeta(3) +
	+ 3/32 * pi^2 * log(2) - 3/4 * pi * Catalan;


###
integrate(\(x) tan(x) / x, 0, pi/4)
# TODO: ???


###
integrate(\(x) x / atan(x), 0, 1)
integrate(\(x) (tan(x) + tan(x)^3)/x, 0, pi/4)

# TODO: ???


#######################
#######################

############
### ATAN ###

### I( atan(x^n) )
n = sqrt(5)
integrate(\(x) atan(x^n), 0, 1)
pi/4 - integrate(\(x) n*x^n / (x^(2*n) + 1), 0, 1)$value
pi/4 - (digamma(((n+1)/(2*n) + 1)/2) - digamma((n+1)/(4*n))) / 4


### I( x^p * atan(x^n) )
n = sqrt(5); p = sqrt(7);
integrate(\(x) x^p * atan(x^n), 0, 1)
pi/(4*(p+1)) - integrate(\(x) n/(p+1) * x^(n+p) / (x^(2*n) + 1), 0, 1)$value
pi/(4*(p+1)) - (digamma(((n+p+1)/(2*n) + 1)/2) - digamma((n+p+1)/(4*n))) / (4*(p+1))
pi/(4*(p+1)) - (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1))


### Lim: p -> -1
n = 1/sqrt(3)
integrate(\(x) atan(x^n) / x, 0, 1)
Catalan / n;


### I( x^p * atan(x^n) * log(x) )
n = sqrt(5); p = sqrt(3);
# works only for n > 0;
integrate(\(x) x^p * atan(x^n) * log(x), 0, 1)
- pi/(4*(p+1)^2) + (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1)^2) +
	- (pracma::psi(1, ((p+1)/n + 3)/4) - pracma::psi(1, ((p+1)/n + 1)/4)) / (16*n*(p+1))


########################
########################

### ATAN-Fractions


### I( (atan(x) - pi/4) / (x^2 - 1) )
integrate(function(x) (atan(x) - pi/4) / (x^2 - 1), 0, 1)
Catalan / 2

# Extra-Component: pi/8*log(2);


### I( atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)) )
# Maths 505: Ahmed's monster integral: solution using Feynman's technique
# https://www.youtube.com/watch?v=Rq-bYzyoSwU


###
b = sqrt(5)
integrate(\(x) (2*x^2 + b^2) * atan(sqrt(x^2 + b^2)) /
	((x^2 + 1) * (x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
pi/4 * (3/2*pi + atan(sqrt(b^2 - 1)) - 2*atan(sqrt((b^2+1)/(b^2-1))) +
	- 2*atan(sqrt((b^2+1)*(b^2-1)))) / sqrt(b^2 - 1)

# Note: the I() collapses for b^2 - 1 == 1 to:
integrate(\(x) atan(sqrt(x^2 + 2)) / ((x^2 + 1) * sqrt(x^2 + 2)), 0, 1)
pi^2/8 * (7/4 - 4/3)


###
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)
pi/4 * (3/2*pi + atan(sqrt(b^2 - 1)) - 2*atan(sqrt((b^2+1)/(b^2-1))) +
	- 2*atan(sqrt((b^2+1)*(b^2-1)))) / sqrt(b^2 - 1) +
	- integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)$value

### [redundant]
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
pi/4*(3/2*pi + atan(sqrt(b^2-1)) - 2*atan(sqrt((b^2+1)*(b^2-1))) +
	- 2*atan(sqrt((b^2+1)/(b^2-1)))) / sqrt(b^2 - 1) +
	- integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)$value


# D(I): 1st I
# I( atan(t * sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)) ) on [0,1]
pi/4 / ((b^2-1)*t^2 + 1) +
	- t * atan(t / sqrt(b^2*t^2 + 1)) / sqrt(b^2*t^2 + 1) / ((b^2-1)*t^2 + 1)
# t => 1/t; dt => - dt / t^2;
# atan(1 / sqrt(t^2 + b^2)) / sqrt(t^2 + b^2) / (t^2 + (b^2-1))
# (pi/2 - atan(sqrt(t^2 + b^2))) / sqrt(t^2 + b^2) / (t^2 + (b^2-1))

# D(I): 2nd I
# I( atan(t * sqrt(x^2 + b^2)) / ((x^2 + b^2-1) * sqrt(x^2 + b^2)) ) on [0,1]
# I: 1 / ((x^2 + b^2-1) * (t^2*x^2 + t^2*b^2 + 1))
# ( 1 / (x^2 + b^2-1) - t^2 / (t^2*x^2 + t^2*b^2 + 1) ) / (t^2 + 1)
(pi/2 - atan(sqrt(b^2-1))) / ((t^2 + 1)*sqrt(b^2 - 1)) +
	- t * atan(t / sqrt(b^2*t^2 + 1)) / sqrt(b^2*t^2 + 1) / (t^2 + 1)


### Helper:
b = 1/sqrt(5)
b = sqrt(5)
integrate(\(x) 1 / ((x^2+1) * sqrt(x^2 + b^2)), 0, 1)
integrate(\(x) 1 / (x^2 + b^2 - 1), sqrt(b^2+1), Inf)
# b^2 > 1
pi/2/sqrt(b^2-1) - atan(sqrt((b^2+1)/(b^2-1)))/sqrt(b^2-1)
# b^2 < 1
(log(sqrt(b^2+1) + sqrt(1 - b^2)) - log(b) - log(2)/2) / sqrt(1 - b^2)

#
integrate(\(x) 1 / ((x^2 + b^2-1) * sqrt(x^2 + b^2)), 0, 1)
integrate(\(x) x / (((b^2-1)*x^2 + 1) * sqrt(b^2*x^2 + 1)), 1, Inf)
integrate(\(x) 1 / ((b^2-1)*x^2 + 1), sqrt(b^2+1), Inf)
# b^2 > 1
pi/2/sqrt(b^2-1) - atan(sqrt((b^2+1)*(b^2-1)))/sqrt(b^2-1)


###
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 1) * sqrt(x^2 + b^2)), 0, 1)

# D(b)
# D( atan(sqrt(x^2 + b^2)) * sqrt(x^2 + b^2) / (x^2 + 1) )
# b / ((x^2 + 1) * (x^2 + b^2+1)) + b * atan(sqrt(x^2 + b^2)) / sqrt(x^2 + b^2) / (x^2 + 1)
# ???

### TODO:
integrate(\(x) atan(sqrt(x^2 + b^2)) / ((x^2 + 2) * sqrt(x^2 + b^2)), 0, 1)
# much work;


#################

### I( atan(sqrt(x^2 + b^2)) / (x^2 + b^2)^(3/2) )
lim = c(1/5, 1/3)
integrate(\(x) atan(sqrt(x^2 + 1)) / (x^2 + 1)^(3/2), lim[1], lim[2])
x = lim
diff(x * atan(sqrt(x^2 + 1)) / sqrt(x^2 + 1) + atan(x) - sqrt(2) * atan(x/sqrt(2)))


### Gen:
lim = c(1/5, 1/3)
b = sqrt(5)
integrate(\(x) atan(sqrt(x^2 + b^2)) / (x^2 + b^2)^(3/2), lim[1], lim[2])
x = lim
diff(x * atan(sqrt(x^2 + b^2)) / sqrt(x^2 + b^2) + b*atan(x/b) +
	- sqrt(b^2 + 1) * atan(x/sqrt(b^2 + 1))) / b^2


#######################

### I( atan(tan(x)^2) )
# Maths 505: An outrageous journey of integration: int 0 to pi/4 arctan(cot^2(x))
# https://www.youtube.com/watch?v=VUGlU_dSgPY

###
integrate(\(x) atan(tan(x)^2), 0, pi/4)
integrate(\(x) atan(x^2) / (x^2 + 1), 0, 1)
pi^2/16 - (digamma(3/8 + 1/2) - digamma(3/8)) *
	(digamma(1/8 + 1/2) - digamma(1/8)) / 32;


### I( x * atan(x) / (x^4 + 1) )
integrate(\(x) x * atan(x) / (x^4 + 1), 0, 1)
(digamma(7/8) - digamma(3/8)) * (digamma(5/8) - digamma(1/8)) / 64;

# Note:
# - uses double integral: on [0,1] x [0,1]
#   I(I( u^2 / ((u^2*t^2 + 1) * (u^4 + 1)) ));


### I( x^3 * atan(x) / (x^4 + 1) )
integrate(\(x) x^3 * atan(x) / (x^4 + 1), 0, 1)
(digamma(3/8 + 1/2) - digamma(3/8))^2 / 128 +
	- (digamma(1/8 + 1/2) - digamma(1/8))^2 / 128 +
	+ Catalan/2;

# Helper:
integrate(\(x) atan(x) / x, 0, 1)
Catalan

# I( atan(x^n) / x )
n = sqrt(3);
integrate(\(x) atan(x^n) / x, 0, 1)
Catalan / n;

### Fraction Decomposition:
# u^4 / ((u^2*t^2 + 1) * (u^4 + 1))
# 1 / (u^2*t^2 + 1) - Fr0;
# =>
# Fr0 + (u^2*t^2 - 1) / ((u^4 + 1) * (t^4 + 1));
# where Fr0 = 1 / ((u^2*t^2 + 1) * (u^4 + 1));
# and (u^2*t^2 - 1) / Prod() is decomposable into 2 separate integrals:
# I( u^2 / (u^4 + 1) )^2 -  I( 1 / (u^4 + 1) )^2;


### Note: Lim p -> 0
p = 1E-4
integrate(\(x) atan(x^(1/p)) / p, 0, 1)
pi/(4*p) - integrate(\(x) 1/p^2 * x^(1/p) / (x^(2/p) + 1), 0, 1)$value
pi/(4*p) - integrate(\(x) 1/p * x^p / (x^2 + 1), 0, 1)$value
- (pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16
Catalan

