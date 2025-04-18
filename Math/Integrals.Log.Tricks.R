########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Various Tricks
###
### draft v.0.4b

# TODO: refactor & cleanup;


### various Integral Tricks
# - only simple/basic tricks;

### draft v.0.4a:
# - [refactor] renamed & split;
# - moved log-fractions to new file:
#   Integrals.Log.Fractions.R;
### draft v.0.3d:
# - moved section on Moebius Transform
#   to separate file:
#   Integrals.Moebius.Transform.R;
### draft v.0.3c - v.0.3c-ex2:
# - definite integrals using Moebius Transform:
#   I( log(4*x-3) / (x*(3*x+1)) ) on c(1, 4);
#   I( log(5*x-3) / (x*(x+1)) ) on c(1, 3);
#   I( log(x-3) / (x*(3*x+1)) ) on c(5, 8);
### draft v.0.3a - v.0.3b:
# - Transformations:
#   f((a*x + b) / (c*x + d)) = - f(x);
# - added 2nd example;
# - acos() variant; [v.0.3b]
### draft v.0.2f:
# - started work on:
#   I( log(x-1) / x ) dx;
#   [see v.0.3a-d for related results]
### draft v.0.2d - v.0.2e:
# - more exponential tricks:
#   I( atan(exp(x)) ) dx;
#   I( atan(exp(x)) / (x^2 + 1) ) dx;
### draft v.0.2c:
# - partial results:
#   I( atan(x)^3 / x ) dx;
### draft v.0.2b:
# - just for completion:
#   I( atan(x) / x ) dx;
### draft v.0.2a:
# - sin(log(x)) type:
#   I( sin(log(x)) / (x^2+b*x+1) ) dx;
### draft v.0.1i - v.0.1k:
# - atan() trick:
#   I( atan(x) / (x^2+b*x+1) ) dx;
#   I( atan(sqrt(x)) / (x^2+b*x+1) ) dx;
#   I( atan(x) / sqrt(x^4+b*x^3+b*x+1) ) dx;
### draft v.0.1f - v.0.1h:
# - experiments with complex integrals;
#   I( log(x) / (x^2 - 1) ) dx over or using complex path;
### draft v.0.1e:
# - extension of the log() trick;


##################
##################

##################
### Logarithms ###
##################

### Example 1:
# - based on: "This trick is new to me!"
#   https://www.youtube.com/watch?v=JCuplIQ6JG4


### I (log(x) / (x^2 + b*x + 1)) dx
# lower = 1/a, upper = a, I = 0;
b = 3
lim = 2
# I == 0
integrate(function(x) log(x) / (x^2 + b*x + 1), lower=1/lim, upper=lim)

# log becomes constant;
solve.intFr = function(b, lim) {
	if(b[2] == 2) {
		# Note: Case b[2] == -2 skipped;
		log(b[1]) * (1/(1/lim + 1) - 1/(lim + 1)) / b[1];
	} else {
		d = (b[2]^2 - 4) / 4;
		if(d > 0) {
			d  = sqrt(d);
			lg = (lim + (b[2]/2 - d))/(lim + (b[2]/2 + d)) *
				(1/lim + (b[2]/2 + d))/(1/lim + (b[2]/2 - d));
			log(b[1]) * log(lg) / (2*b[1]*d);
		} else {
			d = sqrt(- d);
			log(b[1]) * (atan((lim + b[2]/2)/d) - atan((1/lim + b[2]/2)/d)) / (b[1]*d);
		}
	}
}
#
b = c(5, 2)
b = c(5, 3/2)
# b = c(5, 3)
lim = 2
integrate(function(x) log(x) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)
integrate(function(x) log(b[1]) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)
solve.intFr(b, lim)


### Higher Powers
a  = 3
b2 = 5
#
b  = c(1, b2, b2*a, a^3)
bf = c(1, b2/a, b2/a, 1)
# Note: polynomial fraction can be simplified;
integrate(function(x) sapply(x,
	function(x) (x+a)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)

integrate(function(x) sapply(x,
	function(x) log(a)/a * (x+1)/sum(bf*x^(3:0))), lower=0, upper=Inf)

###
integrate(function(x)
	sapply(x, function(x) sqrt(x)*log(x)/sum(b*x^(3:0))), lower=0, upper=Inf)
integrate(function(x)
	sapply(x, function(x) log(a)*sqrt(a)/a^2 *sqrt(x)/ sum(bf*x^(3:0))), lower=0, upper=Inf)


#########
### Ex 2:
a = 3
b1n = 4
b1 = 2
# b2 can be added as well;

integrate(function(x) (x^2 + b1n*x + a^2)*log(x)/
	(x^4 + b1*x^3 + b1*a^2*x + a^4), lower=0, upper=Inf)
integrate(function(x) log(a)/a * (x^2 + b1n/a*x + 1)/
	(x^4 + b1/a*x^3 + b1/a*x + 1), lower=0, upper=Inf)

### only x:
integrate(function(x) x*log(x)/
	(x^4 + b1*x^3 + b1*a^2*x + a^4), lower=0, upper=Inf)
integrate(function(x) log(a)/a^2 * x /
	(x^4 + b1/a*x^3 + b1/a*x + 1), lower=0, upper=Inf)


#########
### Ex 3:
integrate(function(x) (x^3+8)*log(x)/(x^5 + 3*x^4 + 24*x + 32), lower=0, upper=Inf)
integrate(function(x) log(2)/2 * (x^3+1)/(x^5 + 3/2*x^4 + 3/2*x + 1), lower=0, upper=Inf)


###############

library(pracma)

# Finite Integrals

### I( log(x) / (x^2 - 1) ) dx
line_integral(function(x)  log(x) / (x^2+1), c((1-1i)/sqrt(2), (1+1i)/sqrt(2)))
line_integral(function(x) -log(-x*1i)*1i / (x^2-1), c((-1+1i)/sqrt(2), (1+1i)/sqrt(2)))
line_integral(function(x) (-1i*log(x) - pi/2) / (x^2-1), c((-1+1i)/sqrt(2), (1+1i)/sqrt(2)))


a = 2
line_integral(function(x)  log(x) / (x^2+1), c(1/a*(1-1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) -log(-x*1i)*1i / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) (-1i*log(x) - pi/2) / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))

a = 2
line_integral(function(x) log(x) / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
line_integral(function(x) 1i*pi/2 / (x^2-1), c(1/a*(-1+1i)/sqrt(2), a*(1+1i)/sqrt(2)))
1i*pi/4 * log((a*(1+1i) - sqrt(2))*(1/a*(-1+1i) + sqrt(2)) /
	((1/a*(-1+1i) - sqrt(2))*(a*(1+1i) + sqrt(2))))

lim = c(1.1, 3) # lim[1] > 1; # strictly greater than 1!
line_integral(function(x) log(x) / (x^2-1), c(lim[1]*1i, lim[2]*1i))
1i*pi/4 * log((lim[2]*1i - 1)*(lim[1]*1i + 1) /
	((lim[2]*1i + 1)*(lim[1]*1i - 1))) +
	+ 1i * integrate(function(x) log(x) / (x^2+1), lower=1/lim[2], upper=lim[1])$value


lim = c(1.1, 3) # lim[1] > 1; # strictly greater than 1!
integrate(function(x) log(x) / (x^2-1), lower=lim[1], upper=lim[2])
1i * line_integral(function(x) log(-1i*x) / (x^2+1), c(lim[1]*1i, lim[2]*1i))
1i * line_integral(function(x) log(x) / (x^2+1), c(lim[1]*1i, lim[2]*1i)) +
	+ pi/2 * (atan(lim[2]*1i) - atan(lim[1]*1i))


### TODO:
# - find full path;


########################

# Michael Penn: I really like this integral
# https://www.youtube.com/watch?v=5sxtrPgHFw8

###
n = 4
a = sqrt(2)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


###
n = sqrt(5)
a = sqrt(2)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


###
n = 3
a = 2^(1/3)
integrate(function(x) x^a * abs(log(x))^n, lower=0, upper=1)
gamma(n+1) / (a+1)^(n+1)


########################
########################

#####################
### Trigonometric ###
#####################

### sin(log(x))
a = 3
b = 3
integrate(function(x) sin(log(x)) / (x^2 + b*x + 1), lower=1/a, upper=a)
# 0
integrate(function(x) cos(log(exp(pi/2)/x)) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) cos(log(exp(pi/2)*x)) / (x^2+b*x+1), lower=1/a, upper=a)
# 0

integrate(function(x) sin(log(x)) / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)
# 0


########################
########################

### ARCTAN

a = 3
integrate(function(x) atan(x) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 / x, lower=1/a, upper=a)
pi/2 * log(a)


a = 3
b = 3
integrate(function(x) atan(x) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / (x^2+b*x+1), lower=1/a, upper=a)


a = 3
b = 3
integrate(function(x) x*atan(x) / (x^4+b*x^3+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 * x / (x^4+b*x^3+b*x+1), lower=1/a, upper=a)


### Sqrt
a = 3
integrate(function(x) atan(sqrt(x)) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 / x, lower=1/a, upper=a)
pi/2 * log(a)

a = 3
b = 3
integrate(function(x) atan(sqrt(x)) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / (x^2+b*x+1), lower=1/a, upper=a)
d = sqrt(b^2 - 4);
pi/(4*d) * log((2*a + b - d)/(2*a + b + d)*(2/a + b + d)/(2/a + b - d))


###
a = 3
integrate(\(x) atan(x) * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(\(x) pi/4 * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(\(x) pi/2 * sqrt(x^4+1) / x^2, lower=1/sqrt(a), upper=sqrt(a))

### I( atan(x) * sqrt(x + 1/x) / x )
integrate(\(x) atan(x) * sqrt(x + 1/x) / x - pi/2 / x^(1/2), 0, Inf)
integrate(\(x) pi/2 * sqrt(x^4+1) / x^2 - pi/2 / x^2 - pi/2, 0, Inf)
pi * gamma(-1/4)*gamma(1/2 + 1/4) / gamma(1/2) / 4;

### on [0, 1]
integrate(\(x) atan(x) * sqrt(x + 1/x) / x, 0, 1)
integrate(\(x) 2^(-1/2) * x / sin(x)^(3/2), 0, pi/2)
# TODO


a = 3
b = 3
integrate(function(x) atan(x) / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / sqrt(x^4+b*x^3+b*x+1), lower=1/a, upper=a)


### Power 3
a = 3
integrate(function(x) atan(x)^3 / x, lower=1/a, upper=a)
integrate(function(x) (3*pi/4*atan(x)^2 - pi^3/32) / x, lower=1/a, upper=a)
integrate(function(x) 3*pi/4*atan(x)^2 / x, lower=1/a, upper=a)$value - pi^3/16*log(a)


### Exp
a = 3
integrate(function(x) (1 - exp(-pi/2))*exp(atan(x)) / x, lower=1/a, upper=a)
integrate(function(x) (exp(pi/2) - 1)*exp(-atan(x)) / x, lower=1/a, upper=a)

integrate(function(x) (exp(atan(x)) - exp(pi/2)*exp(-atan(x))) / x, lower=1/a, upper=a)
# 0

a = 3
integrate(function(x) exp(atan(x)) / (x * (exp(atan(x)) + exp(-atan(x)))), lower=1/a, upper=a)
# TODO


### atan(exp(x))
a = 3
integrate(function(x) atan(exp(x)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)), lower=-a, upper=a)
a * pi/2

integrate(function(x) atan(exp(x^3)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x^3)), lower=-a, upper=a)
a * pi/2


a = 3
integrate(function(x) atan(exp(x)) / (x^2 + 1), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)) / (x^2 + 1), lower=-a, upper=a)
pi/2 * atan(a)

integrate(function(x) atan(exp(x)) * x^2 / (x^2 + 1), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)) * x^2 / (x^2 + 1), lower=-a, upper=a)
pi/2 * (a - atan(a))


####################
####################

a = 2
integrate(function(x) log(x-1) / x, lower=a, upper=a+1)
integrate(function(x) log(x-1) / (x*(x-1)), lower=a/(a-1), upper=(a+1)/a)

i.f = function(x) 1/2*log(x-1)^2
i.f((a+1)/a) - i.f(a/(a-1))
integrate(function(x) log(x-1) / x, lower=a, upper=a+1)$value +
	+ integrate(function(x) log(x-1) / x, lower=a/(a-1), upper=(a+1)/a)$value
1/2*log(a)^2 - 1/2*log(a-1)^2

# TODO: both limits;


########################
########################
