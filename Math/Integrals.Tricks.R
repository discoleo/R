########################
###
### Leonard Mada
### [the one and only]
###
### Integral Tricks
###
### draft v.0.2d


### various Integral Tricks

### draft v.0.2d:
# - more exponential tricks:
#   I( atan(exp(x)) ) dx;
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
integrate(function(x) log(x) / (x^2 + b*x + 1), lower=1/lim, upper=lim)

b = c(5, 2)
lim = 2
integrate(function(x) log(x) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)
integrate(function(x) log(b[1]) / (x^2 + b[2]*b[1]*x + b[1]^2), lower=b[1]/lim, upper=b[1]*lim)


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

### sin(log(x))
a = 3
b = 3
integrate(function(x) sin(log(x)) / (x^2+b*x+1), lower=1/a, upper=a)
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
pi/2*log(a)


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

a = 3
b = 3
integrate(function(x) atan(sqrt(x)) / (x^2+b*x+1), lower=1/a, upper=a)
integrate(function(x) pi/4 / (x^2+b*x+1), lower=1/a, upper=a)


a = 3
integrate(function(x) atan(x) * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(function(x) pi/4 * sqrt(x+1/x) / x, lower=1/a, upper=a)
integrate(function(x) pi/2 * sqrt(x^4+1) / x^2, lower=1/sqrt(a), upper=sqrt(a))


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

a = 3
integrate(function(x) atan(exp(x)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x)), lower=-a, upper=a)
a * pi/2

integrate(function(x) atan(exp(x^3)), lower=-a, upper=a)
integrate(function(x) atan(exp(-x^3)), lower=-a, upper=a)
a * pi/2

