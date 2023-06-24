

Euler   = 0.57721566490153286060651209008240243079;

### I on [0, 1]
# I( x^p / (x^n + 1))
int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}

#############

### Based on:
# Integrals.Exercises.R &
# Integrals.Fractions.Unity.Derived.R
# - but those files need major refactoring;
# - see also:
#   Maths 505: A STUNNING integral from the JEE (main) exam
#   https://www.youtube.com/watch?v=mjM2GVJusU4


###
integrate(function(x) 1 / (cos(x)^4 + sin(x)^4), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1) / (x^4 + 1), lower=0, upper=1)
pi/(2*sqrt(2))


###
integrate(function(x) 1 / (cos(x)^6 + sin(x)^6), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^2 / (x^6 + 1), lower=0, upper=1)
pi/2


###
integrate(function(x) 1 / (cos(x)^8 + sin(x)^8), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^3 / (x^8 + 1), lower=0, upper=1)
sum(c(1,3,3,1) * sapply(c(6,4,2,0), int.FrU01, n=8))
pi*(1/tan(pi/16) + 1/tan(7*pi/16) + 3/tan(pi*3/16) + 3/tan(pi*5/16))/16
# can be expanded explicitly

# TODO: alternative?
integrate(function(x) (x^2 + 4) / (x^4 + 4*x^2 + 2), lower=0, upper=-Inf)
# integrate(function(x) 2 / (x^2 + 4 + 2*sqrt(2)), lower=0, upper=Inf)$value +
pi/sqrt(4+2*sqrt(2)) +
integrate(function(x) (4 - sqrt(2)) / (x^4 + 4*x^2 + 2), lower=0, upper=Inf)$value


###
n = 10; n2 = 2*n;
integrate(function(x) 1 / (cos(x)^n + sin(x)^n), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^(n/2 - 1) / (x^n + 1), lower=0, upper=1)
pi*(1/tan(pi/n2) + 1/tan(9*pi/n2) + 4/tan(pi*3/n2) + 4/tan(pi*7/n2) +
	+ 6/tan(pi*5/n2))/n2

