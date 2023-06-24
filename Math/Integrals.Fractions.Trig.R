

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
# can be expanded explicitly
# TODO: alternative?


###
n = 10
integrate(function(x) 1 / (cos(x)^n + sin(x)^n), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^(n/2 - 1) / (x^n + 1), lower=0, upper=1)

