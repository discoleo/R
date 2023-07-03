


Euler   = 0.57721566490153286060651209008240243079;


polylog = function(x, n=2) {
	sapply(x, \(x) if(x == 1) pracma::zeta(n)
		else if(x > 0.55) sum(x^seq(128)/seq(128)^n)
		else pracma::polylog(x, n));
}

################
################

### I( Polylog(x^p, 2) )
# 1.) Maths 505: Integrating the dilogarithm!
# https://www.youtube.com/watch?v=Al39YFdEnQk
# 2.) Michael Penn: a surprisingly interesting sum -- 2 ways!
# https://www.youtube.com/watch?v=vJVf14KSks0
# (2nd method uses Polylog)

###
integrate(polylog, n=2, lower=0, upper=1)
pi^2/6 - 1

###
integrate(\(x) polylog(x^(1/2), 2), 0, 1)
pi^2/6 - 3/4

###
integrate(\(x) polylog(x^(1/3), 2), 0, 1)
pi^2/6 - (1+1/2+1/3)/3

###
integrate(\(x) polylog(x^(1/4), 2), 0, 1)
pi^2/6 - (1+1/2+1/3+1/4)/4

###
p = 1/5
integrate(\(x) polylog(x^p, 2), 0, 1)
pi^2/6 - p*(digamma(1/p + 1) + Euler)

###
p = 1/sqrt(5)
integrate(\(x) polylog(x^p, 2), 0, 1)
pi^2/6 - p*(digamma(1/p + 1) + Euler)


### Zeta(3)
integrate(\(x) polylog(x^(1/2), 3), 0, 1)
pracma::zeta(3) - pi^2/12 + 3/8

