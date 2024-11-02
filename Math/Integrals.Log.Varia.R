

#############
### Basic ###

### Bernoulli Integrals

### I( x^(-x) )
# on [0, Inf]
integrate(\(x) x^(-x), 0, Inf)
# TODO: ???


### on [0, 1]
integrate(\(x) x^(-x), 0, 1, rel.tol=1E-8)
id = seq(10); sum(1/id^id)
#
integrate(\(x) x^(x), 0, 1, rel.tol=1E-8)
id = seq(10); sum(rep(c(1,-1), 5)/id^id)

# TODO: Closed formula?


### Convoluted

###
integrate(\(x) x^x * (1-x)^(1-x), 0, 1, rel.tol=1E-12)
# TODO

###
integrate(\(x) x^(1-x) * (1-x)^x, 0, 1, rel.tol=1E-12)
# TODO


###########

###########
### EXP ###

### I( x^(-x) * exp(-x) )
integrate(\(x) abs(log(x))^log(x), 0, 1, rel.tol=1E-8)
integrate(\(x) x^(-x) * exp(-x), 0, Inf, rel.tol=1E-8)

# TODO: ???


### I( x^(-x) * exp(x) )
integrate(\(x) x^(-x) * exp(x), 0, Inf)
integrate(\(x) exp(x - x*log(x)), 0, Inf, rel.tol=1E-8)

# TODO: ???

# library(Rmpfr)
integrate(\(x) {
	x = mpfr(x, 240);
	as.numeric(x^(-x) * exp(x)); }, 0, Inf, rel.tol=1E-8)

