############
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.1b

### Double Integrals

### Examples:
# I( (x+y)*log(x+y) / (1+x*y) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

#####################
#####################

### Simple

### I( 1 / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1-x*y), 0, 1)$value), 0, 1)
pi^2 / 6

### I( 1 / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1+x*y), 0, 1)$value), 0, 1)
pi^2/12


### I( log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y), 0, 1)$value), 0, 1)
2 - pi^2/6

### I( log(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y), 0, 1)$value), 0, 1)
pi^2/12 + 2*log(2) - 2


### I( log(1-x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(1+x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (1 - x*y), 0, 1)$value), 0, 1)


### I( log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (1 + x*y), 0, 1)$value), 0, 1)


### I( (x+y)*log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1+x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
- pi^2/2 + 2*log(2)^2 + 4;


### I( (x+y)*log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1-x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2-x*y+y^2) / (x*y+1) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (x*y+1), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2-x*y+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


###############
###############

############
### ATAN ###

### I( atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y), 0, 1)$value), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2


### I( atan(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 - x*y), 0, 1)$value), 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan)


### I( atan(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x/y) / (x*y+1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x*y), 0, 1)$value), 0, 1)
pi^3 / 48


##############
##############

###########
### Exp ###

### I( log(exp(x) + exp(y)) )
# Maths 505: An extremely captivating double integral
# https://www.youtube.com/watch?v=etR1AI3anpU
# - for polylog2, see file: Integrals.Trig.Tan.R;

integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x) + exp(y)), 0, 1)$value), 0, 1)
- 3/2*pracma::zeta(3) - pi^2/6 - 2*polylog2(-exp(1), 3) + 1/3;


###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(1) - exp(x*y)), 0, 1)$value), 0, 1)
# TODO

###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x*y) - 1), 0, 1)$value), 0, 1)
# TODO


#################

### I( log(tan(x) + tan(y)) )
# Maths 505: A brutal iterated integral!
# https://www.youtube.com/watch?v=UwEgcCZ_SqU

# on [0, pi/4]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(x) + tan(y)), 0, pi/4)$value), 0, pi/4)
7/64 * pracma::zeta(3) + pi^2 * log(2)/16 - pi/4 * Catalan


### I( log(sin(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, pi/4)$value), 0, pi/4)
7/64 * pracma::zeta(3) - pi^2 * log(2)/16


### I( log(tan(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(x+y)), 0, pi/4)$value), 0, pi/4)
0


###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(pi/4*x*y)), 0, 1)$value), 0, 1)
# TODO

