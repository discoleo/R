############
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.1.

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
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1-x*y), 0, 1)$value), 0, 1);
pi^2 / 6


### I( (x+y)*log(x+y) / (x*y + 1) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1+x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
- pi^2/2 + 2*log(2)^2 + 4;

