#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: ATAN( SUM )
##
## v.0.1a

# Note:
# - (x+y) does not behave as nicely with ATAN;
# - Difficult to solve;

####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

### ATAN(x+y)

### I( atan(x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y), 0, 1)$value), 0, 1)
3/2*atan(2) - log(5) + log(2);

### I( atan((x+y) / 2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y) / 2), 0, 1)$value), 0, 1)
3*atan(1/2) + log(5) - 6*log(2) + log(5);
3/2*pi + 2*log(5) - 6*log(2) - 3*atan(2)


### Div:

### I( atan(x+y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (x + y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, y+1)$value), 0, 1)
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, 1)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 1, y+1)$value), 0, 1)$value;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 1, y+1)$value), 0, 1)$value +
	+ pi/4 - log(2)/2;
# alternative:
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, y, y+1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x) / x, 0, y)$value), 0, 2)$value +
	- (Catalan - pi/4 + log(2)/2) * 2;
integrate(\(x) atan(x) / x * (2-x), 0, 2)$value - (Catalan - pi/4 + log(2)/2) * 2;
integrate(\(x) 2 * atan(x) / x, 0, 2)$value +
	- (Catalan - pi/4 + log(2)/2) * 2 - 2*atan(2) + log(5)/2;
# TODO

# Note:
integrate(\(x) atan(x) / x, 0, 2)
pi*log(2)/2 + 0.48722235829452235711;
# sum( (-1)^n / 2^(2*n+1) / (2*n+1)^2) );


### I( atan(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( x * atan((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x * atan((x+y)/2) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( x * atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x * atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

