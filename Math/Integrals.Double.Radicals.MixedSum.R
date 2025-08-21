

### I( ((1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/2 + log(2) - 1;


### I( x * ((1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * 3/4 - 2;


### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(7/4) - digamma(1/4)) / 8 - 1/3;


### I( ((1 + x*y) / (1 - x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 + x*y) / (1 - x*y))^(1/2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/2 - log(2) + 1;


#############

### Power = 3

### I( ((1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO

### I( x * ((1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*digamma(2/3) - digamma(1/6) - digamma(7/6)) * 2/9 - 1/2;


### I( ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(4/3, 1/2) * 2 - 3;


### Reverse:

### I( ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(1/3 - 1/2) + digamma(1/3 + 1/2) - 2*digamma(1/3)) / 9 - 1/2;


#############

### Power = 4

### I( ((1 - x*y) / (1 + x*y))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO




