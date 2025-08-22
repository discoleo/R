

### Generalisations


### Gen: I( x * ((1 - x*y) / (1 + x*y))^(1/n) )
n = sqrt(7);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma((3*n+1)/(2*n)) - digamma((2*n+1)/(2*n))) * (n+1)/n^2 + (n-2)/(2*n);


###################

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

### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*digamma(2/3) - digamma(1/6) - digamma(7/6)) * 4/81 - 1/9 + 1/6;


### I( ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x^0.5 * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^0.5 * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( x * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(4/3, 1/2) * 2 - 3;


### I( x*y * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( x^2 * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
beta(4/3, 1/2) / 4 - 3/16;


### I( x^2*y * ((1 - x*y) / (1 + x*y)^2)^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y)^2)^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- beta(4/3, 1/2) * 3/2 + 21/8;
beta(4/3, 1/2) * betan(-2/3, 1) + 21/8; # ?


### Reverse:

### I( ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(1/3 - 1/2) + digamma(1/3 + 1/2) - 2*digamma(1/3)) / 9 - 1/2;


### I( x^2 * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(1/3 - 1/2) + digamma(1/3 + 1/2) - 2*digamma(1/3)) * 4/81 + 1/6 - 2/9;


#############

### Power = 4

### I( ((1 - x*y) / (1 + x*y))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1 - x*y) / (1 + x*y))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( x * ((1 - x*y) / (1 + x*y))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/4), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(13/8) - digamma(9/8)) * 5/16 + 1/4;


##############

### Power = 5

### I(x * ((1 - x*y) / (1 + x*y))^(1/5)  )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/5), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(16/10) - digamma(11/10)) * 6/5^2 + 3/10;


##############

### Power = 7

### I(x * ((1 - x*y) / (1 + x*y))^(1/7)  )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/7), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(22/14) - digamma(15/14)) * 8/7^2 + 5/14;

