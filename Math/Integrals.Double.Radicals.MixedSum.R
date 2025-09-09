

### Generalisations


### Gen: I( x * ((1 - x*y) / (1 + x*y))^(1/n) )
n = sqrt(7);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma((3*n+1)/(2*n)) - digamma((2*n+1)/(2*n))) * (n+1)/n^2 + (n-2)/(2*n);


### Gen: I( x^2 * ((1 - x*y) / (1 + x*y))^(1/n) )
n = sqrt(5) + 2*sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * (n^2-1)/(3*n^3) - (2*n^2+3*n-6) / (6*n^2);


### Gen: I( x^2*y * ((1 - x*y) / (1 + x*y))^(1/n) )
n = sqrt(13)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * (n^2+3*n+2)/(3*n^3) + 5/6 + 2*(n+1) / n^2;


### Composite:

### Gen: I( x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/n) )
n = sqrt(11) + 3/5;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(1/2 + 1/(2*n)) - digamma(1/(2*n)) - n)^2 / (2*n^2);


##################

### Power 2 (SQRT)

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


### Composite:

### I( ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - pi^2/8 + pi/2 - log(2) - 1;


### I( y * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan - log(2);


### I( ((1-x^2) * (1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x^2) * (1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
2*Catalan - log(2) - 1/2


### I( ((1-x) * (1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1-x) * (1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( ((1+x) * (1 - x*y) / (1 + x*y))^(1/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	((1+x) * (1 - x*y) / (1 + x*y))^(1/2), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


#################

### Higher Powers

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
(digamma(10/6) - digamma(7/6)) * 4/9 + 1/6;

### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(2*digamma(2/3) - digamma(1/6) - digamma(7/6)) * 4/81 - 1/9 + 1/6;
(digamma(10/6) - digamma(1/6)) * 8/81 - 7/18;

### I( x^2*y * ((1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (digamma(10/6) - digamma(1/6)) * 20/(3*3^3) + 31/18;


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


### I( y/x * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y/x * ((1 + x*y) / (1 - x*y))^(1/3) - y/x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(8*digamma(5/6) + digamma(1/3) - 18*digamma(2/3) + 9*digamma(1)) / 18 - 1/6;


### I( x * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(5/6) - digamma(1/3)) * 2/9 + 1/6;


### I( x^2 * ((1 + x*y) / (1 - x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 + x*y) / (1 - x*y))^(1/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(5/6) - digamma(1/3)) * 8/81 + 1/6 + 2/27;


### Radical: 2/3

### I( x * ((1 - x*y) / (1 + x*y))^(2/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(2/3), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(5/6) - digamma(1/3)) * 10/9 + 5/6 - 3;


### Composite:

### I( x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(4/6) - digamma(1/6) - 3)^2 / 18;


### I( y * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	y * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


### I( x^2 * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/3), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
# TODO


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

### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/4) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/4), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(13/8) - digamma(1/8)) * 15/(3*4^3) - 19 / (3*16);

### I( x^2*y * ((1 - x*y) / (1 + x*y))^(1/4) )
n = 4
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * 30/(3*n^3) + 1 + 88 / (3*n^3);


### Composite:

### I( x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/4) )
n = 4;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1-x)/(1+x) * (1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(1/2 + 1/(2*n)) - digamma(1/(2*n)) - n)^2 / (2*n^2);


##############

### Power = 5

### I(x * ((1 - x*y) / (1 + x*y))^(1/5)  )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/5), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(16/10) - digamma(11/10)) * 6/5^2 + 3/10;

### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/5) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/5), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma(16/10) - digamma(1/10)) * 24/(3*5^3) - 59 / (6*25);

### I( x^2*y * ((1 - x*y) / (1 + x*y))^(1/5) )
n = 5
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * (n^2+3*n+2)/(3*n^3) + 985 / (6*n^3);


##############

### Power = 7

### I(x * ((1 - x*y) / (1 + x*y))^(1/7)  )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * ((1 - x*y) / (1 + x*y))^(1/7), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
(digamma(22/14) - digamma(15/14)) * 8/7^2 + 5/14;

### I( x^2 * ((1 - x*y) / (1 + x*y))^(1/7) )
n = 7
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2 * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * (n^2-1)/(3*n^3) - (2*n^2+3*n-6) / (6*n^2);

### I( x^2*y * ((1 - x*y) / (1 + x*y))^(1/7) )
n = 7
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x^2*y * ((1 - x*y) / (1 + x*y))^(1/n), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- (digamma((3*n+1)/(2*n)) - digamma(1/(2*n))) * (n^2+3*n+2)/(3*n^3) + 5/6 + 2*(n + 1) / n^2;

