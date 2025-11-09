#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Type: Atan
##
## v.0.2a


### Type: Atan

### Examples:
# I( atan(x/y) / (x+y) )
# I( atan(x*y) / (1 + x*y) )
# I( atan(x*y) / (1 - x*y) ) [not yet]
# I( atan(sqrt(tan(x)^2 + tan(y)^2)) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

############
### ATAN ###

# Note: (x+y) does not behave as nicely with ATAN;

### I( atan((x+y)/2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2), 0, 1)$value), 0, 1)
3/2*pi + 2*log(5) - 6*log(2) - 3*atan(2)


### I( atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y), 0, 1)$value), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2


### I( atan(x/y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y), 0, 1)$value), 0, 1)
pi/4;


### I( atan(2*x / (x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(2*x / (x+y)), 0, 1)$value), 0, 1)
integrate(\(x) atan(2*x/(x+1)) +
	sapply(x, \(y) integrate(\(x) 2*x*y / (4*x^2 + (x+y)^2), 0, 1)$value), 0, 1);
37/100 * pi + 6/5 * log(2) - log(5)/2 + 1/50 * (2*atan(3) - 23*atan(2));


# Helper:
integrate(\(x) atan(2*x/(x+1)), 0, 1)
integrate(\(x) atan(2 - 2/(x+1)), 0, 1)
pi/4 - (log(8)/5 - 1/5 * (atan(3) - atan(1/2)));

#
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / (4*x^2 + (x+y)^2), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) y*(x-y) / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) -1/4 + sapply(x, \(y)
	integrate(\(x) (x*y + x^2/4) / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)) + sapply(x, \(y)
	integrate(\(x) - y^2 / (4*y^2 + x^2), y, y+1)$value), 0, 1)
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)) +
	- x/2 * (atan((x+1)/(2*x)) - atan(1/2)), 0, 1)
2/25 * (atan(1/2) - atan(3)) - 11/200 * pi + 1/4*atan(1/2) +
	+ 3/10 * log(8) - log(5)/4;

#
integrate(\(x) x/2 * (log(5*x^2 + 2*x + 1) - log(5*x^2)), 0, 1)
integrate(\(x) -4/25 / (5*x^2 + 2*x + 1), 0, 1)$value +
	+ 7/25*log(8) - log(5)/4 + 1/10;
2/25 * (atan(1/2) - atan(3)) + 7/25*log(8) - log(5)/4 + 1/10;


### I( atan(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 - x*y), 0, 1)$value), 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan)


### I( atan(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	- x*y / ((x*y + 1)^2 + 1), 0, 1)$value), 0, 1)$value +
	+ 2*atan(2) - 1/2*log(5) + 1/2*log(2) - atan(1);
# TODO

# alternative:
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y) + atan(1-x*y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2*y^2/2), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2/2) / y, 0, y)$value), 0, 1)
#
integrate(\(x) atan(x^2/2), 0, 1)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x) -4*x^2 / (x^4 + 4) / y, 0, y)$value), 0, 1)$value;
#
integrate(\(x) - 2*sqrt(2) * Re(atan(x / sqrt(2i)) / sqrt(1i)) / x, 0, 1)$value +
	+ log(5)/2 - atan(3/4);
#
integrate(\(x) (pi/4 - atan(x+1)) / x, 0, 1)$value +
integrate(\(x) atan(3*tan(x)), 0, pi/4)$value +
	+ (pi*log(2)/8 - Catalan) + pi^2 / 32 + log(5)/2 - atan(3/4);
# TODO


### I( atan( (x+y)/(1+x*y) ) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan((x+y)/(1+x*y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/2 + 2*log(2) - 2*Catalan;


###################

### Type: Atan(x/y)

### I( atan(x/y) / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan + pi/4 + log(2)/2;

### I( atan(x/y) * y / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/2 + 1/4;

### I( atan(x/y) * y^2 / x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y^2 / x, 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
Catalan/3 + (pi/4 - log(2)/2) / 9 + 1/18;


### I( atan(x/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) * abs(x-y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
# trick for numerical accuracy:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y), y, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value;
pi*log(2)/2 - pi/4;


### I( atan(x/y) / (1+x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x+y), 0, 1)$value), 0, 1)
(log(3) * 3/4 - log(2)) * pi

### I( atan(x/y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (2 - x - y), 0, 1)$value), 0, 1)
pi*log(2)/2


### I( atan(x/y) * x / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * x / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
pi/4 +  pi * log(2)/8 - log(2)/2;

### I( atan(x/y) * y / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * y / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
3/8*pi*log(2) - pi/4 + log(2)/2;

### I( atan(x/y) * x*y / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * x*y / (x+y)^2, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-13)
pi*log(2)/4 - pi/8;


### I( atan(x/y) * abs(x-y) / (x+y)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, 0, 1, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)
# trick for numerical accuracy:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, 0, y, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * abs(x-y) / (x+y)^2, y, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value;
pi*(1 - log(2)) / 2;


### I( atan(x/y)^2 / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y)^2 / (x+y), 0, 1, rel.tol=1E-11)$value), 0, 1, rel.tol=1E-12)
pi^2*log(2) * 3/16 + pi*Catalan/2 - 21/16*pracma::zeta(3);


###################

### Type: Atan(X*Y)


### I( atan(x*y) * y/x )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * y/x, 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
Catalan/2 - pi/8 + 1/4;


### I( (atan(x*y) - atan(y)) / (1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x*y) - atan(y)) / (1-x), 0, 1)$value), 0, 1)
pi^2/32 + pi*log(2)/8 - Catalan - log(2)^2/8;


### I( atan(x*y) / (x*(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) / (x*(x+y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
(pi/4)^3;


### I( atan(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y), 0, 1)$value), 0, 1)
# Derivation: z = x/y;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(y^2*x) / (x+1), 0, 1/y)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * atan(y^2*x) / (x+1), 0, 1)$value), 0, 1)
# Factor: *2;
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x^2*y) / (y+1), 0, 1)$value), 0, 1)
integrate(\(x) atan(x)/(x+1) - 2/(x*(x+1)) * sapply(x, \(y)
	integrate(\(x) x^2 / (x^4 + 1/y^2), 0, 1)$value), 0, 1);
isq = exp(-1i * pi/4); isq2 = exp(1i * pi/4);
integrate(\(x) atan(x)/(x+1) +
	- sqrt(x)/(x*(x+1)) * Re(atan(sqrt(x)/isq)/isq + atan(sqrt(x)/isq2)/isq2), 0, 1)
integrate(\(x) -4 * Re(atan(x/isq)/isq) / (x^2+1), 0, 1)$value + pi * log(2)/8;
integrate(\(x) 4 * Re(1/(x^2 + 1i)) * atan(x), 0, 1)$value +
	(digamma(3/8) - digamma(7/8)) * pi / 8 + pi * log(2)/8;
integrate(\(x) 4 * x^2/(x^4 + 1) * atan(x), 0, 1)$value +
	(digamma(3/8) - digamma(7/8)) * pi / 8 + pi * log(2)/8;
# TODO

# Helper:
isq = exp(-1i * pi/4); isq2 = exp(1i * pi/4);
integrate(\(x) Re(atan(x/isq)/isq + atan(x/isq2)/isq2) / (x^2+1), 0, 1)
integrate(\(x) Re(atan(x/isq)/isq) / (x^2+1) * 2, 0, 1)
integrate(\(x) (x^2 - 1) * (log(1-x) - Re(log(x + 1i))) / (x^4 + 1), 0, 1)$value +
	- (digamma(7/8) - digamma(3/8)) * pi / 32 +
	- (pracma::psi(1, 3/8) - pracma::psi(1, 7/8)) / 64;
# see Integrals.Log.Fractions.P4.R;
integrate(\(x) (1 - x^2) * Re(log(x + 1i)) / (x^4 + 1), 0, 1)$value +
	- (pracma::psi(1, 5/8) - pracma::psi(1, 1/8)) / 64 +
	+ (pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 32 +
	- (digamma(5/8) - digamma(1/8)) * log(2) / 32 +
	+ (digamma(7/8) - digamma(3/8)) * log(2) / 32 +
	- (digamma(5/8) - digamma(1/8)) * pi / 64 +
	- (digamma(7/8) - digamma(3/8)) * pi * 3 / 64;

# TODO


### I( atan(x*y) * x / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * x / (x+y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
- pi^2 / 96 + pi/8 - log(2)/4;


### I( atan(x*y) / (x + y + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y+1), 0, 1)$value), 0, 1)
# TODO


### I( atan(x*y) / (2 - x - y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (2 - x - y), 0, 1)$value), 0, 1)
#
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * atan(x*(y-x)) / (2 - y), 0, y/2)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) -2*atan(x*(y-x)) / (2 - y), 0, y-1)$value +
	integrate(\(x) +2*atan(x*(y-x)) / (2 - y), 0, y/2)$value ), 1, 2)$value;
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) -2*atan(x*(y-x)) / (2 - y), 0, y-1)$value +
	integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value ), 1, 2)$value;
#
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*(y-x)) / (2 - y), 0, y)$value), 0, 1)$value +
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 2*atan(x*(y-x)) / (2 - y), 0, 1)$value +
	integrate(\(x) - atan(x*(y-x)) / (2 - y), 0, y)$value ), 1, 2)$value;
# TODO


### I( atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 + x*y), 0, 1)$value), 0, 1)
((digamma(1/4) + Euler) * pracma::psi(1, 1/4) - pracma::psi(2, 1/4)/2) / 16 +
	- 7/4 * pracma::zeta(3) - (pi/4)^3 + 3/16 * pi^2 * log(2) +
	+ Catalan * pi/4 + Catalan * log(2);
(pi/4)^3 - Catalan * log(2)/2;


# see also: I( atan(x*y) / (1 - x*y) )
integrate(\(x) -1/3 * atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) -1/3 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	integrate(\(x) 1/3 * log(x) * log(x^2+1) / (x^2+1), 0, 1)$value;
integrate(\(x) -1/3 * atan(x) * log(1-x) / x, 0, 1)$value +
integrate(\(x) +1/3 * log(x) * log(1-x) / (x^2+1), 0, 1, rel.tol=1E-12)$value - Catalan * log(2)/3;


### I( atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 - x*y), 0, 1)$value), 0, 1)
integrate(\(x) - atan(x) * log(1-x) / x, 0, 1)$value +
	integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1)$value;
integrate(\(x) - atan(x) * log(x) / (1-x), 0, 1, rel.tol=1E-12); # alternative
#
integrate(\(x) -2 * log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	+ pi^3 * 3/64 - Catalan * log(2)/2;
# TODO


### I( atan(x*y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 48 - pi/4 + log(2)/2 + 3/32 * pracma::zeta(3);

### I( atan(x*y) * log(1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * log(1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
-23/64 * pracma::zeta(3) + Catalan * pi / 4 - Catalan +
	+ 5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4;

### I( atan(x*y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x*y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) - log(1-x) * log(x) / (x^2+1), 0, 1)$value +
	+ pi*log(2)/4 + (pi/4)^3 + (pi/4)^2 - Catalan * log(2)/2 - Catalan;
# TODO


### Div

### I( atan(x/y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x*y), 0, 1)$value), 0, 1)
pi^3 / 48

### I( atan(x/y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1-x*y), 0, 1)$value), 0, 1)
pi^3 / 24;

### I( atan(x/y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) / (2 - (x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (2-x-y), 0, 1)$value), 0, 1)
pi*log(2)/2

### I( atan(x/y) / (x+y+1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (x+y+1), 0, 1)$value), 0, 1)
(log(3) * 3/4 - log(2)) * pi

### I( atan(x/y) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	(atan(x/y) - pi/4) / (x-y), 0, y, rel.tol=1E-6)$value), 0, 1, rel.tol=1E-6)
Catalan - pi*log(2)/8;
# Note: on [y, 1] = identical;


### Div: Higher

### I( x*y * atan(x) / (x^2+y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x) / (x^2+y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi*log(2)/8 - pi/8 - log(2)/4 + Catalan/2;

### I( x*y * atan(x) / (x^2*y^2 + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x*y * atan(x) / (x^2*y^2 + 1), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) log(x) * log(1-x) / (x^2+1), 0, 1)$value +
	- (pi/4)^3 + Catalan * log(2)/2;
# TODO


### Prod: w. LOG

### I( atan(x/y) * log(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 96 - pi/4 + log(2)/4;

### I( atan(x/y) * log(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- pi^2 / 96 - pi/4 - log(2)/4;

### I( atan(x/y) * log(1-x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1-x), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
pi^2 / 12 + pi/8 * log(2) - pi * 3/8 - Catalan;

### I( atan(x/y) * log(1-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * log(1-y), 0, 1, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
- (pi^2 / 12 + pi/8 * log(2) + pi/8 - Catalan);


### I( atan(x/y) * atan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) * atan(x), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^2/16 + 3/8 * pi*log(2) - pi/8 - log(2)/4 - Catalan/2;


### Div: SQRT(...)

### I( atan(x/y) / sqrt(1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * (1-log(2));


### I( atan(x/y) / sqrt(1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi * (log(2) - asinh(1)) + (digamma(5/8) - digamma(3/8));


### I( atan(x/y) / sqrt(x^2 + y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / sqrt(x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi/2 * asinh(1);


### ATAN(x+y)

### I( atan(x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y), 0, 1)$value), 0, 1)
3/2*atan(2) - log(5) + log(2)

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

### I( atan(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1+x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan((x+y)/2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan((x+y)/2) / (1+x*y), 0, 1)$value), 0, 1)
# TODO


### Atan( Fraction )

### I( atan(x/y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(x/y) / (1 - x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
pi^3 / 24


################

### Atan( TRIG )

### I( atan(tan(x)*tan(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(tan(x)*tan(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;

### I( atan(tan(x) / tan(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(tan(x) / tan(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;

### I( atan(sqrt(tan(x) * tan(y))) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(tan(x) * tan(y))), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3/16;


### I( atan(sqrt(tan(x)^2 + tan(y)^2)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(tan(x)^2 + tan(y)^2)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(digamma(7/8) - digamma(3/8)) * (digamma(5/8) - digamma(1/8)) * pi / 16;


### I( atan(sin(x) * sin(y)) / sin(x) )
# Maths 505: A beautiful double integral
# https://www.youtube.com/watch?v=tsxLAsiLyq4
# Note: Series expansion atan(x)
# => I(sin(x)^(2*k)) * I(sin(y)^(2*k+1)) => Beta;

integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan * pi / 2;


### I( atan(sqrt(sin(x) * sin(y))) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sqrt(sin(x) * sin(y))) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan * pi;


### I( atan(sin(x) * sin(y)) * cos(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * cos(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) - (digamma(5/8) - digamma(3/8)) / 2;
# Note: asinh(1) == log(tan(3*pi/8));

### I( atan(sin(x) * sin(y)) * sin(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) - (digamma(5/8) - digamma(3/8)) / 2;


### I( atan(sin(x) * sin(y)) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) * pi/16;


### I( atan(sin(x) * sin(y)) / (sin(x) * sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / (sin(x) * sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-13)
# TODO

#
id = 0:100;
sum((-1)^id / (2*id+1) * beta(id+1/2, 1/2)^2) / 4;

# Note:
id = seq(0, 20000)
sum((-1)^id / (id+1) * beta(id+1/2, 1/2) * beta(id+3/2, 1/2)) / 4;
((beta(1/4, 1/4) + 2*beta(3/4, 3/4)) * sqrt(2) / 4 - pi) * pi/2;
#
sum((-1)^id / (id+1) / (2*id+1) * beta(id+1/2, 1/2) * beta(id+3/2, 1/2)) / 4
(pi - beta(3/4, 3/4) * sqrt(2)) * pi/2;


### I( atan(sin(x) * sin(y)) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * cos(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * cos(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
Catalan - pi/4 + log(2)/2;


### I( atan(sin(x) * sin(y)) / tan(x) / tan(y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / tan(x) / tan(y), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi^3 / 32;


### I( atan(sin(x) * sin(y)) / sin(x)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sin(x)^2 - sin(y)/x, 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * sin(y) / sin(x)^2 )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(y) / sin(x)^2 - sin(y)^2/x, 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/4 * (1/2 - log(pi/2));


### Div: Sum

### I( atan(sin(x) * sin(y)) / sqrt(sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) / sqrt(sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
beta(1/8, 3/4)^2 * beta(3/8, 3/4) *  beta(1/8, 1+1/4) * 2 / 8^4;
beta(1/8, 3/4)^2 * beta(1/4, 3/4) * 4/3 / 8^3;


### I( atan(sin(x) * sin(y)) * sin(x) / (sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(x) / (sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
# ==
integrate(\(x) pi/8 * sapply(x, \(y) integrate(\(x)
	sqrt( x / (1-x^2) / (1 - x^2*y^2) ), 0, 1, rel.tol=1E-10)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( atan(sin(x) * sin(y)) * sin(x)/sin(y) / (sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) * sin(y)) * sin(x)/sin(y) / (sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
2.309240197354 / 2;
# see I( atan(sin(x) * sin(y)) / (sin(x) * sin(y)) )
# TODO


###########
### DIV ###

### I( atan(sin(x) / sin(y)) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
integrate(\(x) pi/2 * atan(sqrt(sin(x))) / sin(x), 0, pi/2, rel.tol=1E-13)
# TODO


### I( atan(sin(x) / sin(y)) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(pi + 2*log(2)) * pi / 8;


### I( atan(sin(x) / sin(y)) * cos(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * cos(y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
((digamma(5/8) - digamma(1/8))^2 - (digamma(7/8) - digamma(3/8))^2) / 16;
log(sqrt(2) + 1) * pi;

### I( atan(sin(x) / sin(y)) * sin(2*y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(2*y) / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * asinh(1) + (digamma(5/8) - digamma(3/8)) / 2;


### I( atan(sin(x) / sin(y)) * sin(y)^2 / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y)^2 / sin(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
3.668202541275/2 - (beta(1/4, 1/4) - 2*beta(3/4, 3/4)) * sqrt(2) * pi/32;
# TODO: I( atan(sin(x) / sin(y)) / sin(x) )


### I( atan(sin(x) / sin(y)) * sin(y) / tan(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) * sin(y) / tan(x), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/2 * (asinh(1) - log(2) + 1) - (digamma(5/8) - digamma(3/8)) / 2;

# Note:
# also = pi/2 - 1/2 * I( atan(x/y) / sqrt(1+x*y) ) on [0,1]^2;


### I( atan(sin(x) / sin(y)) / sqrt(sin(x)^2 + sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x) / sin(y)) / sqrt(sin(x)^2 + sin(y)^2), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
gamma(1/8)^2 / gamma(7/8)^2 * gamma(3/4)^3 * gamma(1/4) / 8^2 / 2;
beta(1/8, 3/4)^2 * beta(1/4, 3/4) / 8^2 / 2;


### Power = 2

### I( atan(sin(x)^2 / sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^2 / sin(y)^2), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^3 / 16;

### I( atan(sin(x)^3 / sin(y)^3) )
p = 3; # any meaningful value;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^p / sin(y)^p), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^3 / 16;


### I( atan(sin(x)^2 / sin(y)^2) * sin(y) / sin(x) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(sin(x)^2 / sin(y)^2) * sin(y) / sin(x), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
# see I( (x^2 - 1) * Im(log(x + 1i)) / (x^4 + 1) ) on [0,1];
(pracma::psi(1, 1/8) - pracma::psi(1, 5/8)) / 32 - Catalan / 2 +
- (digamma(7/8) - digamma(3/8)) * pi * 3/32 +
- (digamma(5/8) - digamma(1/8)) * pi / 32 +
- (digamma(7/8) - digamma(3/8)) * log(2) / 16 +
- (digamma(5/8) - digamma(1/8)) * log(2) / 16 +
+ sqrt(2)*pi/8 * (pi/2 + log(2)) + pi*log(2) / 8;


#############
### Varia ###

### I( log(beta(x, y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(beta(x, y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
log(2*pi)/2 + 3/4;


### I( atan(beta(x, y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	atan(beta(x, y)), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO

