###################
##
## Trig: TAN & ATAN


####################

### Helper Constants
Catalan = 0.915965594177219015054603514;

# Note:
# - function polylog2 moved to file:
#   Integrals.Polylog.Helper.R;


##################
##################


### Base: I( tan(x) )
integrate(\(x) tan(x) - 1/(pi/2-x), 0, pi/2)
- log(pi/2)


### I( x * tan(x) )

# - see Integrals.Log.Trig for the Base-Integrals:
#   I( log(cos(x)) );

### on [0, pi/2]
integrate(\(x) (pi/2 - x) * tan(x), 0, pi/2)
pi*log(2)/2

### on [0, pi/4]
integrate(\(x) x * tan(x), 0, pi/4)
Catalan / 2 - 1/8 * pi*log(2)

### [by Parts]
integrate(\(x) x^2 * tan(x)^2, 0, pi/4)
- 1/3*(pi/4)^3 + pi^2 / 16 + 1/4 * pi*log(2) - Catalan


### on [0, pi/6]
integrate(\(x) x * tan(x), 0, pi/6)
- pi*log(2)/6 - pi/6*log(cos(pi/6)) + sqrt(3)/(4*12^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) +
	- pracma::psi(1, 1/6) + pracma::psi(1, 5/6) +
	+ pracma::psi(1, 2/6) - pracma::psi(1, 4/6));

###
integrate(\(x) x / tan(x), 0, pi/6)
pi*log(2)/6 + pi/6 * log(sin(pi/6)) + sqrt(3)/(4*36) *
	(pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3));


#####################

### I( x^2 * tan(x) )

### on [0, pi/2]
integrate(\(x) x * (x - pi/2) * tan(x), 0, pi/2)
- 7/8 * pracma::zeta(3)

# [variant]
integrate(\(x) x^2 * tan(x) - (pi/2)^2/(pi/2-x), 0, pi/2)
- 7/8 * pracma::zeta(3) - (pi/2)^2*log(2) - (pi/2)^2*log(pi/2);
- 7/8 * pracma::zeta(3) - (pi/2)^2*log(pi);

#
integrate(\(x) (x - pi/2)*(x + pi/2) * tan(x), 0, pi/2)
- 7/8 * pracma::zeta(3) - 1/4 * pi^2 * log(2)

### on [0, pi/4]
integrate(\(x) x^2 * tan(x), 0, pi/4)
- 21/64 * zeta(3) - 1/32 * pi^2 * log(2) + 1/4 * pi * Catalan

### [by Parts]
integrate(\(x) x^3 * tan(x)^2, 0, pi/4)
(pi/4)^3 - 1/4*(pi/4)^4 + 63/64 * zeta(3) +
	+ 3/32 * pi^2 * log(2) - 3/4 * pi * Catalan;


### Gen:

### I( x^2 * tan(x) )
n = 5; # ODD Integer
integrate(\(x) x^2 * tan(x), 0, pi/n)
id = seq((n-1)/2); idp = 2*id*pi/n;
cs = cos(idp); cs2 = cos(2*idp);
sn = sin(idp); sn2 = sin(2*idp);
- pracma::zeta(3) * (3/4 - 3/4 / n^3) / 2 +
	- (pi/n)^2 * log(2) - (pi/n)^2 * log(cos(pi/n)) +
	- sum((cs - cs2/4) * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))
		) / (4*n^3) +
	+ sum((sn - sn2/2) * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi / (n^3);

### I( x^2 / tan(x) )
n = 5; # ODD Integer
integrate(\(x) x^2 / tan(x), 0, pi/n)
id = seq((n-1)/2); idp = 2*id*pi/n; cs = cos(idp); sn = sin(idp);
- pracma::zeta(3)/2 - pracma::psi(2, 1) / (4*n^3) +
	+ (pi/n)^2 * log(2) + (pi/n)^2 * log(sin(pi/n)) +
	- sum(cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) / (4*n^3) +
	+ sum(sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))) * pi / n^3;


### I( x^3 * tan(x) )

### on [0, pi/2]
integrate(\(x) x * (pi/2 - x)*(pi/2 + x) * tan(x), 0, pi/2)
integrate(\(x) x^2 * log(cos(x)), 0, pi/2)$value +
	+ pracma::zeta(3) * pi + pi^3*log(2)/24;
3/4 * pi * pracma::zeta(3);

# [variant]
integrate(\(x) x^2 * (pi/2 - x) * tan(x), 0, pi/2)
integrate(\(x) x^2 * log(cos(x)), 0, pi/2)$value +
	+ 9/16 * pracma::zeta(3) * pi + pi^3*log(2)/24;
5/16 * pi * pracma::zeta(3);

# Base-variant:
integrate(\(x) - x^3 * tan(x) / 8 + pi^3/64 / (pi/2 - x), 0, pi/2)
(pi^3/8 * log(pi) + 3/4 * pi * pracma::zeta(3)) / 8


### on [0, pi/4]
integrate(\(x) x^3 * tan(x), 0, pi/4)
- 1/128 * pi^3*log(2) + 3/32 * pi^2*Catalan +
	+ 9/256 * pi * pracma::zeta(3) +
	- (pracma::psi(3, 1/4) - pracma::psi(3, 3/4)) / 8 / 4^4;


################
### on [0, pi/4]

### I( x * tan(x) )
integrate(\(x) x * tan(x), 0, pi/4)
integrate(\(x) log(cos(x)), 0, pi/4)$value +
	- (pi/4) * log(cos(pi/4));
- 1/8 * pi*log(2) + Catalan/2;

### I( x^2 * tan(x) )
integrate(\(x) x^2 * tan(x), 0, pi/4)
integrate(\(x) 2*x * log(cos(x)), 0, pi/4)$value +
	- (pi/4)^2 * log(cos(pi/4));
- 1/32 * pi^2 * log(2) + pi * Catalan / 4 - 21/64 * pracma::zeta(3);

### I( x^3 * tan(x) )
integrate(\(x) x^3 * tan(x), 0, pi/4)
integrate(\(x) 3*x^2 * log(cos(x)), 0, pi/4)$value +
	- (pi/4)^3 * log(cos(pi/4));
- 1/128 * pi^3*log(2) + 3/32 * pi^2*Catalan +
	+ 9/256 * pi * pracma::zeta(3) +
	- (pracma::psi(3, 1/4) - pracma::psi(3, 3/4)) / 8 / 4^4;


##############

###
integrate(\(x) x / tan(x), 0, pi/4)
pi*log(2)/8 + Catalan/2


###
integrate(\(x) tan(x) / x, 0, pi/4)
# TODO: ???


###
integrate(\(x) x / atan(x), 0, 1)
integrate(\(x) (tan(x) + tan(x)^3)/x, 0, pi/4)

# TODO: ???

###
integrate(\(x) (tan(x) - tan(x)^3) / x, 0, pi/4)
# TODO: ???


###
integrate(\(x) tan(x)^2 / x, 0, pi/4)
integrate(\(x) 1/atan(x) - 1/atan(x) / (x^2+1), 0, 1)
integrate(\(x) 1/atan(x) - 1/x, 0, 1)$value - log(pi/4)
# TODO

# I( 1 / ((x^2+1)*atan(x)) )
integrate(\(x) 1/atan(x) / (x^2+1) - 1/x, 0, 1)
log(pi/4)


############

# Maths 505: A MONSTER INTEGRAL!!!
# int from 0 to infty arctan(x)/(x(x+1)(x^2+1))
# https://www.youtube.com/watch?v=u-FKjn_83l8
# Note: Intermediate integral
# - see Integrals.Log.Trig.R for solutions to LOG( Trig );

integrate(\(x) x / (tan(x) + 1), 0, pi/2)
pi^2/16 + pi*log(2)/8 - Catalan/2

# Derivation:
integrate(\(x) 2 * x / (tan(x) + 1), 0, pi/2)
integrate(\(x) - log(tan(x) + 1), 0, pi/2)$value +
	+ pi*log(2)/2 + pi^2 / 8;


#######################
#######################

############
### ATAN ###

### ATAN: Powers
# moved to file: Integrals.Trig.Tan.Atan.R;
### ATAN( Trig )
# moved to file: Integrals.Trig.Tan.Atan.Trig.R;
### ATAN: Fractions of type ATAN() / Poly
# moved to file: Integrals.Trig.Tan.Atan.Fractions.R;
### ATAN * LOG
# moved to file: Integrals.Trig.Tan.Atan.LogCombi.R;


### I( atan(x^n) )
n = sqrt(5)
integrate(\(x) atan(x^n), 0, 1)
pi/4 - integrate(\(x) n*x^n / (x^(2*n) + 1), 0, 1)$value
pi/4 - (digamma(((n+1)/(2*n) + 1)/2) - digamma((n+1)/(4*n))) / 4


### I( x^p * atan(x^n) )
n = sqrt(5); p = sqrt(7);
integrate(\(x) x^p * atan(x^n), 0, 1)
pi/(4*(p+1)) - integrate(\(x) n/(p+1) * x^(n+p) / (x^(2*n) + 1), 0, 1)$value
pi/(4*(p+1)) - (digamma(((n+p+1)/(2*n) + 1)/2) - digamma((n+p+1)/(4*n))) / (4*(p+1))
pi/(4*(p+1)) - (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1))

### Simple:
p = sqrt(3)
integrate(\(x) x^p * atan(x), 0, 1)
(digamma((p+2)/4) - digamma(p/4 + 1) + pi) / (4*(p+1))


### Lim: p -> -1
n = 1/sqrt(3)
integrate(\(x) atan(x^n) / x, 0, 1)
Catalan / n;


###################
###################
### ATAN: Fractions

### Atan of Fractions

# - see also file:
#   Integrals.Trig.Formulas.R;


### I( atan((x^n - 1)/(x^n + 1)) )
integrate(\(x) atan((x^3 - 1)/(x^3 + 1)), 0, 1)
integrate(\(x) -3*x^3 / (x^6 + 1), 0, 1)
(digamma(4/12) - digamma(1/2 + 4/12)) / 4
(digamma(1/3) - digamma(5/6)) / 4


### Gen:
n = sqrt(5)
integrate(\(x) atan((x^n - 1)/(x^n + 1)), 0, 1)
integrate(\(x) -n*x^n / (x^(2*n) + 1), 0, 1)
(digamma(1/4 + 1/(4*n)) - digamma(3/4 + 1/(4*n))) / 4

