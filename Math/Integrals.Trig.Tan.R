###################
##
## Trig: TAN & ATAN


####################

### Helper Constants
Catalan = 0.915965594177219015054603514;

### Helper Functions

polylog2 = function (z, n = 2) {
	# Bug in pracma::polylog;
	stopifnot(is.numeric(n));
	if(is.complex(z) && n == 2) {
		if(Im(z) == 0) { z = Re(z); }
		else {
			warning("Complex not yet implemented!");
			if(abs(z) != 1) warning("Result will be incorrect!");
			# TODO: only Real-part:
			y = - (pi^2/6 + Re(log(-z)^2)/2) / 2;
			return(y);
		}
	}
	if(z == 1) {
		return(pracma::zeta(n));
	}
	if(z < 0) {
		if(n == 2) {
			y = polylog2(z^2, n=n)/2 - polylog2(-z, n=n);
		} else {
			y = polylog2(z^2, n=n) * 2^(1-n) - polylog2(-z, n=n);
		}
		return(y);
	}
	if(z > 1) {
		if(n == 2) {
			y = - pi^2/6 - log(0i - z)^2/2 - polylog2(1/z, n=n);
		} else if(n == 3) {
			x = log(z) / (2*pi);
			y = 4*pi^3/3 * (- x^3 + 3i/2*x^2 + x/2) + polylog2(1/z, n=n);
		} else stop("Not yet implemented!");
		return(y);
	} else if(z >= 0.55) {
		# pracma::polylog FAILS!
		if(n == 2) {
			y = pi^2/6 - log(z)*log(1-z) - pracma::polylog(1 - z, n=n);
		} else if(n == 3) {
			# see Ref;
			y = log(1-z)^3 / 6 - log(z)*log(1-z)^2/2 + pi^2/6 * log(1-z) +
				+ pracma::zeta(3) - polylog(1-z, 3) - polylog2(z/(z-1), n=3);
		} else stop("Not yet implemented!");
		return(y);
	}
	return(pracma::polylog(z, n=n));
}
# Ref:
# https://math.stackexchange.com/questions/942796/compute-polylog-of-order-3-at-frac12


# Test: 1 / (0.6 / (1-0.6)) == 2/3 > 0.55;
# - requires 2 iterations;
li3_06 = 0.6560025136329806832346611928113322291802755981380034005592;
polylog2(0.6, 3) - li3_06;


##################
##################

### I( x * tan(x) )

### on [0, pi/2]
integrate(\(x) (pi/2 - x) * tan(x), 0, pi/2)
pi*log(2)/2

### on [0, pi/4]
integrate(\(x) x * tan(x), 0, pi/4)
Catalan / 2 - 1/8 * pi*log(2)

### [by Parts]
integrate(\(x) x^2 * tan(x)^2, 0, pi/4)
- 1/3*(pi/4)^3 + pi^2 / 16 + 1/4 * pi*log(2) - Catalan


### I( x^2 * tan(x) )

### on [0, pi/2]
integrate(\(x) x * (x - pi/2) * tan(x), 0, pi/2)
- 7/8 * pracma::zeta(3)

# [variant]
integrate(\(x) (x - pi/2)*(x + pi/2) * tan(x), 0, pi/2)
- 7/8 * pracma::zeta(3) - 1/4 * pi^2 * log(2)

### on [0, pi/4]
integrate(\(x) x^2 * tan(x), 0, pi/4)
- 21/64 * zeta(3) - 1/32 * pi^2 * log(2) + 1/4 * pi * Catalan

### [by Parts]
integrate(\(x) x^3 * tan(x)^2, 0, pi/4)
(pi/4)^3 - 1/4*(pi/4)^4 + 63/64 * zeta(3) +
	+ 3/32 * pi^2 * log(2) - 3/4 * pi * Catalan;


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

### I( x^p * log(x) * atan(x) )
# Maths 505: A cool result: int (0,1) arctan(x)ln(x)
# https://www.youtube.com/watch?v=zL169U85SZc
# Note: uses series expansion of atan,
#   but method based on D(digamma) is more powerful;


###
p = sqrt(3)
integrate(\(x) x^p * log(x) * atan(x), 0, 1)
(pracma::psi(1, (p+2)/4) - pracma::psi(1, p/4 + 1)) / (16*(p+1)) +
	- (digamma((p+2)/4) - digamma(p/4 + 1) + pi) / (4*(p+1)^2);

### Example:
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2 / 48 - pi/4 + log(2)/2


### Varia:
p = sqrt(3)
integrate(\(x) x^p * log(x) / (x^2 + 1), 0, 1)
(pracma::psi(1, (p+1)/4 + 1/2) - pracma::psi(1, (p+1)/4)) / 16

###
p = sqrt(3)
integrate(\(x) x^p * log(x)^2 / (x^2 + 1), 0, 1)
(pracma::psi(2, (p+1)/4 + 1/2) - pracma::psi(2, (p+1)/4)) / 64

###
p = sqrt(3)
integrate(\(x) x^p / (x^2 + 1)^2, 0, 1)
(digamma((p+1)/4) - digamma((p-1)/4)) * (p-1)/8 - 1/4


### Full:

### I( x^p * atan(x^n) * log(x) )
n = sqrt(5); p = sqrt(3);
# works only for n > 0;
integrate(\(x) x^p * atan(x^n) * log(x), 0, 1)
- (pracma::psi(1, ((p+1)/n + 3)/4) - pracma::psi(1, ((p+1)/n + 1)/4)) / (16*n*(p+1)) +
	+ (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1)^2) +
	- pi/(4*(p+1)^2);


### Special Case:
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2 / 48 - pi/4 + log(2)/2


integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32


####################

### Atan * Log

### I( atan(x) * log(x) )
integrate(\(x) atan(x) * log(x), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2

### I( atan(x) * log(1-x) )
integrate(\(x) atan(x) * log(1-x), 0, 1)
integrate(\(x) ((1-x)*log(1-x) + x) / (x^2+1) - pi/4, 0, 1)
integrate(\(x) (1-x)*log(1-x) / (x^2+1) - pi/4 + log(2)/2, 0, 1)
5/96 * pi^2 - log(2)^2 / 8 + pi*log(2)/8 - Catalan - pi/4 + log(2)/2;

### I( atan(x) / (1-x) )
integrate(\(x) atan(x) / (1-x) - pi/4/(1-x), 0, 1)
pi*log(2)/8 - Catalan

### I( log(x^2+1) * log(1-x) )
integrate(\(x) log(x^2+1) * log(1-x), 0, 1)
integrate(\(x) (1-x)*log(1-x) * 2*x/(x^2+1) + 2 - pi/2 - log(2), 0, 1)
integrate(\(x) 2*(x+1)*log(1-x) / (x^2+1) + 4 - pi/2 - log(2), 0, 1)
-5/48 * pi^2 + log(2)^2 / 4 + pi*log(2)/4 - 2*Catalan - pi/2 - log(2) + 4;

# Complex:
integrate(\(x) Re(log(x + 1i)) * log(1-x), 0, 1)
-5/96 * pi^2 + log(2)^2 / 8 + pi*log(2)/8 - Catalan - pi/4 - log(2)/2 + 2;
#
integrate(\(x) Im(log(x + 1i)) * log(1-x), 0, 1)
- (5/96 * pi^2 - log(2)^2 / 8 + pi*log(2)/8 - Catalan + pi/4 + log(2)/2);


###
integrate(\(x) atan(x)^2 * log(x), 0, 1)
integrate(\(x) - 2*x * atan(x) * log(x) / (x^2+1), 0, 1)$value +
	- (pi/4)^2 + Catalan - 1/4 * pi * log(2)

# TODO

###
integrate(\(x) atan(x)^2 * log(1 - x), 0, 1)
# TODO

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


### Simple Fractions

###
integrate(\(x) atan(x) / (x+1), 0, 1)
pi * log(2)/8

###
integrate(\(x) atan(x) / (x+1)^2, 0, 1)
log(2)/4


###
integrate(\(x) atan(x) / (x+1)^3, 0, 1)
- pi/32 + log(2)/8 + 1/8


###
integrate(\(x) (atan(x) - pi/4) / (x - 1), 0, 1)
- pi*log(2)/8 + Catalan

###
integrate(\(x) (atan(x) - pi/4) / (x - 1)^2 + 1/(2*x), 0, 1)
integrate(\(x) 1 / ((x^2 + 1)*(x - 1)) + 1/(2*x) + pi/4 - 1/2, 0, 1)
-log(2)/4 + pi/8 - 1/2


### Other

### I( (atan(x) - pi/4) / (x^2 - 1) )
integrate(\(x) (atan(x) - pi/4) / (x^2 - 1), 0, 1)
Catalan / 2

# Extra-Component: pi/8*log(2);


#######################

### ATAN: Fractions

### I( atan(tan(x)^2) )
# Maths 505: An outrageous journey of integration: int 0 to pi/4 arctan(cot^2(x))
# https://www.youtube.com/watch?v=VUGlU_dSgPY

### I( atan(x^2) / (x^2 + 1) )
integrate(\(x) atan(tan(x)^2), 0, pi/4)
integrate(\(x) atan(x^2) / (x^2 + 1), 0, 1)
pi^2/16 - (digamma(3/8 + 1/2) - digamma(3/8)) *
	(digamma(1/8 + 1/2) - digamma(1/8)) / 32;


### I( x * atan(x) / (x^4 + 1) )
integrate(\(x) x * atan(x) / (x^4 + 1), 0, 1)
(digamma(7/8) - digamma(3/8)) * (digamma(5/8) - digamma(1/8)) / 64;

# Note:
# - uses double integral: on [0,1] x [0,1]
#   I(I( u^2 / ((u^2*t^2 + 1) * (u^4 + 1)) ));


### I( x^3 * atan(x) / (x^4 + 1) )
integrate(\(x) x^3 * atan(x) / (x^4 + 1), 0, 1)
(digamma(3/8 + 1/2) - digamma(3/8))^2 / 128 +
	- (digamma(1/8 + 1/2) - digamma(1/8))^2 / 128 +
	+ Catalan/2;

### on [0, Inf]
integrate(\(x) x^3 * atan(x) / (x^4 + 1) - pi/2/(x+1), 0, Inf)
- pi * log(2*cos(pi/8)) / 2

### on [0, Inf]
integrate(\(x) x * atan(x) / (x^4 + 1), 0, Inf)
integrate(\(x) x * (pi/2 - atan(x)) / (x^4 + 1), 0, Inf)
pi^2/16

### I( atan(x) / (x^4 + 1) )
integrate(\(x) atan(x) / (x^4 + 1), 0, Inf)
(pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 64 +
	+ pi/32 * (digamma(5/8) - digamma(1/8));

### I( x^2 * atan(x) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x) / (x^4 + 1), 0, Inf)
integrate(\(x) (pi/2 - atan(x)) / (x^4 + 1), 0, Inf)
(pracma::psi(1, 3/8) - pracma::psi(1, 7/8)) / 64 +
	+ pi/32 * (digamma(1/8) - digamma(5/8)) + pi^2/8 / sin(pi/4);


### Pow: 6

###
# x => 1/x;
integrate(\(x) x^2 * atan(x) / (x^6 + 1), 0, Inf)
integrate(\(x) x^2 * (pi/2 - atan(x)) / (x^6 + 1), 0, Inf)
pi^2 / 24;

###
integrate(\(x) (x^4 + 1) * atan(x) / (x^6 + 1), 0, Inf)
pi^2 / 6

###
integrate(\(x) (x^3 + x) * atan(x) / (x^6 + 1), 0, Inf)
pi^2 / 12 / sin(2*pi/3)
pi^2 / (6*sqrt(3))


# Helper:
integrate(\(x) atan(x) / x, 0, 1)
Catalan

# I( atan(x^n) / x )
n = sqrt(3);
integrate(\(x) atan(x^n) / x, 0, 1)
Catalan / n;

### Fraction Decomposition:
# u^4 / ((u^2*t^2 + 1) * (u^4 + 1))
# 1 / (u^2*t^2 + 1) - Fr0;
# =>
# Fr0 + (u^2*t^2 - 1) / ((u^4 + 1) * (t^4 + 1));
# where Fr0 = 1 / ((u^2*t^2 + 1) * (u^4 + 1));
# and (u^2*t^2 - 1) / Prod() is decomposable into 2 separate integrals:
# I( u^2 / (u^4 + 1) )^2 -  I( 1 / (u^4 + 1) )^2;


### Note: Lim p -> 0
p = 1E-4
integrate(\(x) atan(x^(1/p)) / p, 0, 1)
pi/(4*p) - integrate(\(x) 1/p^2 * x^(1/p) / (x^(2/p) + 1), 0, 1)$value
pi/(4*p) - integrate(\(x) 1/p * x^p / (x^2 + 1), 0, 1)$value
- (pracma::psi(1, 3/4) - pracma::psi(1, 1/4)) / 16
Catalan


##############

### Simple Fractions

### I( atan(x) / (x*(x^2+1)) )
integrate(\(x) atan(x) / (x*(x^2+1)), 0, Inf)
pi*log(2)/2

### on [0, 1]
integrate(\(x) atan(x) / (x*(x^2+1)), 0, 1)
pi*log(2)/8 + Catalan/2


### I( atan(x) / ((x+1) * (x^2+1)) )
integrate(\(x) atan(x) / ((x+1) * (x^2+1)), 0, Inf)
1/8 * pi*log(2) + pi^2/16 - Catalan/2

### on [0, 1]
integrate(\(x) atan(x) / ((x+1) * (x^2+1)), 0, 1)
pi*log(2)/8 + pi^2/64 - 1/4 * Catalan

# Note: see next section;


### I( atan(x) / (x*(x+1)*(x^2+1)) )
# Maths 505: A MONSTER INTEGRAL!!!
# int 0 to infty arctan(x)/(x(x+1)(x^2+1))
# https://www.youtube.com/watch?v=u-FKjn_83l8

### on [0, Inf]
integrate(\(x) atan(x) / (x*(x+1)*(x^2+1)), 0, Inf)
integrate(\(x) x / (tan(x) * (tan(x)+1)), 0, pi/2)
integrate(\(x) x / tan(x) - x / (tan(x)+1), 0, pi/2)
3/8 * pi*log(2) - pi^2/16 + Catalan/2

### on [0, 1]
integrate(\(x) atan(x) / (x*(x+1)*(x^2+1)), 0, 1)
- pi^2/64 + 3/4 * Catalan


# Helper: I( x * cos(x) / (cos(x) + sin(x)) )
x = pi/5
integrate(\(x) x * cos(x) / (cos(x) + sin(x)), 0, x)
x*log(cos(x) + sin(x)) / 2 + x^2/4 +
	- integrate(\(x) 1/2 * log(cos(x) + sin(x)), 0, x)$value
x*log(cos(x) + sin(x)) / 2 + x^2/4 - x*log(2)/4 +
	- integrate(\(x) 1/2 * log(sin(x)), pi/4, x+pi/4)$value

# on [0, pi/4]
integrate(\(x) x * cos(x) / (cos(x) + sin(x)), 0, pi/4)
pi^2/64 + pi*log(2)/8 - Catalan/4

#
x = pi/5
integrate(\(x) cos(x) / (cos(x) + sin(x)), 0, x)
(log(cos(x) + sin(x)) + x) / 2;

