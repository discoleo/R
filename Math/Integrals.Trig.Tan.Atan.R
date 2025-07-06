##################
##
## Leonard Mada
##
## Integrals: Trig
## Variants: ATAN


### Atan-Fractions


### Examples:

### on [0, Inf]
# I( atan(k * x^2) / (x^4 + 1) )
# I( x^2 * atan(k * x^2) / (x^4 + 1) )
# I( atan(k * x^4) / (x^4 + 1) )
# I( x^2 * atan(k * x^4) / (x^4 + 1) )

# Note:
# - various examples are still in file:
#   Integrals.Trig.Tan.R;


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

#####################
#####################

### Basic Integrals

### I( atan(x) / x )
integrate(\(x) atan(x) / x - pi/2 / (x+1), 0, Inf)
0;

### on [0, 1]
integrate(\(x) atan(x) / x, 0, 1)
Catalan

### on [0, tan(pi/3)]
integrate(\(x) atan(x) / x, 0, tan(pi/3))
integrate(\(x) log(cos(x)) - log(sin(x)), 0, pi/3)$value +
	+ pi/3 * log(tan(pi/3));
pi/3 * log(tan(pi/3)) + sqrt(3)/(8*36) *
	( 7 * pracma::psi(1, 1/6) - 7 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) +
	- pracma::psi(1, 1/12) + pracma::psi(1, 11/12) +
	+ pracma::psi(1, 5/12) - pracma::psi(1, 7/12) );

### on [0, tan(pi/6)]
integrate(\(x) atan(x) / x, 0, tan(pi/6))
integrate(\(x) log(cos(x)) - log(sin(x)), 0, pi/6)$value +
	+ pi/6 * log(tan(pi/6));
pi/6 * log(tan(pi/6)) + sqrt(3)/(4*12^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) +
	+ 3 * pracma::psi(1, 1/6) - 3 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) );

### on [0, tan(pi/5)]
integrate(\(x) atan(x) / x, 0, tan(pi/5))
integrate(\(x) log(cos(x)) - log(sin(x)), 0, pi/5)$value +
	+ pi/5 * log(tan(pi/5));
sn  = sin(2*pi/5 * c(1,2));
sn3 = sin(6*pi/10 * c(1,2,3,4));
pi/5 * log(tan(pi/5)) +
	+ (sn[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10))) / (2*100) +
	+ (sn3[1] * (pracma::psi(1, 1/20) - pracma::psi(1, 19/20) +
		- pracma::psi(1, 9/20) + pracma::psi(1, 11/20)) +
	- sn3[2] * (pracma::psi(1, 2/20) - pracma::psi(1, 18/20) +
		- pracma::psi(1, 8/20) + pracma::psi(1, 12/20)) +
	+ sn3[3] * (pracma::psi(1, 3/20) - pracma::psi(1, 17/20) +
		- pracma::psi(1, 7/20) + pracma::psi(1, 13/20)) +
	- sn3[4] * (pracma::psi(1, 4/20) - pracma::psi(1, 16/20) +
		- pracma::psi(1, 6/20) + pracma::psi(1, 14/20)) ) / (2*20^2);
# TODO: simplify;

### on [0, tan(2/5 * pi)]
integrate(\(x) atan(x) / x, 0, tan(pi*2/5))
integrate(\(x) log(cos(x)) - log(sin(x)), 0, pi*2/5)$value +
	+ pi*2/5 * log(tan(pi*2/5));
sn2 = sin(4*pi/5 * c(1,2));
sn1 = sin(2*pi/10 * c(1,2,3,4));
pi*2/5 * log(tan(pi*2/5)) +
	+ (sn2[1] * (pracma::psi(1, 1/10) - pracma::psi(1, 9/10) +
		+ pracma::psi(1, 4/10) - pracma::psi(1, 6/10)) +
	- sn2[2] * (pracma::psi(1, 2/10) - pracma::psi(1, 8/10) +
		+ pracma::psi(1, 3/10) - pracma::psi(1, 7/10)) ) / (2*100) +
	+ (sn1[1] * (pracma::psi(1, 1/20) - pracma::psi(1, 19/20) +
		- pracma::psi(1, 9/20) + pracma::psi(1, 11/20)) +
	- sn1[2] * (pracma::psi(1, 2/20) - pracma::psi(1, 18/20) +
		- pracma::psi(1, 8/20) + pracma::psi(1, 12/20)) +
	+ sn1[3] * (pracma::psi(1, 3/20) - pracma::psi(1, 17/20) +
		- pracma::psi(1, 7/20) + pracma::psi(1, 13/20)) +
	- sn1[4] * (pracma::psi(1, 4/20) - pracma::psi(1, 16/20) +
		- pracma::psi(1, 6/20) + pracma::psi(1, 14/20)) ) / (2*20^2);


### Other

### I( atan(x) / (1 - x) )
integrate(\(x) atan(x) / (1 - x) - pi/4 / (1-x), 0, 1)
pi*log(2)/8 - Catalan


# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


### I( atan(x) / x^2 )
integrate(\(x) - (atan(x) / x^2 - 1/x), 0, 1)
pi/4 + log(2)/2 - 1


### I( atan(1-x) / x )
integrate(\(x) atan(1-x) / x - pi/4 / x, 0, 1)
pi*log(2)/8 - Catalan

### I( atan(x) / (x+1) )
integrate(\(x) atan(x) / (x+1), 0, 1)
pi * log(2)/8

### I( atan(x+1) / x )
integrate(\(x) atan(x+1) / x - pi/4/x, 0, 1)
integrate(\(x) -1/2 * log(x) / (x^2 + 1), 1/2, 1)$value +
	- pi*log(2)/8 + (Catalan + atan(1/2)*log(2)) / 2;
integrate(\(x) 1/2 * log(x) / (x^2 + 1), 0, 1/2)$value +
	- pi*log(2)/8 + Catalan + atan(1/2)*log(2) / 2;
# TODO

# Note:
integrate(\(x) atan(x+1) / x - pi/4/x, 0, 1)
(- pi*log(2)/8 + Catalan) - 0.48722235829452235711 / 2;
# where:
# sum( (-1)^n / 2^(2*n+1) / (2*n+1)^2 ) # ==
# 0.48722235829452235711;


### I( atan(x/(1-x)) / x )
integrate(\(x) atan(x/(1-x)) / x, 0, 1)
pi*log(2)/4 + Catalan;

### I( atan(x*(1-x)) / x )
integrate(\(x) atan(x*(1-x)) / x, 0, 1)
# TODO


######################

### I( atan(x)^2 / x )

### on [0, 1]
# Note: log(tan(pi/4)) == 0;
integrate(\(x) atan(x)^2 / x, 0, 1)
- 7/8 * pracma::zeta(3) + pi * Catalan / 2;


### on [0, tan(pi/3)]
integrate(\(x) atan(x)^2 / x, 0, tan(pi/3))
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
(pi/3)^2 * log(tan(pi/3)) - pracma::zeta(3) * (7/8 - 1 / (2*3^3)) +
	- sum(cos(2*pi/3) * (pracma::psi(2, 1/3) + pracma::psi(2, 1 - 1/3))) / (4*3^3) +
	+ pi * sin(2*pi/3) * (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) / (3^3) +
	+ pi * sin(2*pi/3)/(2*36) *
	( pracma::psi(1, 1/6) - pracma::psi(1, 5/6) +
	+ pracma::psi(1, 1/3) - pracma::psi(1, 2/3) ) +
	+ sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (4*6^3) +
	- sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))) * pi / (6^3);


### on [0, tan(pi/6)]
integrate(\(x) atan(x)^2 / x, 0, tan(pi/6))
ic = 1:3; id = 1:2; cs = cos(2*pi*ic/6); sn = sin(2*pi*id/6);
(pi/6)^2 * log(tan(pi/6)) +
- pracma::zeta(3) * (7/8 + 1 / (8*3^3)) +
	- sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (2*6^3) +
	+ sum(sn * (pracma::psi(1, id/6) - pracma::psi(1, 1 - id/6))) * 2*pi / (6^3) +
	+ sum(cos(2*pi/3) * (pracma::psi(2, 1/3) + pracma::psi(2, 2/3)) / 4 +
		- sin(2*pi/3) * (pracma::psi(1, 1/3) - pracma::psi(1, 2/3)) * pi
	) / (4*3^3);


### Other Fractions:

### I( atan(x)^2 / (1-x) )
integrate(\(x) atan(1-x)^2 / x - pi^2/16/x, 0, 1)
- 7/32 * pracma::zeta(3) + pi^2 * log(2) / 32 - Catalan * pi / 4;


#################
#################

### ASIN

### I( asin(x)^n / x^m )
# Maths 505: This integral is actually one of your favorite constants
# https://www.youtube.com/watch?v=wCpN0aE3aU0


### I( asin(x) * acos(x) / x )
integrate(\(x) asin(x) * acos(x) / x, 0, 1)
integrate(\(x) asin(x) * (pi/2 - asin(x)) / x, 0, 1)
7/8 * pracma::zeta(3)

### I( asin(x) * acos(x) / x^2 )
integrate(\(x) asin(x) * acos(x) / x^2 - pi/2/x, 0, 1)
integrate(\(x) asin(x) * (pi/2 - asin(x)) / x^2 - pi/2/x, 0, 1)
pi/2 * log(2) + pi/2 - 4*Catalan


### I( asin(x)^2 / x )
integrate(\(x) asin(x)^2 / x, 0, 1)
1/4 * pi^2 * log(2) - 7/8 * pracma::zeta(3)

### I( asin(x) / x^2 )
integrate(\(x) asin(x) / x^2 - 1/x, 0, 1)
- pi/2 + log(2) + 1

### I( asin(x)^2 / x^2 )
integrate(\(x) asin(x)^2 / x^2, 0, 1, rel.tol=1E-9)
integrate(\(x) 2*x / sin(x), 0, pi/2)$value - (pi/2)^2;
- (pi/2)^2 + 4*Catalan


#################
#################

##############
### Powers ###

### I( atan(x)^2 )
integrate(\(x) atan(x)^2, 0, 1)
(pi/4)^2 + pi*log(2)/4 - Catalan

### I( atan(x)^3 )
integrate(\(x) atan(x)^3, 0, 1)
(pi/4)^3 + 63/64 * zeta(3) + 3/32 * pi^2 * log(2) - 3/4 * pi * Catalan

# see x^2 * tan(x):
integrate(\(x) x * atan(x)^2 / (x^2 + 1), 0, 1)
- 21/64 * zeta(3) - 1/32 * pi^2 * log(2) + pi/4 * Catalan


####################
####################

### Fractions

# Note:
# - [refactor] Fractions of Power = 4 moved to file:
#   Integrals.Trig.Tan.Atan.Fractions.P4.R;
# - see also file: Integrals.Log.Fractions.Complex.R;


### I( atan(x) / (x + 1) ) on [0, Inf]
integrate(function(x) atan(x) / (x + 1) - pi/2 / (x+1), 0, Inf)
- pi*log(2)/4 - Catalan;


### Poly: x^2 + 1

### I( atan(x^2) / (x^2 + 1) )
integrate(function(x) atan(x^2) / (x^2 + 1), 0, Inf)
pi^2 / 8

### I( atan(x^n) / (x^2 + 1) )
n = sqrt(3) # independent of n;
integrate(function(x) atan(x^n) / (x^2 + 1), 0, Inf)
pi^2 / 8


### I( x^p * atan(x) / (x^2 + 1) )
p = - 1/sqrt(5);
integrate(\(x) x^p * atan(x) / (x^2 + 1), 0, Inf)
pi/4 / sin(pi*p/2) * (digamma(-p/2) - 2*digamma(-p) - Euler);


##################
##################

### Power = 6

### I( atan(x^3) / (x^6 + 1) )
integrate(\(x) atan(x^3) / (x^6 + 1), 0, Inf)
- pi/12 / sin(pi/3) * (digamma(1/3) - 2*digamma(2/3) - Euler);

### I( x^3 * atan(x^3) / (x^6 + 1) )
integrate(\(x) x^3 * atan(x^3) / (x^6 + 1), 0, Inf)
pi/12 / sin(pi/6) * (digamma(-1/6) - 2*digamma(-1/3) - Euler);


####################
####################

### Log-Combinations

### I( atan(x) * log((1-x)/(1+x)) )
# Maths 505: Another INSANE integral!
# https://www.youtube.com/watch?v=KEDEzVqlAYU

integrate(\(x) atan(x) * log((1-x)/(1+x)), 0, 1)
pi^2/16 - pi/4*log(2) - Catalan;


### I( atan(x) * log(1-x) )
integrate(\(x) atan(x) * log(1-x), 0, 1)
integrate(\(x) ((1-x)*log(1-x) + x) / (x^2+1) - pi/4, 0, 1)
5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan;

### I( atan(x) * log(1+x) )
integrate(\(x) atan(x) * log(1+x), 0, 1)
-1/96*pi^2 + 3/8*pi*log(2) - log(2)^2 / 8 + log(2)/2 - pi/4;


### I( atan(x) * atan(1-x) )
integrate(\(x) atan(x) * atan(1-x), 0, 1)
integrate(\(x) (1-2*x) * atan(1-x) / (x^2+1), 0, 1)
# TODO

### I( x * atan(1-x) / (x^2+1) )
integrate(\(x) x * atan(1-x) / (x^2+1), 0, 1)
# using Im(Li2((3+1i)/5)) = 0.598496790125 / 2;
atan(1/2) * log(2) / 2 - 0.598496790125 / 2 +
	+ log(5) * (pi/4 + atan(1/3) - atan(1/2)) / 4;
# TODO


### I( atan(x) * log(1-x) / x )
integrate(\(x) atan(x) * log(1-x) / x, 0, 1)
integrate(\(x) 1/2 * log(x) * log(x^2+1) / (x^2+1), 0, 1)$value +
	+ Catalan * log(2) - pi^3 * 3/64;
# TODO


#
integrate(\(x) atan(1-x) / (x^2+1), 0, 1)
1.393582/2 + # (Li2((3-1i)/5) + Li2((3+1i)/5)) / 2 +
	+ (log(5)*log(5/4)/4 - pi^2/8  + pi/4*atan(1/2) + atan(1/2)*atan(1/3))/2;

#
integrate(\(x) x * atan(1-x) / (x^2+1), 0, 1)
-0.59849679 * 1i/2i + # (Li2((3-1i)/5) - Li2((3+1i)/5)) / 2i +
	+ (2*log(2)*atan(1/2) + log(5)*(pi/4 + atan(1/3) - atan(1/2))) / 4;


# Gen 1: TODO
k = 2
integrate(\(x) atan(x) * log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


###
integrate(\(x) atan(x) * log((1-x^2)/(1+x^2)) / (x^2+1), 0, tan(pi/8))
(pi*Catalan - 1/4 * pi^2 * log(2) - 21/16 * pracma::zeta(3)) / 32

# Weierstrass substitution, see file Integrals.Log.Trig.R;

####################
####################

### I( log(x) * atan(x) / x )
# 1.) Maths 505: a superb integral sprinkled with some fourier analysis
# https://www.youtube.com/watch?v=aHCLS4l65VE
# 2.a) see also: Sums.Fractions.Higher.R;
# 2.b) Flammable Maths: An AMAZING Journey of Series Evaluation!
# Calculating Euler's Sum! [ Series pi^3/32 (-1)^k/(2k+1)^3 ]

integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32


### I( log(x) * atan(x) )
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2/48 - pi/4 + log(2)/2

# Derivation:
id = seq(0, 20000)
- sum((-1)^id * (1/(2*id+1) - 1/(2*id+2))) + pi^2/48
- sum((-1)^id / (2*id+1)) + pi^2/48 + log(2)/2
pi^2/48 - pi/4 + log(2)/2


### I( log(x) * atan(x^p) / x )
integrate(\(x) log(x) * atan(x^2) / x, 0, 1)
- pi^3 / 128

### Gen:
p = sqrt(5)
integrate(\(x) log(x) * atan(x^p) / x, 0, 1)
- pi^3 / (32 * p^2)

###############

### I( log(x^2 + 1) * atan(x) )
integrate(\(x) log(x^2+1) * atan(x), 0, 1)
3/48*pi^2 + pi*log(2)/4 - log(2)^2 / 4 + log(2) - pi/2;

### I( log(x^2 + 1) * atan(x) / (x^2+1) )
integrate(\(x) log(x^2+1) * atan(x) / (x^2+1), 0, 1)
(pi/4)^2 * log(2) + 21/64 * pracma::zeta(3) - pi * Catalan/4;

### I( x * log(x^2 + 1) * atan(x) / (x^2+1) )
integrate(\(x) x * log(x^2+1) * atan(x) / (x^2+1), 0, 1)
integrate(\(x) -1/4 * log(x^2+1)^2 / (x^2+1), 0, 1)$value +
	+ log(2)^2 * pi/16;
integrate(\(x) -1/24 * log(x^2+1)^3 / x^2, 0, 1)$value +
	+ log(2)^2 * pi/16 - log(2)^3 / 24;
# TODO


### I( log(x^2+1)^2 * atan(x) / x^2 )
integrate(\(x) log(x^2+1)^2 * atan(x) / x^2, 0, 1)
pi^2 * log(2) / 4 - log(2)^2 * pi/4 - log(2)^3/6 +
	+ 23/16 * pracma::zeta(3) - pi * Catalan;


# Helper
integrate(\(x) log(x+1)^2 / x, 0, 1)
pracma::zeta(3) / 4;
#
integrate(\(x) log(x)^2 * (1/(x-1) - 1/x), 1, 2)
pracma::zeta(3) / 4 - log(2)^3/3;


### Other:
integrate(\(x) log(x) * x / (x^2 + 1), 0, 1)
- pi^2/48


### Pow = 2
p = sqrt(5)
integrate(\(x) log(x)^2 * x^(p-1) / (x^(2*p) + 1), 0, 1)
pi^3 / (16 * p^3)

###
integrate(\(x) log(x)^2 / (x^2 + 1), 0, 1)
pi^3 / 16


###
integrate(\(x) log(x^2 + 1) / x, 0, 1)
pi^2/24


################

### I( asin(x)^2 * log(x) )
# Maths 505: I don't have a title for this integral: int (0,1) (arcsin(x))^2 ln(x)
# https://www.youtube.com/watch?v=ZhcYWV9imrU
# Series-Expansion log(sin(x));

integrate(\(x) asin(x)^2 * log(x), 0, 1)
integrate(\(x) x^2 * log(sin(x)) * cos(x), 0, pi/2)
# Alternative: Integration by parts
integrate(\(x) - (x^2*sin(x) + 2*x*cos(x) - 2*sin(x)) / tan(x), 0, pi/2);
- pi^2/4 - 4*Catalan + 6;


### I( asin(x)^2 * log(1 - x^2) )
integrate(\(x) asin(x)^2 * log(1-x^2), 0, 1)
-7*pracma::zeta(3) + pi^2*log(2) / 2 - pi^2/2 - 4*log(2) + 12;

### I( x^2 * log(cos(x)) * cos(x) )
integrate(\(x) x^2 * log(cos(x)) * cos(x), 0, pi/2)
-7/2 * pracma::zeta(3) + pi^2*log(2) / 4 - pi^2/4 - 2*log(2) + 6;

# Derivation:
x = pi/2 - 1E-5;
integrate(\(x) (x^2*sin(x) + 2*x*cos(x) - 2*sin(x)) * tan(x), 0, x)$value +
	+ (x^2*sin(x) + 2*x*cos(x) - 2*sin(x)) * log(cos(x));
integrate(\(x) (x^2*sin(x) - 2*sin(x)) * tan(x), 0, x)$value +
	+ (x^2*sin(x) - 2*sin(x)) * log(cos(x)) + 2;
integrate(\(x) x^2*sin(x) * tan(x), 0, x)$value +
	+ (x^2*sin(x) - 2*sin(x)) * log(cos(x)) +
	+ 2*sin(x) + log(1-sin(x)) - log(1+sin(x)) + 2;
integrate(\(x) x*(x-pi/2)*sin(x) * tan(x), 0, pi/2)$value +
	+ pi^2*log(2) / 4 - pi^2/4 + pi/2 - pi*Catalan - 2*log(2) + 4;
-7/2 * pracma::zeta(3) + pi^2*log(2) / 4 - pi^2/4 - 2*log(2) + 6;


# Lim: x -> pi/2
x = pi/2 - 1E-6;
# x = Const("pi", 200) / 2 - mpfr("1E-24", 200);
log(1-sin(x)) - 2*sin(x) * log(cos(x))
log1p(-sin(x)) - 2*sin(x) * log(cos(x))
- log(2);

# Lim: x -> pi/2
x = pi/2 - 1E-6;
pi^2/4 * (-sin(x) + (log(1+sin(x)) - log(1-sin(x)))/2) +
	+ x^2*sin(x) * log(cos(x)) + 2*sin(x) - log(1+sin(x));
pi^2*log(2) / 4 - pi^2/4 - log(2) + 2;


### Helper:

###
integrate(\(x) x*(x-pi/2)*sin(x) * tan(x), 0, pi/2)
integrate(\(x) (2*x-pi/2)*(sin(x) - (log(1+sin(x)) - log(1-sin(x)))/2), 0, pi/2);
integrate(\(x) x*(log(1-sin(x)) - log(1+sin(x))), 0, pi/2)$value +
	- pi/2 * (- 2*Catalan + 1) + 2;
-7/2 * pracma::zeta(3) + pi*Catalan - pi/2 + 2;

#
integrate(\(x) x*log(1-sin(x)), 0, pi/2)
integrate(\(x) (pi/2-x)*log(2 - 2*cos(x/2)^2), 0, pi/2)
integrate(\(x) (pi/2-x)*log(2*sin(x/2)^2), 0, pi/2)
integrate(\(x) pi/2*log(2) + pi*log(sin(x/2)) - x*log(2) - 2*x*log(sin(x/2)), 0, pi/2)
integrate(\(x) - 8*x*log(sin(x)), 0, pi/4)$value +
	- 3/8* pi^2*log(2) - pi*Catalan;
8*(pi^2 * log(2)/32 + pi * Catalan / 8 - 35/128 * pracma::zeta(3)) +
	- 3/8* pi^2*log(2) - pi*Catalan;
-35/16 * pracma::zeta(3) - 1/8 * pi^2*log(2);

#
integrate(\(x) x*log(1+sin(x)), 0, pi/2)
21*pracma::zeta(3)/16 - 1/8 * pi^2*log(2);


###
x = pi/2 - 1E-6;
integrate(\(x) (x-pi/2)*sin(x) * tan(x), 0, pi/2)
integrate(\(x) sin(x) - (log(1+sin(x)) - log(1-sin(x)))/2, 0, x)$value +
	+ (x - pi/2) * (-sin(x) + (log(1+sin(x)) - log(1-sin(x)))/2);
integrate(\(x) sin(x) - (log(1+sin(x)) - log(1-sin(x)))/2, 0, pi/2);
- 2*Catalan + 1;

###
z = pi/5;
integrate(\(x) sin(x) * tan(x), 0, z)
integrate(\(x) x^2 / (1-x^2), 0, sin(z))
-sin(z) + (log(1+sin(z)) - log(1-sin(z)))/2;


########################

### I( log(atan(x)) )
integrate(\(x) log(atan(x)), 0, 1)
integrate(\(x) - tan(x) / x, 0, pi/4)$value + log(pi/4);
# TODO


###############
###############

### Experiments

### ATAN * ATAN

# Related Materials:
# 1. Dr. Peyam: What is convolution?
#    https://www.youtube.com/watch?v=MkdPzDxUkz0


FUN = function(x, normalize = FALSE) {
	y = sapply(x, \(lim) {
			integrate(\(x) atan(x) * atan(lim - x), 0, lim)$value;
	});
	if(normalize) y = y / x;
	return(y);
}

#
up = 3; # up = 10; # up = 50;
curve(FUN(x), 0, up)
curve(exp(x) - 1, add=TRUE, col="red")
curve(tan(pi/2 * x/up), add=TRUE, col="green")
curve(eval(x), add=TRUE, col="purple")
curve(eval(pi/2*x), add=TRUE, col="pink")

#
up = 10; # up = 30; # up = 500;
curve(FUN(x, normalize = T), 0, up)
curve(atan(x), add=TRUE, col="pink")

# Limit: lim -> Inf
# library(Rmpfr)
lim = mpfr("1E+8", 240)
integrate(\(x) {
	x = mpfr(x, 240);
	y = atan(x) * atan(lim - x) / lim;
	as.numeric(y);
}, 0, as.numeric(lim))
pi^2 / 4

# Note: the proof is left as an exercise for the reader;

