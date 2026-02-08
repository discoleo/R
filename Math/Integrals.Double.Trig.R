#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Trig Variants
##
## v.0.1a

### Double Integrals

### Examples:
# I( log(cos(x-y)) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

#####################
#####################

############
### Trig ###


### I( sin(pi/2*x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x*y), 0, 1)$value), 0, 1)
(sum(pracma::expint(pi/2*1i*c(-1,1)))/2 + log(pi/2) + Euler) * 2/pi;


### Simple Fractions

### I( sin(x) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(x) / (x + y), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
# TODO


### I( sin(pi/2*(x+y)) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (x + y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) 4/pi*sin(x)*cos(y) / (x + y), 0, pi/2)$value), 0, pi/2)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	2*sin(x) / (x + y), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12) $value +
	- (sum(pracma::expint(pi/2*1i*c(-1,1))) + 2*log(pi/2) + 2*Euler);
# TODO


### I( sin(pi/2*x) * sin(pi/2*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(pi/2*x) * sin(pi/2*y) / (x + y), 0, 1)$value), 0, 1)


### I( sin(pi/2*x) * cos(pi/2*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(pi/2*x) * cos(pi/2*y) / (x + y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(x) / (x + y), 0, pi/2)$value), 0, pi/2) $value +
	- (sum(pracma::expint(pi/2*1i*c(-1,1)))/2 + log(pi/2) + Euler);
# TODO


### I( sin(pi*(x+y)) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*(x+y)) / (x + y), 0, 1)$value), 0, 1)


### Fractions: x*y-Type


### I( sin(pi/2*x) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x) / (1 + x*y), 0, 1)$value), 0, 1)

### I( cos(pi/2*x) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(pi/2*x) / (1 + x*y), 0, 1)$value), 0, 1)


### I( sin(pi/2*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)

### I( cos(pi/2*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(pi/2*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)


### I( sin(pi*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)

### I( cos(pi*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(pi*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)


### I( sin(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(x+y) / (1 + x*y), 0, pi/2)$value), 0, pi/2)

### I( cos(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(x+y) / (1 + x*y), 0, pi/2)$value), 0, pi/2)


### I( sin(pi/4*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/4*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)

### I( cos(pi/4*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(pi/4*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)


### Div: (1 - x*y)

### I( sin(pi/2*(x+y)) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (1 - x*y), 0, 1)$value), 0, 1)


### I( sin(pi/2*x) * sin(pi/2*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(pi/2*x) * sin(pi/2*y) / (1 - x*y), 0, 1)$value), 0, 1)
	
### I( cos(pi/2*x) * cos(pi/2*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	cos(pi/2*x) * cos(pi/2*y) / (1 - x*y), 0, 1)$value), 0, 1)
	

### Trig( Prod )

### I( sin(pi/2*x*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x*y) / (x + y), 0, 1)$value), 0, 1)

### I( sin(pi*x*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*x*y) / (x + y), 0, 1)$value), 0, 1)

### I( cos(pi/2*x*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(pi/2*x*y) / (x + y), 0, 1)$value), 0, 1)


### Fraction: 1 +/- x*y

### I( sin(pi/2*x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x*y) / (1 - x*y), 0, 1)$value), 0, 1)

### I( sin(pi/2*x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x*y) / (1 + x*y), 0, 1)$value), 0, 1)


######################
### Trig-Fractions ###

### I( (sin(x) + sin(y)) / (cos(x) + cos(y)) )
# Maths 505: A beautifully symmetric double integral
# resulting in an important constant
# https://www.youtube.com/watch?v=MVgDeJaLcd0

###
integrate(\(x) sapply(x, \(y) integrate(
	\(x) (sin(x) + sin(y)) / (cos(x) + cos(y)), 0, pi/2)$value), 0, pi/2)
integrate(\(x) sapply(x, \(y) integrate(
	\(x) 2*sin(x) / (cos(x) + cos(y)), 0, pi/2)$value), 0, pi/2)
4 * Catalan


### Diff-type
up = pi/2 - 1E-12; # Numeric instability;
integrate(\(x) sapply(x, \(y) integrate(\(x) {
	x = mpfr(x, 200); y = mpfr(y, 200);
	as.numeric((sin(x) - sin(y)) / (cos(x) - cos(y))); }, 0, up)$value), 0, pi/2)
-4 * Catalan


### I( sin(x+y) / (cos(x) + cos(y)) )
integrate(\(x) sapply(x, \(y) integrate(
	\(x) sin(x+y) / (cos(x) + cos(y)), 0, pi/2)$value), 0, pi/2)
pi - 2*log(2)


### I( sin(x)*cos(y) / (cos(x) + cos(y)) )
integrate(\(x) sapply(x, \(y) integrate(
	\(x) sin(pi/2*x)*cos(pi/2*y) / (cos(pi/2*x) + cos(pi/2*y)), 0, 1)$value), 0, 1)
(pi/2 - log(2)) * 4/pi^2

# Helper:
integrate(\(x) cos(x)^2 * log(cos(x)), 0, pi/4)
- (2*pi*log(2) - pi - 4*Catalan + 2*log(2) + 2) / 16;



###
integrate(\(x) sapply(x, \(y) integrate(
	\(x) sin(pi/2*x*y) / (cos(pi/2*x) + cos(pi/2*y)), 0, 1)$value), 0, 1)
# TODO


################
### Radicals ###

### I( sin(x) / sqrt(1 - sin(x)^2*sin(y)^2*sin(z)^2) )
# 1. Hmath: Расправимся с тройным интегралом с помощью рядов
#    [in Russian]
# 2. Hmath: What if we find the average length of an ellipse?
#    https://www.youtube.com/watch?v=542frBtmBgA

# TODO


### I( sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
gamma(1/4)^4 / gamma(1/2)^2 / 32 + gamma(3/4)^4 / gamma(3/2)^2 / 8;
beta(1/4, 1/4)^2 / 32 + beta(3/4, 3/4)^2 / 8;

### I( 1 / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2)$value), 0, pi/2)
gamma(1/4)^4 / gamma(1/2)^2 / 16;


### Trig * SQRT

### I( sin(x) * sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(1 - sin(x)^2*sin(y)^2) * sin(x), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
pi^2 / 8;

### I( cos(x) * sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sqrt(1 - sin(x)^2*sin(y)^2) * cos(x), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
Catalan + 1/2;


### I( sin(x) / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(x) / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2)$value), 0, pi/2)
pi^2 / 4;


### I( cos(x) / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	cos(x) / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2)$value), 0, pi/2)
2 * Catalan;


### I( sin(2*x) / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(2*x) / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2)$value), 0, pi/2)
2;

### I( sin(4*x) / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(4*x) / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2)$value), 0, pi/2)
- 4/9;

### I( sin(5*x) / sqrt(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(5*x) / sqrt(1 - sin(x)^2*sin(y)^2), 0, pi/2, rel.tol=1E-10)$value), 0, pi/2, rel.tol=1E-10)
pi^2 / 16;


##################
### Log & Trig ###

### I( log(x)*log(y) * cos(x+y) / sqrt(x*y) ) on [0, Inf] x [0, Inf]
# Maths 505: The best double integral you'll see this week
# https://www.youtube.com/watch?v=9TjOahihJuE
# Note: separation of x vs y;

# upper = Inf
p = 1/2; up = 2000;
# p = 2/5; up = 3200; # p = 3/5; up = 2000*pi;
integrate(\(x)   x^(p-1) * log(x) * cos(x), 0, up, subdivisions = 1024)$value^2 +
integrate(\(x) - x^(p-1) * log(x) * sin(x), 0, up, subdivisions = 1024)$value^2
gamma(p)^2 * (digamma(p)*cos(pi*p/2) - pi/2*sin(pi*p/2))^2 +
- gamma(p)^2 * (digamma(p)*sin(pi*p/2) + pi/2*cos(pi*p/2))^2;
#
gamma(p)^2 * (cos(pi*p)*(digamma(p)^2 - pi^2/4) - pi*digamma(p)*sin(pi*p));


library(Rmpfr)
# fails:
fx = \(x, y) {
	x = mpfr(x, 128); y = mpfr(y, 128);
	as.numeric(log(x)*log(y) * cos(x+y) / sqrt(x*y)); }
sc = 400; integrate(\(x) sapply(x, \(y)
	integrate(fx, 0, sc*pi, y=x, subdivisions=80024)$value), 0, sc*pi, subdivisions=80024)



###################
### Log( Trig ) ###

### I( log(sin(x)/cos(y) + cos(x)/sin(y)) )
# Maths 505: A beautiful iterated integral
# https://www.youtube.com/watch?v=Cyvgh9RDaV8

###
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(sin(x)/cos(y) + cos(x)/sin(y)), 0, pi/2)$value), 0, pi/2)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(2 * cos(x-y) / sin(y)), 0, pi/2)$value), 0, pi/2)
7/8 * pracma::zeta(3) + (pi/2)^2 * log(2);


### I( log(cos(x-y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x-y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
7/8 * pracma::zeta(3) - (pi/2)^2 * log(2);

### I( x * log(cos(x-y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(cos(x-y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
pi/4 * (7/8 * pracma::zeta(3) - (pi/2)^2 * log(2));

### I( x * log|sin(x-y)| ) on [0, y]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(sin(y-x)), 0, y, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
- pi/8 * (pracma::zeta(3) + pi^2*log(2) / 6);

### I( log(cos(x-y)) / (x-y) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x-y)) / (x-y), 0, y, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
- pi/2 * (log(2) + integrate(\(x) log(cos(x)) / x, 0, pi/2, rel.tol=1E-13)$value);
# TODO


### I( log(cos(x) + cos(y)) )
# see next Section: cos(x-y)*cos(x+y);
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x) + cos(y)), 0, pi/2)$value), 0, pi/2)
7/4 * pracma::zeta(3) - (pi/2)^2 * log(2);

### I( log(cos(x) - cos(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x) - cos(y)), 0, y, rel.tol=1E-10)$value), 0, pi/2, rel.tol=1E-10)
- (7/8 * pracma::zeta(3) + 1/8 * pi^2 * log(2));

### Gen: on [0, k/n * pi]
up = 3*pi/17; # up = sqrt(2) / pi;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x) + cos(y)), 0, up, rel.tol=1E-10)$value), 0, up, rel.tol=1E-10)
integrate(\(x) 4*(up - x) * log(cos(x)), 0, up)$value + log(2) * up^2;
# TODO: explicit formula;

###
up = 3*pi/17; # up = sqrt(2) / pi;
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x) - cos(y)), 0, y, rel.tol=1E-10)$value), 0, up, rel.tol=1E-10)
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(2*sin((x+y)/2) * sin((y-x)/2)), 0, y, rel.tol=1E-10)$value), 0, up, rel.tol=1E-10)
# TODO


### I( log(1 - sin(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 - sin(x)*sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
- (7/2 * pracma::zeta(3) + pi^2 * log(2) / 2 - 2*pi * Catalan);

### I( log(1 + sin(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(1 + sin(x)*sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
7/2 * pracma::zeta(3) - pi^2 * log(2) / 2;

### I( x * log(1 - sin(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(1 - sin(x)*sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 256 + pi^2 * Catalan +
	- pi/2 * (21/16 * pracma::zeta(3) + (pi/2)^2 * log(2));

### I( x * log(1 + sin(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(1 + sin(x)*sin(y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 256 +
	- pi/2 * (21/16 * pracma::zeta(3) + (pi/2)^2 * log(2));

# Derivation:

### I( x * log(1 - sin(x)^2*sin(y)^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * log(1 - sin(x)^2*sin(y)^2), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
- pi * (21/16 * pracma::zeta(3) + (pi/2)^2 * log(2) - pi * Catalan);

### Diff:
integrate(\(x) sapply(x, \(y) integrate(\(x)
	x * (log(1 - sin(x)*sin(y)) - log(1 + sin(x)*sin(y))), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 + pi^2 * Catalan;


### Trig * Log( Trig )

### I( sin(x) * log(cos(x) + cos(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(x) * log(cos(x) + cos(y)), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
2*Catalan - log(2) - pi*log(2)/2;

### I( cos(x) * log(cos(x) + cos(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	cos(x) * log(cos(x) + cos(y)), 0, pi/2, rel.tol=1E-13)$value), 0, pi/2, rel.tol=1E-13)
- pi*log(2)/2 + pi/2;


### I( log(tan(x) + tan(y)) )
# Maths 505: A brutal iterated integral!
# https://www.youtube.com/watch?v=UwEgcCZ_SqU

# on [0, pi/4]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(x) + tan(y)), 0, pi/4)$value), 0, pi/4)
7/64 * pracma::zeta(3) + pi^2 * log(2)/16 - pi/4 * Catalan


### I( log(sin(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, pi/4)$value), 0, pi/4)
7/64 * pracma::zeta(3) - pi^2 * log(2) / 16

### I( log(cos(x-y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/4)$value), 0, pi/4)
21/64 * pracma::zeta(3) - pi^2 * log(2) / 16;

### I( log(cos(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x+y)), 0, pi/4)$value), 0, pi/4)
7/64 * pracma::zeta(3) - pi^2 * log(2) / 16;


# Derivation: I( log(cos(x-y)) )
integrate(\(x) pi/4 * log(cos(x - pi/4)) +
	+ sapply(x, \(y) integrate(\(x) x * tan(x-y), 0, pi/4)$value), 0, pi/4);
integrate(\(x) sapply(x, \(y) integrate(\(x) x * tan(x-y), 0, pi/4)$value), 0, pi/4)$value +
	- pi/4 * (pi * log(2) / 4 - Catalan / 2);
integrate(\(x) x * log(cos(x-pi/4)) - x*log(cos(x)), 0, pi/4)$value +
	- pi/4 * (pi * log(2) / 4 - Catalan / 2);
integrate(\(x) (pi/4 - x) * log(cos(x)), 0, pi/4)$value +
	- pi^2/32 * log(2) + 21/128 * pracma::zeta(3);
- pi^2/16 * log(2) + 21/64 * pracma::zeta(3);

# I( log(cos(x+y)) )
up = pi/4; # up = pi/6; # up = 3*pi/13;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x+y)), 0, up)$value), 0, up)
integrate(\(x) up * log(cos(x + up)) +
	+ sapply(x, \(y) integrate(\(x) x * tan(x+y), 0, up)$value), 0, up);
integrate(\(x) up * (2*log(cos(2*x)) - log(cos(x))) +
	+ sapply(x, \(y) integrate(\(x) x * tan(x+y), 0, up)$value), 0, up)
integrate(\(x) 2*up * (2*log(cos(2*x)) - log(cos(x))) +
	+ (x*log(cos(x)) - (x+up)*log(cos(x+up))), 0, up);
integrate(\(x) up * (2*log(cos(2*x)) - log(cos(x))) +
	+ x*(log(cos(x)) - log(cos(x+up))), 0, up);
integrate(\(x) 2*up * (2*log(cos(2*x)) - log(cos(x))) +
	- 2*x * (2*log(cos(2*x)) - log(cos(x))), 0, up);

# Simplification for up = pi/4;
integrate(\(x) pi/4 * log(sin(x)) +
	+ x * (log(cos(x)) - log(cos(x+pi/4))), 0, pi/4);
integrate(\(x) x * log(sin(x)) + x * log(cos(x)), 0, pi/4);


### on [0, pi/3]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/3)$value), 0, pi/3)
integrate(\(x) 2*(pi/3 - x) * log(cos(x)), 0, pi/3);
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
-(pi/3)^2 * log(2) + 3/8 * pracma::zeta(3) +
	- sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (4*6^3);
-(pi/3)^2 * log(2) + (1/6 + 3/8) * pracma::zeta(3);

### on [0, pi/6]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/6)$value), 0, pi/6)
(11/2 * pracma::zeta(3) - pi^2 * log(2)) / 36;

### on [0, pi/6] * [pi/6, pi/3]
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/6)$value), pi/6, pi/3)
(17/4 * pracma::zeta(3) - pi^2 * log(2)) / 36;

### on [0, pi/6] * [0, pi/3]
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/6)$value), 0, pi/3)
(39/4 * pracma::zeta(3) - 2*pi^2 * log(2)) / 36;

### on [0, pi/6] * [0, pi/2]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x-y)), 0, pi/6, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
(35/4 * pracma::zeta(3) - 3*pi^2 * log(2)) / 36;


### I( log(cos(x+y)) ) on [0, pi/6]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x+y)), 0, pi/6)$value), 0, pi/6)
integrate(\(x) pi/3 * (2*log(cos(2*x)) - log(cos(x))) +
	- 2*x * (2*log(cos(2*x)) - log(cos(x))), 0, pi/6);
# TODO: full computation;

# Note:
# - for sub-integrals, see: Integrals.Log.Trig.R;
# - alternatively: Clausen function;


### on ANY Fraction of pi
n = 6; k = 1; # n = 7; k = 3; # n = 8; k = 3;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi*k/n)$value), 0, pi*k/n)
integrate(\(x) 2*(pi*k/n - x) * log(cos(x)), 0, pi*k/n);
#
id = seq(2*n); idn = id / (4*n);
sn = (-1)^id * sin(2*k*id*pi/n); cs = (-1)^id * cos(2*k*id*pi/n);
id1 = seq(n); idn1 = id1 / (2*n);
sn1 = -(-1)^id1 * sin(2*k*id1*pi/n);
2*pi*k/n * (sum(sn1 * (pracma::psi(1, idn1) + (-1)^n * pracma::psi(1, 1/2 + idn1))) / (8*n^2) +
	- k/n * pi * log(2)) +
- 2*(-2*(k/n)^2 * (pi/2)^2 * log(2) - 3/16 * pracma::zeta(3) + (
	- sum(sn * (pracma::psi(1, idn) + pracma::psi(1, 1/2 + idn))) * k*pi +
	+ sum(cs * (pracma::psi(2, idn) + pracma::psi(2, 1/2 + idn))) / 16 ) / (32*n^3));


### I( log(cos(x+y)) )
up = pi/6; # up = 3*pi/13;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x+y)), 0, up)$value), 0, up)
integrate(\(x) 2*up * (2*log(cos(2*x)) - log(cos(x))) +
	- 2*x * (2*log(cos(2*x)) - log(cos(x))), 0, up);
# TODO: Explicit formula;


### Other Intervals

### on [0, up] * [0, pi/2]
up = pi/5;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, up)$value), 0, pi/2)
integrate(\(x) up * log(cos(x-up)), 0, pi/2)$value +
integrate(\(x) x * log(cos(x-pi/2)) - x * log(cos(x)), 0, up)$value;
integrate(\(x) up * log(sin(x)), 0, pi/2)$value +
integrate(\(x) (x-up) * log(sin(x)) - (x-up) * log(cos(x)), 0, up)$value;
# TODO: explicit formula;

### on [0, up] * [up, pi/2]
up = pi/5;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, up)$value), up, pi/2)
integrate(\(x) up * log(sin(x)), 0, pi/2)$value +
integrate(\(x) (x-up) * log(sin(x)) + (x-up) * log(cos(x)), 0, up)$value;
# TODO: explicit formula;


### [old]
### on [0, pi/8]^2
n = 8; # n = 12; # EVEN: 4*k;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/n)$value), 0, pi/n)
integrate(\(x) 2*(pi/n - x) * log(cos(x)), 0, pi/n);
p  = 1; id0 = seq(n); n2 = n/2;
id = seq(n/2 - 1); id2 = seq(n/4 - 1);
sn = sin(2*pi*id/n); sn2 = sin(2*pi*id2/n2);
cs = cos(2*pi*id/n); cs2 = cos(2*pi*id2/n2);
idn = id0 / (2*n); snsg = -(-1)^id0 * sin(2*p*id0*pi/n);
pi * sum(snsg * (pracma::psi(1, idn) + (-1)^n * pracma::psi(1, 1/2 + idn))) / (4*n^3) +
	- (2*p-1)*log(2)*(pi/n)^2 + pracma::zeta(3) * (3/8 - 3/n^3) +
	- 1 / (4*n^3) * (
	+ sum(cs2 * (pracma::psi(2, id2/n2) + pracma::psi(2, 1 - id2/n2))) * 2 +
	- sum(sn2 * (pracma::psi(1, id2/n2) - pracma::psi(1, 1 - id2/n2))) * 8*pi +
	- sum( cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) +
	+ sum( sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))) * 4*pi );

# TODO: generalize [DONE] & simplify;


### Log( SIN )

### on [0, up]^2
up = 2*pi/7;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, up)$value), 0, up)
integrate(\(x) up*log(sin(x+up)) - sapply(x, \(y) integrate(\(x) x / tan(x+y), 0, up)$value), 0, up);
integrate(\(x) (up-x)*log(sin(x+up)) + x*log(sin(x)), 0, up);
integrate(\(x) 4*(up-x)*log(sin(2*x)) - 2*(up-x)*log(sin(x)), 0, up);
# TODO: explicit formula;

### on [0, up] * [0, pi/2-up]
up = 2*pi/7;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, up)$value), 0, pi/2-up)
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, up)$value), up, pi/2);
integrate(\(x) up * log(sin(x)), 0, pi/2)$value +
integrate(\(x) (x-up) * log(sin(x)) + (x-up) * log(cos(x)), 0, up)$value;
# TODO: explicit formula;

### on [0, up] * [up, pi/2-up]
up = 2*pi/7;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, up)$value), up, pi/2-up)
integrate(\(x) up * log(sin(x)), 0, pi/2)$value +
integrate(\(x) (x-up) * log(cos(x)) - (x-up) * log(sin(x)) +
	- 4*(up-x) * log(sin(2*x)), 0, up)$value;
integrate(\(x) up * log(sin(x)), 0, pi/2)$value + # - up/2 * pi*log(2);
integrate(\(x) (x-up) * log(cos(x)) - (x-up) * log(sin(x)), 0, up)$value +
integrate(\(x) (x-2*up) * log(sin(x)), 0, 2*up)$value;
# Alternative: sin(2*x) = 2*sin(x)*cos(x);
integrate(\(x) 5*(x-up) * log(cos(x)) + 3*(x-up) * log(sin(x)), 0, up)$value +
	- 1/2 * up*pi*log(2) - 2*up^2*log(2);
# TODO: explicit formula;

# on [up, pi/2-up]^2
# - using Diff(I[up], I[pi/2-up]),
#   where I is on [0, up] * [up, pi/2-up];
up = 3*pi/7;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), up, pi/2-up)$value), up, pi/2-up)
integrate(\(x) 5*(pi/2-up-x) * log(cos(x)) + 3*(pi/2-up-x) * log(sin(x)), 0, pi/2-up)$value +
integrate(\(x) -5*(x-up) * log(cos(x)) - 3*(x-up) * log(sin(x)), 0, up)$value +
	+ 3/4 * pi^2*log(2) - 2*up*pi*log(2) + 4*up^2*log(2);


# Helper: on [0, pi/2-up] * [up, pi/2-up]
up = 3*pi/7;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(x+y)), 0, pi/2-up)$value), up, pi/2-up)
integrate(\(x) 5*(pi/2-up-x) * log(cos(x)) + 3*(pi/2-up-x) * log(sin(x)), 0, pi/2-up)$value +
	+ 1/2 * (pi/2-up)*pi*log(2) + 2*(pi/2-up)^2*log(2);


### I( log(sin(x-y)) )
up = 5*pi/7; # up = pi/4;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(y-x)), 0, y)$value), 0, up)
integrate(\(x) sapply(x, \(y) integrate(\(x) x / tan(y-x) - y/(y-x), 0, y)$value), 0, up)$value +
	+ up^2*log(up)/2 - up^2/4;
# TODO

### on [0, pi/4]
integrate(\(x) sapply(x, \(y) integrate(\(x) log(sin(y-x)), 0, y)$value), 0, pi/4)
- (35/4 * pracma::zeta(3) + pi^2 * log(2)) / 32;

### on [0, pi/3]
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(sin(y-x)), 0, y, rel.tol=1E-13)$value), 0, pi/3, rel.tol=1E-13)
- (13/2 * pracma::zeta(3) + pi^2 * log(2)) / 18;


### Log( TAN )

### I( log(tan(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(x+y)), 0, pi/4)$value), 0, pi/4)
0


### Prod: x*y

### I( log(tan(pi/2*x*y)) )
# Note: as split integral;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(pi/2*x*y)), 0, 1)$value), 0, 1/2)$value +
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(tan(pi/2*x*y)), 0, 1, rel.tol=1E-12)$value), 1/2, 1, rel.tol=1E-12)$value
# TODO


### I( log(tan(pi/4*x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(pi/4*x*y)), 0, 1)$value), 0, 1)
# TODO

