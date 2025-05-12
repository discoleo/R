#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
##
## v.0.1b

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
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1-x*y), 0, 1)$value), 0, 1)
pi^2 / 6

### I( 1 / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) 1 / (1+x*y), 0, 1)$value), 0, 1)
pi^2/12


### Radicals

### I( (x^p + y^p) / (1-x*y) / (x*y)^q )
# Maths 505: A surprisingly interesting integral
# https://www.youtube.com/watch?v=BNzVNcZrbfA
# Note: Series expansion of 1/(1 - x*y)

integrate(\(x) sapply(x, \(y)
	integrate(\(x) (sqrt(x) + sqrt(y)) / (1-x*y) / (x*y)^(1/4), 0, 1)$value), 0, 1)
(digamma(5/4) - digamma(3/4)) * 4


### Gen: I( (x^p + y^p) / (1-x*y) / (x*y)^q )
p = sqrt(3); q = sqrt(5) - 2;
integrate(\(x) sapply(x, \(y)
	integrate(\(x) (x^p + y^p) / (1-x*y) / (x*y)^q, 0, 1)$value), 0, 1)
(digamma(1+p-q) - digamma(1-q)) * 2 / p


###########

###########
### Log ###

### Simple Logs

### I( x * log(x) * log(y)^2 / (1 - x*y) )
# Maths 505: A very interesting double integral with a beautiful result
# https://www.youtube.com/watch?v=lWVBnWCz63I
# - series expansion of 1/(1 - x*y);

integrate(\(x) sapply(x, \(y)
	integrate(\(x) x * log(x) * log(y)^2 / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^3 * 1/(j+1)^2); j >= 1;
- 2*pracma::zeta(3) + 6*pracma::zeta(2) - 8;


### I( x * log(x)^2 * log(y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x * log(x)^2 * log(y) / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^2 * 1/(j+1)^3); j >= 1;
- 2*pracma::zeta(3) - 6*pracma::zeta(2) + 12;


### I( x^2 * log(x)^2 * log(y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x^2 * log(x)^2 * log(y) / (1 - x*y), 0, 1)$value), 0, 1)
# - 2 * sum(1/j^2 * 1/(j+2)^3); j >= 1;
- pracma::zeta(3)/2 - 3/4*pracma::zeta(2) + 3/8 * (1+1/2) + 1/2*(2+1/4+1/8);
- pracma::zeta(3)/2 - 3/4*pracma::zeta(2) + 1 + 3/4;


# TODO: more variants;


### I( log(x)*log(y) / (1 - x*y)^3 )
# Maths 505: A lovely surprise awaits at the end of the video
# https://www.youtube.com/watch?v=24CCdCBF5ck
# Series expansion of 1/(1 - x*y)^3;

integrate(\(x) sapply(x, \(y) integrate(\(x) log(x)*log(y) / (1 - x*y)^3, 0, 1)$value), 0, 1)
(pracma::zeta(2) + pracma::zeta(3)) / 2


### Mixed Logs

### I( log(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y), 0, 1)$value), 0, 1)
pi^2/6 - 2

### I( log(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y), 0, 1)$value), 0, 1)
pi^2/12 + 2*log(2) - 2


### I( log(1 + x + x*y) )
# Maths 505: 
# https://www.youtube.com/watch?v=xiSw7Ou4RIY

### I( log(1 + x + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 + x + x*y), 0, 1)$value), 0, 1)
3*log(3) - 2*log(2) + polylog2(-1) - polylog2(-2) - 2;
3*log(3) - 2*log(2) - pi^2/12 - polylog2(-2) - 2;


### I( log(1 - x + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 - x + x*y), 0, 1)$value), 0, 1)
pi^2/6 - 2


### I(x, y, z)
# Fails numerically: but trivial extension of above integral;
# library(Rmpfr)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) sapply(x, \(z)
	integrate(\(x) {
		x = mpfr(x, 240); y = mpfr(y, 240);
		as.numeric(log(x) + log(1 + y + y*z));
	}, 0, 1, rel.tol=1E-12)$value),
	0, 1, rel.tol=1E-10)$value), 0, 1)
3*log(3) - 2*log(2) - pi^2/12 - polylog2(-2) - 3;


### I( log(1 - x^2*y^2) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x^2*y^2), 0, 1)$value), 0, 1)
pi^2/4 + 2*log(2) - 4


### Log-Fractions

### I( log(1-x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1-x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(1+x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1+x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (1 - x*y), 0, 1)$value), 0, 1)

# I( log((x+y)/2) / (1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log((x+y)/2) / (1 - x*y), 0, 1)$value), 0, 1)


### I( log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (1 + x*y), 0, 1)$value), 0, 1)


### I( (x+y) * log(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1+x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
- pi^2/2 + 2*log(2)^2 + 4;


### I( (x+y) * log(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) {
	integrate(\(x) (x+y)*log(x+y) / (1-x*y), 0, 1, rel.tol=1E-12)$value
	}), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y) integrate(\(x) 2 * log(1-x*y) / (x+y),
	0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)$value +
	- 2*log(2)^2 + 4*log(2) + pi^2/3 - 4;
# TODO


### I( log(x+y+1) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y+1) / (1 + x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(x+y+1) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y+1) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(1 + x*y) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 + x*y) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(1 - x*y) / (x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(1 - x*y) / (x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)


### I( log(x+y) / (1 + x^2 + y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y) / (1 + x^2 + y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x+y) / (2 - x^2 - y^2) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x+y) / (2 - x^2 - y^2), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### Pow = 2 (inside Log)

### I( log(x^2+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2-x*y+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2-x*y+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2-x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2+x*y+y^2) / (1 - x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+x*y+y^2) / (1-x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### I( log(x^2+x*y+y^2) / (1 + x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x^2+x*y+y^2) / (1+x*y), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
# TODO


### Prod( LOG )

### I( log(x) * log(y) * log(1 - x*y) )
# Maths 505: A mesmerizing result
# https://www.youtube.com/watch?v=QqVhd_xfjnc
# Note: series expansion of log(1 - x*y);

integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(x) * log(y) * log(1 - x*y), 0, 1, rel.tol=1E-12)$value
	), 0, 1, rel.tol=1E-12)
pracma::zeta(2) + pracma::zeta(3) + pracma::zeta(4) - 4


### I( log(x) * log(y) * log(z) * log(1 - x*y*z) )
# Note: rather slow (5 - 10 s);
integrate(\(x) sapply(x, \(y)
	integrate(\(x) sapply(x, \(z)
	integrate(\(x) log(x) * log(y) * log(z) * log(1 - x*y*z),
		0, 1, rel.tol=1E-12)$value
	), 0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
6 - sum(pracma::zeta(2:6));


###############
###############

############
### ATAN ###

### I( atan(x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y), 0, 1)$value), 0, 1)
- pi^2 / 48 + pi/4 - log(2)/2


### I( atan(1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 - x*y), 0, 1)$value), 0, 1)
- (5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan)


### I( atan(1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(1 + x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x*y) / (x + y + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (x+y+1), 0, 1)$value), 0, 1)
# TODO

### I( atan(x*y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 + x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x*y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x*y) / (1 - x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x/y) / (x*y+1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x/y) / (1+x*y), 0, 1)$value), 0, 1)
pi^3 / 48


### I( atan(x+y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (x + y), 0, 1)$value), 0, 1)

### I( atan(x+y) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1-x*y), 0, 1)$value), 0, 1)
# TODO

### I( atan(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) atan(x+y) / (1+x*y), 0, 1)$value), 0, 1)
# TODO


##############
##############

###########
### Exp ###

### I( log(exp(x) + exp(y)) )
# Maths 505: An extremely captivating double integral
# https://www.youtube.com/watch?v=etR1AI3anpU
# - for polylog2, see file: Integrals.Trig.Tan.R;

integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x) + exp(y)), 0, 1)$value), 0, 1)
- 3/2*pracma::zeta(3) - pi^2/6 - 2*polylog2(-exp(1), 3) + 1/3;


### I( log(exp(1) - exp(x*y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(1) - exp(x*y)), 0, 1)$value), 0, 1)
# TODO

###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(exp(x*y) - 1), 0, 1)$value), 0, 1)
# TODO


### I( 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x) + exp(y))) )
# Maths 505: This double integral will make you love calculus
# https://www.youtube.com/watch?v=-Nm5BofLiS4


### on [0, Inf]^2
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x) + exp(y))),
		0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((x+1) * (y+1) * (1/x + 1/y)),
		0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12);
2*log(2) - log(2)^2 - pi^2/12;


### on [0,1]^2
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((exp(x)+1) * (exp(y)+1) * (exp(x)+exp(y))),
		0, 1, rel.tol=1E-12)$value), 0, 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) 1 / ((x+1) * (y+1) * (1/x + 1/y)),
		exp(-1), 1, rel.tol=1E-12)$value), exp(-1), 1, rel.tol=1E-12)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) x / ((x+1) * (x+y)),
		exp(-1), 1, rel.tol=1E-12)$value), exp(-1), 1, rel.tol=1E-12)$value +
	- log(2/(1+exp(-1)))^2 / 2;
integrate(\(x) - x * log(x + exp(-1)) / (x+1), exp(-1), 1)$value +
2*log(2) - log(2)^2 + exp(-1) - 1 +
	- (1+exp(-1) - log(2))*log(1+exp(-1));
polylog2(-1 - 2/(exp(1)-1)) - polylog2(-2/(exp(1)-1)) +
	+ log(exp(2)-1) * log(2) - 2*log(exp(1)+1)*exp(-1) +
	- log(exp(1)-1) * log(exp(1)+1) +
	+ 2*log(2)*exp(-1) - log(2)^2 + 1;
# TODO: closed formula for Li2?


# Helper:
integrate(\(x) x * log(x + 1) / (x+1), exp(-1), 1)
2*log(2) - log(2)^2/2 - 1 - (exp(-1)+1)*log(exp(-1)+1) +
	+ exp(-1) + log(exp(-1)+1)^2/2

# Note:
# - for function polylog2, see file:
#   Integrals.Trig.Tan.R;
integrate(\(x) x * log(x + exp(-1)) / (x+1), exp(-1), 1)
polylog2(-2/(exp(1)-1)) - polylog2(-1 - 2/(exp(1)-1)) +
	+ 2*(1 - log(2))*exp(-1) - 1 + log(exp(1)+1)*exp(-1) +
	+ (log(exp(1)-1) - 1) * log((exp(1)+1)/2);


###############

### I( log(x+y) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO

### I( log(1 - x*y) / (exp(x*y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1 - x*y) / (exp(x*y) - 1), 0, 1)$value), 0, 1)
# TODO

### I( log(1 - x*y) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(1 - x*y) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO


##############

### Other Exp/Base:

### I( log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y)) )
# Dr Peyam: A surprisingly elegant double integral
# https://www.youtube.com/watch?v=q17ezsqKRis
# Note: => triple integral:
# I( 1 / 2^(x+y) / (z*(2^x + 2^y) + 1) ) with z on [0,1] =>
# I( x*y / (x*y + x*z + y*z) ) / log(2)^2 on [0,1]^3;


# Note: numeric instability: upper = Inf;
integrate(\(x) sapply(x, \(y)
	integrate(\(x) log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y)),
		0, 100, rel.tol=1E-10)$value), 0, 100, rel.tol=1E-8)
1/3 / log(2)^2


# library(Rmpfr)
integrate(\(x) sapply(x, \(y)
	integrate(\(x) {
		x = mpfr(x, 240); y = mpfr(y, 240);
		v = log(2^x + 2^y + 1) / (2^(x+y) * (2^x + 2^y));
		as.numeric(v);
		}, 0, Inf, rel.tol=1E-10)$value), 0, Inf, rel.tol=1E-10)


############

############
### Trig ###

### I( sin(pi/2*(x+y)) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (x + y), 0, 1)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(\(x) 4/pi*sin(x)*cos(y) / (x + y), 0, pi/2)$value), 0, pi/2)
# TODO


### I( sin(pi*(x+y)) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*(x+y)) / (x + y), 0, 1)$value), 0, 1)


### I( sin(pi/2*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)

### I( sin(pi*(x+y)) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*(x+y)) / (1 + x*y), 0, 1)$value), 0, 1)

### I( sin(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(x+y) / (1 + x*y), 0, pi/2)$value), 0, pi/2)


### I( cos(x+y) / (1 + x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) cos(x+y) / (1 + x*y), 0, pi/2)$value), 0, pi/2)


### I( sin(pi/2*(x+y)) / (1 - x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*(x+y)) / (1 - x*y), 0, 1)$value), 0, 1)


### Trig( Prod )

### I( sin(pi/2*x*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi/2*x*y) / (x + y), 0, 1)$value), 0, 1)

### I( sin(pi*x*y) / (x + y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) sin(pi*x*y) / (x + y), 0, 1)$value), 0, 1)



### Trig-Fractions

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
- (2*pi*log(2) - pi - 4*Catalan + 2*log(2) + 2) / 16



###
integrate(\(x) sapply(x, \(y) integrate(
	\(x) sin(pi/2*x*y) / (cos(pi/2*x) + cos(pi/2*y)), 0, 1)$value), 0, 1)
# TODO


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
7/8 * pracma::zeta(3) + (pi/2)^2 * log(2)


### I( log(cos(x-y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(cos(x-y)), 0, pi/2, rel.tol=1E-12)$value), 0, pi/2, rel.tol=1E-12)
7/8 * pracma::zeta(3) - (pi/2)^2 * log(2)


### I( log(cos(x) + cos(y)) )
# see next Section: cos(x-y)*cos(x+y);
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x) + cos(y)), 0, pi/2)$value), 0, pi/2)
7/4 * pracma::zeta(3) - (pi/2)^2 * log(2);


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


### on [0, pi/3]^2
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/3)$value), 0, pi/3)
integrate(\(x) 2*(pi/3 - x) * log(cos(x)), 0, pi/3);
id = 1:2; ic = 1:3; sn = sin(2*pi*id/6); cs = cos(2*pi*ic/6);
-(pi/3)^2 * log(2) + 3/8 * pracma::zeta(3) +
	- sum(cs * (pracma::psi(2, ic/6) - pracma::psi(2, 1/2 + ic/6))) / (4*6^3);

# Note:
# - for sub-integrals, see: Integrals.Log.Trig.R;
# - alternatively: Clausen function;

# TODO: simplify the sum?


### on [0, pi/8]^2
n = 8; # EVEN: 4*k;
integrate(\(x) sapply(x, \(y) integrate(\(x) log(cos(x-y)), 0, pi/n)$value), 0, pi/n)
integrate(\(x) 2*(pi/n - x) * log(cos(x)), 0, pi/n);
id = seq(n/2 - 1); id2 = seq(n/4 - 1); n2 = n/2;
sn = sin(2*pi*id/n); sn2 = sin(2*pi*id2/n2);
cs = cos(2*pi*id/n); cs2 = cos(2*pi*id2/n2);
dd = 512; sg = c(1,-1,1); # NOT generalized;
2*pi/n * (- pi*log(2) / 8 + sum(
		+ sg*sn * (pracma::psi(1, id/16) - pracma::psi(1, 1 - id/16)) +
		- sg*sn * (pracma::psi(1, 1/2 - id/16) - pracma::psi(1, 1/2 + id/16)) ) / dd) +
	- 1/2 * (pracma::zeta(3) * (1/2 + 3/n2^3)/2 - (pi/n2)^2 * log(2)/2 +
	+ sum(cs2 * (pracma::psi(2, id2/n2) + pracma::psi(2, 1 - id2/n2))) / (8*n2^3) +
	- sum(sn2 * (pracma::psi(1, id2/n2) - pracma::psi(1, 1 - id2/n2))
		) * pi / (2*n2^3) ) + log(2)*(pi/n)^2 +
	+ 2*(pracma::zeta(3) * (1/2 + 3/n^3)/2 - (pi/n)^2 * log(2)/2 +
		+ sum(cs * (pracma::psi(2, id/n) + pracma::psi(2, 1 - id/n))) / (8*n^3) +
		- sum(sn * (pracma::psi(1, id/n) - pracma::psi(1, 1 - id/n))
		) * pi / (2*n^3) );

# TODO: generalize & simplify;


### I( log(tan(x+y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(x+y)), 0, pi/4)$value), 0, pi/4)
0


### Prod:

###
integrate(\(x) sapply(x, \(y) integrate(\(x) log(tan(pi/4*x*y)), 0, 1)$value), 0, 1)
# TODO

