#####################
##
## Leonard Mada
## (the one and only)
##
## Integrals: Double Integrals
## Exp-Type
##
## v.0.1a

### Double Integrals

### Examples:
# I( x^p * y^q / (exp(x+y) - 1) )
# I( x^p * y^q / (exp(x+y) + 1) )
# I( log(exp(x) + exp(y)) )


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;
dzeta2  = - 0.937548254316;

#####################
#####################

###########
### Exp ###

### Fractions

### I( exp(x+y) / (1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(-x-y) / (1+x*y), 0, Inf)$value), 0, Inf)
# TODO


### on [0,1]^2

### I( exp(x+y) / (1+x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x+y) / (1+x*y), 0, 1)$value), 0, 1)

### I( exp(x+y) / (1-x*y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x+y) / (1-x*y), 0, 1)$value), 0, 1)

### I( (exp(1) - exp((x+y)/2)) / (1-x*y) )
integrate(\(x) sapply(x, \(y)
	integrate(\(x) (exp(1) - exp((x+y)/2)) / (1-x*y), 0, 1)$value), 0, 1)
# TODO


### Exp( x*y )

### I( exp(x*y) / (x+y) )
integrate(\(x) sapply(x, \(y) integrate(\(x) exp(x*y) / (x+y), 0, 1)$value), 0, 1)
# TODO


### Div: Exp

### I( (x+y) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
2 * pracma::zeta(3);

### I( x*y / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) x*y / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
pracma::zeta(4); pi^4/90;


### Gen: I( (x+y) / (exp(k*(x+y)) - 1) )
k = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(k*(x+y)) - 1), 0, Inf)$value), 0, Inf)
2 * pracma::zeta(3) / k^3;

### Gen: I( x^p / (exp(x+y) - 1) )
p = sqrt(2);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^p / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
gamma(p+1) * pracma::zeta(p+2);

### Gen: I( x^p * y^q / (exp(x+y) - 1) )
px = sqrt(2); py = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(x+y) - 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2);

### Gen: I( x^p * y^q / (exp(k1*x + k2*y) - 1) )
px = sqrt(2); py = sqrt(3);
k1 = 1 / sqrt(3); k2 = sqrt(5) - sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(k1*x+k2*y) - 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2) / k1^(px+1) / k2^(py+1);


### I( (x+y) / (exp(x+y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(x+y) + 1), 0, Inf)$value), 0, Inf)
3/2 * pracma::zeta(3)

### Gen: I( (x+y) / (exp(k*(x+y)) + 1) )
k = sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) (x+y) / (exp(k*(x+y)) + 1), 0, Inf)$value), 0, Inf)
3/2 * pracma::zeta(3) / k^3;

### Gen: I( x^p * y^q / (exp(k1*x + k2*y) + 1) )
px = sqrt(2); py = sqrt(3);
k1 = 1 / sqrt(3); k2 = sqrt(5) - sqrt(3);
integrate(\(x) sapply(x, \(y) integrate(\(x) x^px * y^py / (exp(k1*x+k2*y) + 1), 0, Inf)$value), 0, Inf)
gamma(px+1) * gamma(py+1) * pracma::zeta(px+py+2) / k1^(px+1) / k2^(py+1) * (1 - 1/2^(px+py+1));


### Other

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


#####################

### Mixed: Log( Exp )

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


####################

### Mixed: Log & Exp

### I( log(x) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x) / (exp(x+y) - 1), 0, Inf, rel.tol=1E-8)$value), 0, Inf, rel.tol=1E-6)
# gamma(1) * digamma(1) * zeta(2) + gamma(1) * dzeta(2);
- Euler * pi^2 / 6 + dzeta2;


### I( log(x+y) / (exp(x+y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (exp(x+y) - 1), 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-11)
pi^2/6 - Euler*pi^2/6 + dzeta2;


### I( log(x+y) / (exp(x+y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y) / (exp(x+y) + 1), 0, Inf, rel.tol=1E-12)$value), 0, Inf, rel.tol=1E-11)
(pi^2/6 * (log(2) - Euler + 1) + dzeta2) / 2;


### EXP(x*y)

### I( log(x+y) / (exp(x*y) - 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	log(x+y)*(1/(exp(x*y) - 1) - 1/(x*y)), 0, 1)$value), 0, 1)
# TODO


### I( log(x+y) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
# TODO

### I( log(x+y+1) / (exp(x*y) + 1) )
integrate(\(x) sapply(x, \(y) integrate(\(x) log(x+y+1) / (exp(x*y) + 1), 0, 1)$value), 0, 1)
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


##################
##################

### Exp & Trig

### I( sin(y) * exp(sin(y) * (cos(x) - sin(x))) )
# Maths 505: Impossible integral? Solution using symmetry
# https://www.youtube.com/watch?v=OkuYTH0S_dM
# Note: Inverse of Spherical coords;

integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(y) * exp(sin(y) * (cos(x) - sin(x))), 0, 2*pi, rel.tol=1E-12)$value), 0, pi, rel.tol=1E-12)
pi * sinh(sqrt(2)) * sqrt(2) * 2;


### Gen: I( sin(y) * exp(k * sin(x)*sin(y)) )
k = 1/sqrt(5);
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(y) * exp(k * sin(x)*sin(y)), 0, 2*pi, rel.tol=1E-12)$value), 0, pi, rel.tol=1E-12)
4*pi * sinh(k) / k;


### I( sin(y) * exp(sin(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(y) * exp(sin(x)*sin(y)), 0, 2*pi, rel.tol=1E-12)$value), 0, pi, rel.tol=1E-12)
4*pi * sinh(1);


### I( sin(y) * exp(cos(x)*sin(y)) )
integrate(\(x) sapply(x, \(y) integrate(\(x)
	sin(y) * exp(cos(x)*sin(y)), 0, 2*pi, rel.tol=1E-12)$value), 0, pi, rel.tol=1E-12)
4*pi * sinh(1);

