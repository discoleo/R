


####################

### Helper Functions

source("Integrals.Tools.mpfr.R")

### Helper Constants
Euler   = mpfr("0.577215664901532860606512090", 144);
Catalan = mpfr("0.915965594177219015054603514", 144);

####################

### Examples:


### Example 1:
f1 = \(x) exp(x^2)
integrate(f1, 0, 1)
int.mpfr(f1, 0, 1)


### Example 2:
f2  = \(x) log(1 + log(1+x));
f2e = \(x) log(1+x) * exp(x);
integrate(f2, 0, 1)
int.mpfr(f2, 0, 1)
int.mpfr(f2e, 0, log(mpfr(2, 128)))
# Error Quantitation:
integrate(f2, 0, 1)$value - int2.mpfr(f2, 0, 1)$value;
int2.mpfr(f2, 0, 1)$value - int.mpfr(f2, 0, 1)$value


### Example 3:
# Still fails:
FUN = \(x, y) sqrt( y/x * (1-x) / (1 - x*y) );
int.mpfr(\(x) sapplyMpfr(x, \(y) int.mpfr(FUN, 0, 1, y=y)$value), 0, 1)
3 - 2*Catalan;

#
integrate(\(x) as.numeric(sapplyMpfr(x, \(y) int.mpfr(FUN,
	0, 1, y=y, p=p61, w = list(w30, w61))$value)), 0, 1, rel.tol=1E-9)
3 - 2*Catalan;

#
integrate(\(x) as.numeric(sapplyMpfr(x, \(y) int.mpfr(FUN,
	0, 1/4, y=y, p=p61, w = list(w30, w61))$value)), 0, 1, rel.tol=1E-9)$value +
integrate(\(x) as.numeric(sapplyMpfr(x, \(y) int.mpfr(FUN,
	1/4, 1, y=y, p=p61, w = list(w30, w61))$value)), 0, 1, rel.tol=1E-9)$value;
3 - 2*Catalan;

#
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, 1/2, y=y, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)$value +
integrate(\(x) sapply(x, \(y) integrate(FUN, 1/2, 1, y=y, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)$value;
3 - 2*Catalan;
#
FUN = \(x, y) sqrt( x * (1-y) / (1 - x)) / y^2;
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, y/2, y=y, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)$value +
integrate(\(x) sapply(x, \(y) integrate(FUN, y/2, y, y=y, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9)$value;
#
FUN = \(x, y) 2 * sin(x)^2 * sqrt(1-y)/y^2;
FUN = \(x, y) (1 - cos(2*x)) * sqrt(1-y)/y^2;
integrate(\(x) sapply(x, \(y) integrate(FUN, 0, asin(sqrt(y)), y=y, rel.tol=1E-9)$value), 0, 1, rel.tol=1E-9);
#
integrate(\(x)  (asin(sqrt(x)) - sqrt(x*(1-x))) * sqrt(1-x)/x^2, 0, 1);
integrate(\(x)  2*(asin(x) - x*sqrt(1-x^2)) * sqrt(1-x^2)/x^3, 0, 1);
integrate(\(x)  2*(x - sin(x)*cos(x)) * cos(x)^2/sin(x)^3, 0, pi/2);
integrate(\(x)  2*(x - sin(x)) * cos(x)^2/sin(x)^3, 0, pi/2)$value + 4 - pi;
3 - 2*Catalan;


# Helper:
integrate(\(x) (cos(x) - 1) * cos(x)^3/sin(x)^2, 0, pi/2);
integrate(\(x) - sin(x) * (1+sin(x)^2) / sin(x), 0, pi/2)$value + 2;
2 - 3/4 * pi;
# =>
integrate(\(x)  2*(sin(x) - sin(x)*cos(x)) * cos(x)^2/sin(x)^3, 0, pi/2);
integrate(\(x)  2*(sin(x)^2 + cos(x)^2 - cos(x)) * cos(x)^2/sin(x)^2, 0, pi/2);
pi/2 + 2*(2 - 3/4 * pi); 4 - pi;


###################
###################

### I( x * sqrt( (1-x)/(1+y) / (1 - x^2*y^2) ) )
# Precision: only 1E-8
FUN = \(x, y) x * sqrt( (1-x)/(1+y) / (1 - x^2*y^2) );
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN, 0, 1, y=y)$value), 0, 1)

#
FUN = \(x, y) x * sqrt( (1-x)/(1+y) / (1 - x^2*y^2) );
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN, 0, 1, y=y)$value), 0, 1)$value
#
FUN2 = \(x, y) x * sqrt( (1-x/y)/(1+y) / (1 - x^2) ) / y^2;
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN2, 0, y, y=y)$value), 0, 1)
# FAILS at 1E-7;
FUN2 = \(x, y) sin(x) * sqrt( (1-sin(x)/y)/(1+y) ) / y^2;
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN2, 0, asin(y), y=y)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(FUN2, 0, asin(y), y=y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)

# Alternative:
# FAILS at 1E-7;
FUN2 = \(x, y) sqrt( y*(1-y) / (x+y) / (1 - x^2) );
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN2, 0, y, y=y)$value), 0, 1)
integrate(\(x) sapply(x, \(y) integrate(FUN2, 0, y, y=y, rel.tol=1E-13)$value), 0, 1, rel.tol=1E-13)
#
FUN2 = \(x, y) sqrt( y*(1-y) / (sin(x) + y) );
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN2, 0, asin(y), y=y)$value), 0, 1)
# TODO: ???

# Alternative:
FUN2 = \(x, y) y * sqrt( x*(1-x*y) / (1+x) / (1 - y^2) );
int2.mpfr(\(x) sapplyMpfr(x, \(y) int2.mpfr(FUN2, 1, 1/y, y=y)$value), 0, 1)

