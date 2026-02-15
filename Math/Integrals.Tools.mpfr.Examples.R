


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
f2 = \(x) log(1 + log(1+x))
integrate(f2, 0, 1)
int.mpfr(f2, 0, 1)
integrate(f2, 0, 1)$value +
	- int.mpfr(f2, 0, 1, p=p61, w = list(w30, w61))$value;


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
3 - 2*Catalan;

