##################
##
## Leonard Mada
##
## Integrals: Trig
## Variants: ATAN


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

#####################
#####################

### Basic Integrals

### I( atan(x) / x )
integrate(\(x) atan(x) / x, 0, 1)
Catalan

# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


### I( atan(x) / x^2 )
integrate(\(x) - (atan(x) / x^2 - 1/x), 0, 1)
pi/4 + log(2)/2 - 1


####################
####################

### Fractions

# Note:
# - see also file: Integrals.Log.Fractions.Complex.R;


### I( atan(x^2) / (x^2 + 1) )
integrate(function(x) atan(x^2) / (x^2 + 1), 0, Inf)
pi^2 / 8

### I( atan(x^n) / (x^2 + 1) )
n = sqrt(3) # independent of n;
integrate(function(x) atan(x^n) / (x^2 + 1), 0, Inf)
pi^2 / 8


### I( atan(x^2) / (x^4 + 1) )
integrate(\(x) atan(x^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(pi/4) / 8 - pi * (digamma(1/2) - digamma(1/4)) / sin(pi/4) / 8
pi*(pi - 2*log(2)) * sqrt(2) / 16


### Gen: I( atan(k * x^2) / (x^4 + 1) )
k = 5^(1/3)
integrate(\(x) atan(k*x^2) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
x = 1/k^(1/2);
pi^2 / sin(pi/4) / 8 +
	+ (log(x^2 + 1) - 2*log(x + 1) - 2*atan(x)) * pi / sin(pi/4) / 8;

# Solution: Explicit formula
id = c(2,4) * pi/4; x = 1/k^(1/2);
cs = cos(id); cs1 = cs[1]; csp = cos(pi);
sn = sin(id); sn1 = sn[1]; snp = sin(pi);
pi^2 / sin(pi/4) / 8 +
- (sum(csp*log(x^2 - 2*cs1*x + 1) +
	- 2*snp * (atan((x - cs1)/sn1) + atan(cs1/sn1))) +
- (sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))))) * pi / sin(pi/4) / 8;



### ATAN-Pow: 4

### I( atan(x^4) / (x^4 + 1) )
integrate(\(x) atan(x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(pi/4) / 8 - pi * (
	(digamma(13/16) - digamma(5/16)) / sin(pi/8) +
	(digamma(9/16) - digamma(1/16)) / sin(5*pi/8) +
	- (digamma(3/4) - digamma(1/4)) / sin(pi/4) * 2) / 32;


### I( x^2 * atan(x^4) / (x^4 + 1) )
integrate(\(x) x^2 * atan(x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
pi^2 / sin(3*pi/4) / 8 - pi * (
	(digamma(15/16) - digamma(7/16)) / sin(3*pi/8) +
	(digamma(11/16) - digamma(3/16)) / sin(7*pi/8) +
	- (digamma(3/4) - digamma(1/4)) / sin(3*pi/4) * 2) / 32;


### Gen: I( atan(k * x^4) / (x^4 + 1) )
# see Integrals.Fractions.Unity.R;
k = 3^(1/4)
integrate(\(x) atan(k * x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
# Solution:
pfr = seq(1, 8, by=2) * pi / 8; cs = cos(pfr); sn = sin(pfr);
cs4 = cos(5*pfr); sn4 = sin(5*pfr); x = 1 / k^(1/4);
cs3 = cos(4*pfr); sn3 = sin(4*pfr);
pi^2 / sin(pi/4) / 8 - pi*(
	+ sum(cs4*log(x^2 + 2*cs*x + 1) +
		+ 2*sn4 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(pi/8) +
	+ sum(cs*log(x^2 + 2*cs*x + 1) +
		+ 2*sn * (atan((x + cs)/sn) - atan(cs/sn))) / sin(5*pi/8) +
	+ sum(cs3*log(x^2 + 2*cs*x + 1) +
		+ 2*sn3 * (atan((x + cs)/sn) - atan(cs/sn))) * 2 / sin(pi/4)) / 16;


### Gen: I( x^2 * atan(k * x^4) / (x^4 + 1) )
k = 3^(1/4)
integrate(\(x) x^2 * atan(k * x^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
# Solution:
pfr = seq(1, 8, by=2) * pi / 8; cs = cos(pfr); sn = sin(pfr);
cs6 = cos(7*pfr); sn6 = sin(7*pfr); x = 1 / k^(1/4);
cs3 = cos(4*pfr); sn3 = sin(4*pfr);
cs2 = cos(3*pfr); sn2 = sin(3*pfr);
pi^2 / sin(pi/4) / 8 - pi*(
	+ sum(cs6*log(x^2 + 2*cs*x + 1) +
		+ 2*sn6 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(3*pi/8) +
	+ sum(cs2*log(x^2 + 2*cs*x + 1) +
		+ 2*sn2 * (atan((x + cs)/sn) - atan(cs/sn))) / sin(7*pi/8) +
	+ sum(cs3*log(x^2 + 2*cs*x + 1) +
		+ 2*sn3 * (atan((x + cs)/sn) - atan(cs/sn))) * 2 / sin(3*pi/4)) / 16;

# [alternative] Compact formula:
pi^2 / sin(pi/4) / 8 - pi*(
	sum((cs6/sin(3*pi/8) + cs2/sin(7*pi/8) + 2*cs3/sin(3*pi/4)) *
		log(x^2 + 2*cs*x + 1)) / 16 +
	sum((sn6/sin(3*pi/8) + sn2/sin(7*pi/8) + 2*sn3/sin(3*pi/4)) *
		(atan((x + cs)/sn) - atan(cs/sn))) / 8);


# Derivation:

# from atan(x^2 / b^2) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 2*b
integrate(\(x) x^2/(x^4 + b^4) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) x^2 * (1/(x^4+1) - 1/(x^4 + b^4)) / (b^4 - 1), 0, Inf)
pi * (1/sin(3*pi/4) - 1/sin(pi/4) / b) / 4 / (b^4 - 1)


# from atan(x^4 / b^4) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 4*b^3
integrate(\(x) x^4/(x^8 + b^8) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) ((x^4+b^8)/(x^8 + b^8) - 1/(x^4+1)) / (b^8 + 1), 0, Inf)
pi * (1/sin(pi/8) * b + 1/sin(5*pi/8) / b^3 - 2/sin(pi/4)) / 8 / (b^8+1)

# from x^2 * atan(x^4 / b^4) / (x^4 + 1)
b = 5^(1/8) # does NOT include factor * 4*b^3
integrate(\(x) x^6/(x^8 + b^8) / (x^4 + 1), 0, Inf, rel.tol=1E-9)
integrate(\(x) x^2 * ((x^4+b^8)/(x^8 + b^8) - 1/(x^4+1)) / (b^8 + 1), 0, Inf)
pi * (1/sin(3*pi/8) * b^3 + 1/sin(7*pi/8) / b - 2/sin(3*pi/4)) / 8 / (b^8+1)

# Note:
# (digamma(3/4) - digamma(1/4)) == pi;
# (allows slight simplification)


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


###
integrate(\(x) atan(x) * atan(1-x), 0, 1)
integrate(\(x) (1-2*x) * atan(1-x) / (x^2+1), 0, 1)
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
integrate(function(x) atan(x) * log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


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

