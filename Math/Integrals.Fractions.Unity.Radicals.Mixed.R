########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration:
##   Polynomial Fractions:
##   Mixed Radicals
##
## draft v.0.1i


### Mixed Radicals

### Examples
# I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
# I( x^2 * sqrt(1 + x^2) / (x^4 + 1) )
# I( x^3 * (1 + x^3)^(2/3) / (x^6 + 1) )
# I( x^3 * (1 - x^3)^(2/3) / (x^6 + 1) )
# I( x^4 * (1 - x^4)^(3/4) / (x^8 + 1) )
# I( x^5 * (1 - x^5)^(4/5) / (x^10 + 1) )
# I( x^6 * (1 - x^5)^(3/5) / (x^5 + 1) )
# I( x^7 * (1 - x^5)^(2/5) / (x^5 + 1) )
# I( x^8 * (1 - x^5)^(1/5) / (x^5 + 1) )
### Pow: 3*n
# I( x^3 * (1 + x^3)^(2/3) / (x^9 + 1) )
# I( x^5 * (1 + x^5)^(4/5) / (x^15 + 1) )
### Simple:
# I( 1 / (x^n + 1)^(1/n) )


### History

# - moved "Mixed Radicals" to this file;
# TODO: move the remaining section;


########################
########################

### I( x^2 * sqrt(1 + x^2) / (x^4 + 1) )
# Maths 505: Why this integral is MUCH harder than it looks
# https://www.youtube.com/watch?v=vPWTmwjYM8s
# Initial: I( sqrt(x * (1 - x)) / (x^2 + 1) ) on [0, 1]
# Euler subst: sqrt(x - x^2) = x*y => x = 1/(1 + y^2);

### I( x^2 * sqrt(1 - x^2) / (x^4 + 1) )
integrate(\(x) x^2 * sqrt(1 - x^2) / (x^4 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(x * (1 - x)) / (x^2 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(tan(x) * (1 - tan(x))), 0, pi/4); # initial I;
pi/2 * (-1 + 2^(1/4) * sin(3*pi/8))

### I( x^2 * sqrt(1 + x^2) / (x^4 + 1) ) on [0, 1]
integrate(\(x) x^2 * sqrt(1 + x^2) / (x^4 + 1), 0, 1)
integrate(\(x) 1/2 * sqrt(x + 1) / x / (x^2 + 1), 1, Inf)
integrate(\(x) sin(x)/(1 - cos(x)^2) / cos(x)^2 / (tan(x)^4 + 1), pi/4, pi/2)
integrate(\(x) x^2 / (1 - x^2) / (2*x^4 - 2*x^2 + 1), 0, 1/sqrt(2))
integrate(\(x) x^2 / (x^2 - 1) / (x^4 - 2*x^2 + 2), sqrt(2), Inf)
integrate(\(x) 1/2 * Im((1+1i)/(x^2 - 1) - (1+1i)/(x^2 - 1 + 1i) +
	- (1-1i)/(x^2 - 1) + (1-1i)/(x^2 - 1 - 1i)), sqrt(2), Inf)
integrate(\(x) 1/(x^2 - 1) + 1/2 * Im(
	- (1+1i)/(x^2 - 1 + 1i) +
	+ (1-1i)/(x^2 - 1 - 1i) ), sqrt(2), Inf)
1/2 * log((sqrt(2)+1)/(sqrt(2)-1)) - 1/4 * Im(
	(1+1i)/sqrt(1-1i) * log((sqrt(2) + sqrt(1-1i))/(sqrt(2) - sqrt(1-1i))) +
	+ (1-1i)/sqrt(1+1i) * log((sqrt(2) - sqrt(1+1i))/(sqrt(2) + sqrt(1+1i))));


# Derivation:
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 0, 1)
integrate(\(x) 2*x^2 / ((x^2 + 1)*(x^4 + 2*x^2 + 2)), 0, Inf)
integrate(\(x) 1 / (x^2 + 1) *
	Re((1-1i)/(x^2 + 1 + 1i) + (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
integrate(\(x) Im((1-1i)/(x^2 + 1) - (1-1i)/(x^2 + 1 + 1i) +
	- (1+1i)/(x^2 + 1) + (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
integrate(\(x) -2/(x^2 + 1) +
	- Im((1-1i)/(x^2 + 1 + 1i) - (1+1i)/(x^2 + 1 - 1i)), 0, Inf)
- pi - pi/2 * Im((1-1i)/sqrt(1+1i) - ((1+1i)/sqrt(1-1i)))
- pi - pi/2 * 2^(1/4) * Im(exp(-3i*pi/8) - exp(3i*pi/8))
- pi + pi * 2^(1/4) * sin(3*pi/8)


### I( x^2 * sqrt(1 - x^2) / (x^4 + 1) ) on [0, 1/2]
integrate(\(x) x^2 * sqrt(1 - x^2) / (x^4 + 1), 0, 1/2)
integrate(\(x) 1/2 * sqrt(x * (1 - x)) / (x^2 + 1), 0, 1/4)
integrate(\(x) x^2 / ((x^2 + 1)*(x^4 + 2*x^2 + 2)), sqrt(3), Inf)
- pi/2 + pi * 2^(-3/4) * sin(3*pi/8) + atan(sqrt(3)) +
	+ 2^(-3/4) * Im(exp(-3i*pi/8) * atan(sqrt(3)/sqrt(1+1i)) +
	- exp(3i*pi/8) * atan(sqrt(3)/sqrt(1-1i)));


### I( sqrt(x * (1 - x)) / (x^2 + 1) ) on [0, 1/2]
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 0, 1/2)
integrate(\(x) sqrt(tan(x) * (1 - tan(x))), 0, atan(1/2)) # as initial I;
- pi/2 + pi * 2^(1/4) * sin(3*pi/8) +
	+ Im((1-1i)/sqrt(1+1i)*atan(1/sqrt(1+1i)) +
	- (1+1i)/sqrt(1-1i)*atan(1/sqrt(1-1i)));

###
integrate(\(x) sqrt(tan(x) * (1 - tan(x))), atan(1/2), pi/4)
integrate(\(x) sqrt(x * (1 - x)) / (x^2 + 1), 1/2, 1)
- pi/2 - Im((1-1i)/sqrt(1+1i)*atan(1/sqrt(1+1i)) +
	- (1+1i)/sqrt(1-1i)*atan(1/sqrt(1-1i)))
- pi/2 - 2^(1/4) * Im(exp(-3i*pi/8) * atan(1/sqrt(1+1i)) +
	- exp(3i*pi/8) * atan(1/sqrt(1-1i)))
- pi/2 + pi * 2^(1/4) * sin(3*pi/8) +
	- 2^(-3/4) * exp(-3i*pi/8) * log((1 + 1i*sqrt(1+1i)) / (1 - 1i*sqrt(1+1i))) +
	+ 2^(-3/4) * exp(3i*pi/8) * log((1 + 1i*sqrt(1-1i)) / (1 - 1i*sqrt(1-1i)));


##################
##################

### Radical Pow: 3

### I( x^3 * (x^3 + 1)^(2/3) / (x^6 + 1) )
integrate(\(x) x^3 * (x^3 + 1)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) x^4 / (x^3 - 1) / (x^6 - 2*x^3 + 2), 2^(1/3), Inf)
# Fraction decomposition: see below;


### I( x^3 * (1 - x^3)^(2/3) / (x^6 + 1) )
integrate(\(x) x^3 * (1 - x^3)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (x*(1-x)^2)^(1/3) / (x^2 + 1), 0, 1)
# y = x^3 / (1+x^3) => x^3 = 1/(1-y) - 1;
integrate(\(x) x^3 / (1 + x^3)^3 / ((x^3/(1+x^3))^2 + 1), 0, Inf)
integrate(\(x) x^3 / (1 + x^3) / (2*x^6 + 2*x^3 + 1), 0, Inf)
integrate(\(x) x^4 / (x^3 + 1) / (x^6 + 2*x^3 + 2), 0, Inf)
# Fraction decomposition:
bp = (1+1i)^(1/3); bn = (1-1i)^(1/3);
integrate(\(x) 1/3 * (1/(x+1) - (x+1)/(x^2 - x + 1)) +
	- 1/6 * Re(bp^2 / (x+bp) + bn^2 / (x+bn) +
	- bp^2 * (x+bp) / (x^2 - bp*x + bp^2) +
	- bn^2 * (x+bn) / (x^2 - bn*x + bn^2) ), 0, Inf)
beta(2/3,1/3) / 3 * (-1 + 1/2*(bp^2 + bn^2))
beta(2/3,1/3) / 3 * (-1 + 2^(1/3) * cos(pi/6))


# Derivation:

###
integrate(\(x) x^3 * (x^3 + 1)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (x*(x+1)^2)^(1/3) / (x^2 + 1), 0, 1)
# rev of: y = x^3 / (1-x^3) => x^3 = 1 - 1/(y+1);
integrate(\(x) x^3 / (1 - x^3)^3 / ((x^3/(1-x^3))^2 + 1), 0, 2^(-1/3))
integrate(\(x) x^3 / (1 - x^3) / (2*x^6 - 2*x^3 + 1), 0, 2^(-1/3))
integrate(\(x) x^4 / (x^3 - 1) / (x^6 - 2*x^3 + 2), 2^(1/3), Inf)
# Fraction decomposition;
integrate(\(x) 1/2 * x^4 / (x^3 - 1) *
	Im(1/(x^3 - (1+1i)) - 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) 1/2 * x^4 * Re(2/(x^3 - 1) +
	- 1/(x^3 - (1+1i)) - 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) x^4/(x^3 - 1) - 1/2 * x^4 *
	Re(1/(x^3 - (1+1i)) + 1/(x^3 - (1-1i))), 2^(1/3), Inf)
integrate(\(x) x/(x^3 - 1) - 1/2 * x *
	Re((1+1i)/(x^3 - (1+1i)) + (1-1i)/(x^3 - (1-1i))), 2^(1/3), Inf)
bp = (1+1i)^(1/3); bn = (1-1i)^(1/3);
integrate(\(x) 1/3 * (1/(x-1) - (x-1)/(x^2 + x + 1)) +
	- 1/6 * Re((1+1i)/bp / (x-bp) + (1-1i)/bn / (x-bn) +
	- (1+1i)/bp * (x-bp) / (x^2 + bp*x + bp^2) +
	- (1-1i)/bn * (x-bn) / (x^2 + bn*x + bn^2) ), 2^(1/3), Inf)
# TODO


x = sqrt(7)
x / (x^3 - 1)
x * (1/(x-1) - (x+2)/(x^2 + x + 1)) / 3
(1/(x-1) - (x-1)/(x^2 + x + 1)) / 3


x = sqrt(7); b = (sqrt(2) + 1i)^(1/3);
x / (x^3 - b^3)
x * (1/(x-b) - (x + 2*b)/(x^2 + b*x + b^2)) / (3*b^2)
(1/(x-b) - (x-b)/(x^2 + b*x + b^2)) / (3*b)
# Note: x = b*y is an alternative;


###############

### I( x^4 * (1 - x^4)^(3/4) / (x^8 + 1) )
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^8 + 1), 0, 1)
integrate(\(x) x^6 / ((x^4 + 1)*(x^8 + 2*x^4 + 2)), 0, Inf)
beta(1/4, 3/4) / 4 * (-1 + 2^(3/8) * sin(5*pi/16))


###############

### I( x^6 * (1 - x^6)^(5/6) / (x^12 + 1) )
integrate(\(x) x^6 * (1 - x^6)^(5/6) / (x^12 + 1), 0, 1)
integrate(\(x) x^10 / ((x^6 + 1)*(x^12 + 2*x^6 + 2)), 0, Inf)
beta(1/6, 5/6) / 6 * (-1 + 2^(5/12) * sin(7*pi/24))


###############

### I( x^5 * (1 - x^5)^(4/5) / (x^10 + 1) )
integrate(\(x) x^5 * (1 - x^5)^(4/5) / (x^10 + 1), 0, 1)
integrate(\(x) x^8 / (x^5 + 1) / (x^10 + 2*x^5 + 2), 0, Inf)
beta(1/5, 4/5) / 5 * (-1 + 2^(2/5) * cos(2*pi/10))


### I( x^5 * (1 + x^5)^(4/5) / (x^10 + 1) )
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^10 + 1), 0, 1)
integrate(\(x) x^8 / (x^5 - 1) / (x^10 - 2*x^5 + 2), 2^(1/5), Inf)
# see Solution in the Derivation section;

# Derivation:
integrate(\(x) x^8 / (x^5 - 1) / (x^10 - 2*x^5 + 2), 2^(1/5), Inf)
integrate(\(x) 1/2 * x^3 * Re(2/(x^5 - 1) +
	- (1+1i)/(x^5 - (1+1i)) - (1-1i)/(x^5 - (1-1i))), 2^(1/5), Inf)
#
cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
integrate(\(x) 1/5 /(1-x) - 2/5 * (cs[1]*x-1) / (x^2 - 2*cs[1]*x + 1) +
	- 2/5*(cs[2]*x-1) / (x^2 - 2*cs[2]*x + 1), 0, 2^(-1/5))$value +
integrate(\(x) - 1/2 * x^3 * Re(
	(1+1i)/(x^5 - (1+1i)) + (1-1i)/(x^5 - (1-1i))), 2^(1/5), Inf)$value;
# Solution:
x1 = 2^(-1/5) * (1+1i)^(1/5); x2 = 2^(-1/5) * (1-1i)^(1/5);
- log(1 - 2^(-1/5))/5 +
	- 1/5 * sum(cs*log(2^(-2/5) - 2^(4/5)*cs + 1) +
	- 2*sn * (atan((2^(-1/5) - cs)/sn) + pi/2*c(1,-3)/5)) +
	+ 1/2 * Re(log(1 - x1) * (1+1i)^(4/5)/5 +
		+ (1+1i)^(4/5)/5 * sum(cs*log(x1^2 - 2*cs*x1 + 1) +
		- 2*sn * (atan((x1 - cs)/sn) + atan(cs/sn))) +
		+ log(1 - x2) * (1-1i)^(4/5)/5 +
		+ (1-1i)^(4/5)/5 * sum(cs*log(x2^2 - 2*cs*x2 + 1) +
		- 2*sn * (atan((x2 - cs)/sn) + atan(cs/sn))));

# Sub-integral 1:
# Note: atan(cs/sn) does NOT transform as nicely;
# atan(cs/sn) = pi/2 * c(1,-1) - atan(sn/cs)
#  = pi * c(1,-3)/10;
integrate(\(x) 1/(1 - x^5), 0, 2^(-1/5))
- log(1 - 2^(-1/5))/5 +
	- 1/5 * sum(cs*log(2^(-2/5) - 2^(4/5)*cs + 1) +
	- 2*sn * (atan((2^(-1/5) - cs)/sn) + atan(cs/sn)));
- log(1 - 2^(-1/5))/5 +
	- 1/5 * sum(cs*log(2^(-2/5) - 2^(4/5)*cs + 1) +
	- 2*sn * (atan((2^(-1/5) - cs)/sn) + pi/2*c(1,-3)/5));

#
integrate(\(x) x^3 * Re((1+1i)/(x^5 - (1+1i))), 2^(1/5), Inf)
# (1+1i)^(4/5) * x^3 / (x^5 - 1) on: [2^(1/5), Inf] / (1+1i)^(1/5);
# (1+1i)^(4/5) / (1 - x^5) on: [0, 2^(-1/5)] * (1+1i)^(1/5);
Re(pracma::line_integral(
	\(x) (1+1i)^(4/5)/(1 - x^5), c(0, 2^(-1/5)) * (1+1i)^(1/5)) )
x = 2^(-1/5) * (1+1i)^(1/5);
cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
- Re(log(1 - x) * (1+1i)^(4/5)/5 +
	+ (1+1i)^(4/5)/5 * sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))) );


# Fraction decomposition:
x = 3^(1/7)
x^8 / (x^5 - 1) / (x^10 - 2*x^5 + 2)
x^8 / (x^5 - 1) * (1/(x^5 - (1+1i)) - 1/(x^5 - (1-1i))) / 2i
1/2 * x^3 * (2/(x^5 - 1) - (1+1i)/(x^5 - (1+1i)) - (1-1i)/(x^5 - (1-1i)))


cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
1/(x^5 - 1)
(1/(x-1) + sum((2*cs*x-2)/(x^2 - 2*cs*x + 1))) / 5
#
x^3/(x^5 - 1)
(x^2 + x + 1 + 1/(x-1) + 2 * x^3 * sum((cs*x - 1)/(x^2 - 2*cs*x + 1))) / 5;
(x^2 + x + 1 + 1/(x-1) + 2 * sum(cs*x^2 - x + 2*cs^2*x - 3*cs + 4*cs^3 +
	+ (x - 8*cs^2*x + 8*cs^4*x + 3*cs - 4*cs^3)/(x^2 - 2*cs*x + 1))) / 5;
(1/(x-1) + 2 * sum(
	+ (cs*x + 3*cs - 4*cs^3)/(x^2 - 2*cs*x + 1))) / 5;

# source("Polynomials.Helper.R")

p1 = as.pm("cs*x^4 - x^3")
p2 = as.pm("x^2 - 2*cs*x + 1")
tmp = div.pm(p1, p2, by = "x")
print.pm(tmp$Rez, lead="x")
print.pm(tmp$Rem, lead="x")


###################
### Generalization:

### I( x^n * (1 + x^n)^(1 - 1/n) / (x^(2*n) + 1) )
n = 7; # ODD Integer
integrate(\(x) x^n * (1 + x^n)^(1 - 1/n) / (x^(2*n) + 1), 0, 1)
integrate(\(x) 1/2 * x^(n-2) * Re(2/(x^n - 1) +
	- (1+1i)/(x^n - (1+1i)) - (1-1i)/(x^n - (1-1i))), 2^(1/n), Inf)
# Solution:
id = seq(2, n-1, by=2); cs = cos(id*pi/n); sn = sin(id*pi/n);
x = 2^(-1/n); x1 = x * (1+1i)^(1/n); x2 = x * (1-1i)^(1/n);
- log(1 - x)/n +
	- 1/n * sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))) +
	+ 1/2 * Re(log(1 - x1) * (1+1i)^(1-1/n)/n +
		+ (1+1i)^(1-1/n)/n * sum(cs*log(x1^2 - 2*cs*x1 + 1) +
		- 2*sn * (atan((x1 - cs)/sn) + atan(cs/sn))) +
		+ log(1 - x2) * (1-1i)^(1-1/n)/n +
		+ (1-1i)^(1-1/n)/n * sum(cs*log(x2^2 - 2*cs*x2 + 1) +
		- 2*sn * (atan((x2 - cs)/sn) + atan(cs/sn))));


### Gen: Arbitrary Interval
lim = 4/5;
n = 7; # ODD Integer
ylim = (1 - 1/(lim^n + 1))^(-1/n);
integrate(\(x) x^n * (1 + x^n)^(1 - 1/n) / (x^(2*n) + 1), 0, lim)
integrate(\(x) 1/2 * x^(n-2) * Re(2/(x^n - 1) +
	- (1+1i)/(x^n - (1+1i)) - (1-1i)/(x^n - (1-1i))), ylim, Inf)
# Solution:
id = seq(2, n-1, by=2); cs = cos(id*pi/n); sn = sin(id*pi/n);
x = 1 / ylim; x1 = x * (1+1i)^(1/n); x2 = x * (1-1i)^(1/n);
- log(1 - x)/n +
	- 1/n * sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))) +
	+ 1/2 * Re(log(1 - x1) * (1+1i)^(1-1/n)/n +
		+ (1+1i)^(1-1/n)/n * sum(cs*log(x1^2 - 2*cs*x1 + 1) +
		- 2*sn * (atan((x1 - cs)/sn) + atan(cs/sn))) +
		+ log(1 - x2) * (1-1i)^(1-1/n)/n +
		+ (1-1i)^(1-1/n)/n * sum(cs*log(x2^2 - 2*cs*x2 + 1) +
		- 2*sn * (atan((x2 - cs)/sn) + atan(cs/sn))));


####################
####################

### Pow: 3*n

### Diff-Type

### I( x^3 * (1 - x^3)^(2/3) / (x^9 + 1) )
integrate(\(x) x^3 * (1 - x^3)^(2/3) / (x^9 + 1), 0, 1)
beta(1/3, 2/3) / 9 * (- 2^(2/3) + 2*cos(pi/9));
# see below for technique;

### I( x^5 * (1 - x^5)^(4/5) / (x^15 + 1) )
integrate(\(x) x^5 * (1 - x^5)^(4/5) / (x^15 + 1), 0, 1)
beta(1/5, 4/5) / 15 * (- 2^(4/5) + 2*cos(pi/15));

### Gen: I( x^n * (1 - x^n)^(1-1/n) / (x^(3*n) + 1) )
n = sqrt(5)
integrate(\(x) x^n * (1 - x^n)^(1-1/n) / (x^(3*n) + 1), 0, 1)
beta(1/n, 1-1/n) / (3*n) * (- 2^(1-1/n) + 2*cos(pi/(3*n)));


### Sum-Type

### I( x^3 * (1 + x^3)^(2/3) / (x^9 + 1) )
integrate(\(x) x^3 * (1 + x^3)^(2/3) / (x^9 + 1), 0, 1)
# Solution:
csn = cos(2*pi/3); snn = sin(2*pi/3);
m  = cos(2*pi/3) + c(1i, -1i)*sin(2*pi/3); x = 2^(-1/3) * (1 - m)^(1/3);
1i/(6*sn) * (diff((1 - m)^(-1/3) * log(1-x)) +
	- (1 - m[1])^(-1/3) * sum(csn*log(x[1]^2 - 2*csn*x[1] + 1) +
		- 2*snn * (atan((x[1] - csn)/snn) + atan(csn/snn))) +
	+ (1 - m[2])^(-1/3) * sum(csn*log(x[2]^2 - 2*csn*x[2] + 1) +
		- 2*snn * (atan((x[2] - csn)/snn) + atan(csn/snn))) );


### I( x^5 * (1 + x^5)^(4/5) / (x^15 + 1) )
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^15 + 1), 0, 1)
# Solution:
id = c(2,4); cs5 = cos(id*pi/5); sn5 = sin(id*pi/5);
m  = cos(2*pi/3) + c(1i, -1i)*sin(2*pi/3); x = 2^(-1/5) * (1 - m)^(1/5);
1i/(10*sn) * (diff((1 - m)^(-1/5) * log(1-x)) +
	- (1 - m[1])^(-1/5) * sum(cs5*log(x[1]^2 - 2*cs5*x[1] + 1) +
		- 2*sn5 * (atan((x[1] - cs5)/sn5) + atan(cs5/sn5))) +
	+ (1 - m[2])^(-1/5) * sum(cs5*log(x[2]^2 - 2*cs5*x[2] + 1) +
		- 2*sn5 * (atan((x[2] - cs5)/sn5) + atan(cs5/sn5))) );


### Gen: I( x^n * (1 + x^n)^(1-1/n) / (x^(3*n) + 1) )
n = 7; # ODD Integer;
integrate(\(x) x^n * (1 + x^n)^(1-1/n) / (x^(3*n) + 1), 0, 1)
# Solution:
id = seq(2, n-1, by=2); csn = cos(id*pi/n); snn = sin(id*pi/n);
m  = cos(2*pi/3) + c(1i, -1i)*sin(2*pi/3); x = 2^(-1/n) * (1 - m)^(1/n);
1i/(2*n*sn) * (diff((1 - m)^(-1/n) * log(1-x)) +
	- (1 - m[1])^(-1/n) * sum(csn*log(x[1]^2 - 2*csn*x[1] + 1) +
		- 2*snn * (atan((x[1] - csn)/snn) + atan(csn/snn))) +
	+ (1 - m[2])^(-1/n) * sum(csn*log(x[2]^2 - 2*csn*x[2] + 1) +
		- 2*snn * (atan((x[2] - csn)/snn) + atan(csn/snn))) );


# Derivation: Variant (1 - x^5)
integrate(\(x) x^5 * (1 - x^5)^(4/5) / (x^15 + 1), 0, 1)
integrate(\(x) 1/5 * (x * (1 - x)^4)^(1/5) / (x^3 + 1), 0, 1)
integrate(\(x) x^5 / (1+x^5)^3 / ((x^5/(1+x^5))^3 + 1), 0, Inf)
integrate(\(x) x^5 / ((x^5 + 1)^3 + x^15), 0, Inf)
integrate(\(x) x^8 / ((x^5 + 1)^3 + 1), 0, Inf)
#
cs = cos(2*pi/3); sn = sin(2*pi/3);
m  = cs + c(1i, -1i)*sn;
integrate(\(x) -2/3 * x^3 / (x^5 + 2) +
	+ 1/3 * x^3 * Re(1/(x^5 + 1 + m[1]) +
	+ 1/(x^5 + 1 + m[2])), 0, Inf)
integrate(\(x) - 2^(4/5)/3 * x^3 / (x^5 + 1), 0, Inf)$value +
integrate(\(x) + Re((1+m[1])^(-1/5))/3 * x^3 / (x^5 + 1), 0, Inf)$value +
integrate(\(x) + Re((1+m[2])^(-1/5))/3 * x^3 / (x^5 + 1), 0, Inf)$value;
beta(1/5,4/5) / 15 * (- 2^(4/5) + (1+m[1])^(-1/5) + (1+m[2])^(-1/5));
beta(1/5,4/5) / 15 * (- 2^(4/5) + 2*cos(pi/15));


# Derivation: Variant (1 + x^5)
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^15 + 1), 0, 1)
integrate(\(x) 1/5 * (x * (1 + x)^4)^(1/5) / (x^3 + 1), 0, 1)
integrate(\(x) x^5 / (1-x^5)^3 / ((x^5/(1-x^5))^3 + 1), 0, 2^(-1/5))
integrate(\(x) x^5 / ((1 - x^5)^3 + x^15), 0, 2^(-1/5))
integrate(\(x) x^8 / ((x^5 - 1)^3 + 1), 2^(1/5), Inf)
#
cs = cos(2*pi/3); sn = sin(2*pi/3);
id = c(2,4); cs5 = cos(id*pi/5); sn5 = sin(id*pi/5);
m  = cs + c(1i, -1i)*sn; lim = 2^(-1/5); lz = lim * (1 - m)^(1/5);
integrate(\(x) -1/(2*sn) * x^3 *
	Im(1/(x^5 - 1 + m[1]) - 1/(x^5 - 1 + m[2])), 2^(1/5), Inf)
pracma::line_integral(\(x) - 1i/(2*sn) * (1 - m[1])^(-1/5) /(x^5 - 1), c(0, lz[1])) +
pracma::line_integral(\(x) 1i/(2*sn) * (1 - m[2])^(-1/5) /(x^5 - 1), c(0, lz[2]))
#
x = lz;
1i/(10*sn) * (diff((1 - m)^(-1/5) * log(1-x)) +
	- (1 - m[1])^(-1/5) * sum(cs5*log(x[1]^2 - 2*cs5*x[1] + 1) +
		- 2*sn5 * (atan((x[1] - cs5)/sn5) + atan(cs5/sn5))) +
	+ (1 - m[2])^(-1/5) * sum(cs5*log(x[2]^2 - 2*cs5*x[2] + 1) +
		- 2*sn5 * (atan((x[2] - cs5)/sn5) + atan(cs5/sn5))) );
# Variant:
1i/(10*sn) * (diff((1 - m)^(-1/5) * log(1-x)) +
	+ 1i * 2^(9/5) * sin(pi/30) / sin(pi/3)^(1/5) * sum(sn5 * atan(cs5/sn5)) +
	- (1 - m[1])^(-1/5) * sum(cs5*log(x[1]^2 - 2*cs5*x[1] + 1) +
		- 2*sn5 * atan((x[1] - cs5)/sn5)) +
	+ (1 - m[2])^(-1/5) * sum(cs5*log(x[2]^2 - 2*cs5*x[2] + 1) +
		- 2*sn5 * atan((x[2] - cs5)/sn5)) );


# Fraction Decomposition:
x = 1/5^(1/7)
cs = cos(2*pi/3); sn = sin(2*pi/3);
m  = cs + c(-1i, 1i)*sn;
1 / ((x^5 + 1)^3 + 1)
1/(x^5 + 2) * (1/(x^5 + 1 + m[1]) - 1/(x^5 + 1 + m[2])) / (2i*sn)
diff((1/(x^5 + 2) - 1/(x^5 + 1 + m)) / (2i*sn*(1-m)))
sum(-1/3 * m * (1/(x^5 + 2) - 1/(x^5 + 1 + m)))

# Variant: (x^5 - 1)
1 / ((x^5 - 1)^3 + 1)
1/x^5 * (1/(x^5 - 1 + m[1]) - 1/(x^5 - 1 + m[2])) / (2i*sn)
diff((1/x^5 - 1/(x^5 - 1 + m)) / (2i*sn*(1-m)))
sum(-1/3 * m * (1/x^5 - 1/(x^5 - 1 + m))) # not used;

#
x^8 / ((x^5 + 1)^3 + 1)
1/3 * x^8 / (x^5 + 2) +
	+ 1/3 * x^8 * Re(m[1]/(x^5 + 1 + m[1]) + m[2]/(x^5 + 1 + m[2]));
-2/3 * x^3 / (x^5 + 2) - 1/3 * x^3 *
	Re(m[1]*(1+m[1])/(x^5 + 1 + m[1]) +
	+ m[2]*(1+m[2])/(x^5 + 1 + m[2]));


#################
#################

### Simple

### Pow = 5

### I( x^5 * (1 + x^5)^(4/5) / (x^5 + 1) )
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^5 + 1), 0, 1)
# Solution:
x = 2^(-1/5); id = c(2,4); cs = cos(id*pi/5); sn = sin(id*pi/5);
1/5 * x / (1 - x^5) + 1/25*(log(1-x) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

# Simplified formula for [0,1] is available in file:
# Integrals.Fractions.Unity.Radicals.R;


### Arbitrary Interval:
lim = 4/5
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^5 + 1), 0, lim)
# Solution:
x = lim/(lim^5 + 1)^(1/5);
id = c(2,4); cs = cos(id*pi/5); sn = sin(id*pi/5);
1/5 * x / (1 - x^5) + 1/25*(log(1-x) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

### I( 1 / (x^5 + 1)^(1/5) )
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) 1 / (x^5 + 1)^(1/5), 0, lim)
x  = lim / (lim^5 + 1)^(1/5);
id = c(2,4); cs = cos(id*pi/5); sn = sin(id*pi/5);
lim*(1 + lim^5)^(4/5) - x / (1 - x^5) - 1/5*(log(1-x) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));


### Pow 6

### I( x^6 * (1 + x^6)^(5/6) / (x^6 + 1) )
integrate(\(x) x^6 * (1 + x^6)^(5/6) / (x^6 + 1), 0, 1)
# Solution:
x = 2^(-1/6); id = c(2,4); cs = cos(id*pi/6); sn = sin(id*pi/6);
1/6 * x / (1 - x^6) + 1/36*(log((1-x)/(1+x)) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

### I( x^6 * (1 + x^6)^(5/6) / (x^6 + 1) )
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) x^6 * (1 + x^6)^(5/6) / (x^6 + 1), 0, lim)
# Solution:
x  = lim / (lim^6 + 1)^(1/6);
id = c(2,4); cs = cos(id*pi/6); sn = sin(id*pi/6);
1/6 * x / (1 - x^6) + 1/36*(log((1-x)/(1+x)) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

### I( 1 / (x^6 + 1)^(1/6) )
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) 1 / (x^6 + 1)^(1/6), 0, lim)
x  = lim / (lim^6 + 1)^(1/6);
id = c(2,4); cs = cos(id*pi/6); sn = sin(id*pi/6);
- 1/6*(log((1-x)/(1+x)) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));


### Generalization:

### I( 1 / (x^n + 1)^(1/n) )
n = 7; # ODD Integer
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) 1 / (x^n + 1)^(1/n), 0, lim)
x  = lim / (lim^n + 1)^(1/n);
id = seq(2, n-1, by=2); cs = cos(id*pi/n); sn = sin(id*pi/n);
- 1/n*(log(1-x) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

### I( 1 / (x^n + 1)^(1/n) )
n = 8; # EVEN Integer
lim = 6/7; # Arbitrary Interval: can be > 1;
integrate(\(x) 1 / (x^n + 1)^(1/n), 0, lim)
x  = lim / (lim^n + 1)^(1/n);
id = seq(2, n-2, by=2); cs = cos(id*pi/n); sn = sin(id*pi/n);
- 1/n*(log((1-x)/(1+x)) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

# Note:
lim*(1 + lim^n)^(1-1/n) - x / (1 - x^n) # == 0!

# Derivation: Variant (1 + x^5)
integrate(\(x) x^5 * (1 + x^5)^(4/5) / (x^5 + 1), 0, 1)
integrate(\(x) 1/5 * (x * (1 + x)^4)^(1/5) / (x + 1), 0, 1)
integrate(\(x) x^5 / (1-x^5)^3 / ((x^5/(1-x^5)) + 1), 0, 2^(-1/5))
integrate(\(x) x^5 / (1 - x^5)^2 / ((1 - x^5) + x^5), 0, 2^(-1/5))
integrate(\(x) 1 / (1 - x^5)^2 - 1 / (1 - x^5), 0, 2^(-1/5))
x = 2^(-1/5); id = c(2,4); cs = cos(id*pi/5); sn = sin(id*pi/5);
1/5 * x / (1 - x^5) + 1/25*(log(1-x) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));


# Derivation: Variant (1 + x^6)
integrate(\(x) x^6 * (1 + x^6)^(5/6) / (x^6 + 1), 0, 1)
integrate(\(x) x^6 / (1 - x^6)^2 / ((1 - x^6) + x^6), 0, 2^(-1/6))
integrate(\(x) 1 / (1 - x^6)^2 - 1 / (1 - x^6), 0, 2^(-1/6))
x = 2^(-1/6); id = c(2,4); cs = cos(id*pi/6); sn = sin(id*pi/6);
1/6 * x / (1 - x^6) + 1/36*(log((1-x)/(1+x)) +
	+ sum(cs*log(x^2 - 2*cs*x + 1) +
	- 2*sn * (atan((x - cs)/sn) + atan(cs/sn))));

# Reduction: Pow = 2
x = 2^(-1/5);
integrate(\(x) 1 / (1 - x^5), 0, x)
x / (1 - x^5) - integrate(\(x) 5*x^5 / (1 - x^5)^2, 0, x)$value
x / (1 - x^5) - integrate(\(x) 5/(1 - x^5)^2 - 5/ (1 - x^5), 0, x)$value
# =>
integrate(\(x) 1/(1 - x^5)^2, 0, x)
1/5 * x / (1 - x^5) + integrate(\(x) 4/5 / (1 - x^5), 0, x)$value

# Reduction: Pow = 2
x = 2^(-1/5);
integrate(\(x) 1/(1 + x^5)^2, 0, x)
1/5 * x / (1 + x^5) + integrate(\(x) 4/5 / (1 + x^5), 0, x)$value

#
x = 2^(-1/5);
integrate(\(x) (1 + x^5)^(4/5), 0, x)
x*(1 + x^5)^(4/5) - integrate(\(x) 4*x^5 / (1 + x^5)^(1/5), 0, x)$value
# =>
integrate(\(x) (1 + x^5)^(4/5), 0, x)
1/5 * x*(1 + x^5)^(4/5) + integrate(\(x) 4/5 / (1 + x^5)^(1/5), 0, x)$value


#
n = 6; lim = 4/5;
x = lim / (lim^n + 1)^(1/n);
lim*(1 + lim^n)^(1-1/n) - x / (1 - x^n) # == 0!


#####################

### Simple: Power Set [p, n-p]

### I( x^6 * (1 - x^5)^(3/5) / (x^5 + 1) )
integrate(\(x) x^6 * (1 - x^5)^(3/5) / (x^5 + 1), 0, 1)
(2 - 2/5 - 2^(3/5)) * beta(2/5, 3/5) / 5

### I( x^7 * (1 - x^5)^(2/5) / (x^5 + 1) )
integrate(\(x) x^7 * (1 - x^5)^(2/5) / (x^5 + 1), 0, 1)
(2 - 3/5 - 2^(2/5)) * beta(2/5, 3/5) / 5

### I( x^8 * (1 - x^5)^(1/5) / (x^5 + 1) )
integrate(\(x) x^8 * (1 - x^5)^(1/5) / (x^5 + 1), 0, 1)
(2 - 4/5 - 2^(1/5)) * beta(1/5, 4/5) / 5

### Simplified

### I( x * (1 - x^5)^(3/5) / (x^5 + 1) )
integrate(\(x) x * (1 - x^5)^(3/5) / (x^5 + 1), 0, 1)
(2^(3/5) - 1) * beta(2/5, 3/5) / 5

### I( x^2 * (1 - x^5)^(2/5) / (x^5 + 1) )
integrate(\(x) x^2 * (1 - x^5)^(2/5) / (x^5 + 1), 0, 1)
(2^(2/5) - 1) * beta(2/5, 3/5) / 5

### I( x^3 * (1 - x^5)^(1/5) / (x^5 + 1) )
integrate(\(x) x^3 * (1 - x^5)^(1/5) / (x^5 + 1), 0, 1)
(2^(1/5) - 1) * beta(1/5, 4/5) / 5


### I( x * (1 - x^5)^(3/5) )
integrate(\(x) x * (1 - x^5)^(3/5), 0, 1)
integrate(\(x) 1/5 * x^(2/5-1) * (1 - x)^(3/5), 0, 1)
beta(2/5, 8/5) / 5

### I( x^2 * (1 - x^5)^(2/5) )
integrate(\(x) x^2 * (1 - x^5)^(2/5), 0, 1)
integrate(\(x) 1/5 * x^(3/5-1) * (1 - x)^(2/5), 0, 1)
beta(3/5, 7/5) / 5

### I( x^3 * (1 - x^5)^(1/5) )
integrate(\(x) x^3 * (1 - x^5)^(1/5), 0, 1)
integrate(\(x) 1/5 * x^(4/5-1) * (1 - x)^(1/5), 0, 1)
beta(4/5, 6/5) / 5


# Derivation:
integrate(\(x) x^6 * (1 - x^5)^(3/5) / (x^5 + 1), 0, 1)
integrate(\(x) 1/5 * (x^2 * (1 - x)^3)^(1/5) / (x + 1), 0, 1)
integrate(\(x) x^6 / (1+x^5)^3 / ((x^5/(1+x^5)) + 1), 0, Inf)
integrate(\(x) x^6 / (1+x^5)^2 / (2*x^5 + 1), 0, Inf)
integrate(\(x) x^7 / (x^5 + 1)^2 / (x^5 + 2), 0, Inf)
integrate(\(x) x^2 * (2/(x^5 + 1) - 1/(x^5 + 1)^2 - 2/(x^5 + 2)), 0, Inf)
# Simplification on [0, Inf];
integrate(\(x) x^2 * ((2 - 2^(3/5))/(x^5 + 1) - 1/(x^5 + 1)^2), 0, Inf)
integrate(\(x) (2 - 2/5 - 2^(3/5)) * x^2 / (x^5 + 1), 0, Inf)
(2 - 2/5 - 2^(3/5)) * beta(2/5, 3/5) / 5


# Derivation:
integrate(\(x) x^7 * (1 - x^5)^(2/5) / (x^5 + 1), 0, 1)
integrate(\(x) 1/5 * (x^3 * (1 - x)^2)^(1/5) / (x + 1), 0, 1)
integrate(\(x) x^7 / (1+x^5)^3 / ((x^5/(1+x^5)) + 1), 0, Inf)
integrate(\(x) x^7 / (1+x^5)^2 / (2*x^5 + 1), 0, Inf)
integrate(\(x) x^6 / (x^5 + 1)^2 / (x^5 + 2), 0, Inf)
integrate(\(x) x * (2/(x^5 + 1) - 1/(x^5 + 1)^2 - 2/(x^5 + 2)), 0, Inf)
# Simplification on [0, Inf];
integrate(\(x) x * ((2 - 2^(2/5))/(x^5 + 1) - 1/(x^5 + 1)^2), 0, Inf)
integrate(\(x) (2 - 3/5 - 2^(2/5)) * x / (x^5 + 1), 0, Inf)
(2 - 3/5 - 2^(2/5)) * beta(2/5, 3/5) / 5


# Derivation:
integrate(\(x) x^8 * (1 - x^5)^(1/5) / (x^5 + 1), 0, 1)
integrate(\(x) 1/5 * (x^4 * (1 - x))^(1/5) / (x + 1), 0, 1)
integrate(\(x) x^8 / (1+x^5)^3 / ((x^5/(1+x^5)) + 1), 0, Inf)
integrate(\(x) x^8 / (1+x^5)^2 / (2*x^5 + 1), 0, Inf)
integrate(\(x) x^5 / (x^5 + 1)^2 / (x^5 + 2), 0, Inf)
integrate(\(x) (2/(x^5 + 1) - 1/(x^5 + 1)^2 - 2/(x^5 + 2)), 0, Inf)
# Simplification on [0, Inf];
integrate(\(x) ((2 - 2^(1/5))/(x^5 + 1) - 1/(x^5 + 1)^2), 0, Inf)
integrate(\(x) (2 - 4/5 - 2^(1/5)) / (x^5 + 1), 0, Inf)
(2 - 4/5 - 2^(1/5)) * beta(1/5, 4/5) / 5


# Reduction: Pow = 2
x = 2^(-1/5);
integrate(\(x) x^2 / (x^5 + 1)^2, 0, x)
1/5 * x^3 / (x^5 + 1) + integrate(\(x) 2/5 * x^2 / (x^5 + 1), 0, x)$value

#
x = 2^(-1/5);
integrate(\(x) x / (x^5 + 1)^2, 0, x)
1/5 * x^2 / (x^5 + 1) + integrate(\(x) 3/5 * x / (x^5 + 1), 0, x)$value

