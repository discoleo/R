########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### Polynomial Radicals
###
### draft v.0.2f


### Types:

# I( x^p / (x^n + 1)^k ) on [0, Inf]
# I( x^p / (x^n + 1)^k ) on [0, 1]
# I( x^[p]  * (x^n - 1)^k1 / (x^n + 1)^k2 ) on [0, 1]
#   where p,n,k = integers (and not yet independent);


### History

# - moved Integrals with Radicals
#   to this file from file:
#   Integrals.Fractions.Unity.Definite.R;
# - moved other variants to new file:
#   Integrals.Fractions.Unity.Radicals.Other.R;
# - moved "Mixed Radicals" to new file:
#   Integrals.Fractions.Unity.Radicals.Mixed.R;


######################

### Helper Functions

Euler = 0.57721566490153286060651209008240243079;

### I on [0, 1]
# I( x^p / (x^n + 1))
int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}
int.FrDU01 = function(n, p=0) {
	if( p!= 0) {
		# r = int.FrDU01(1, p) - int.FrDU01(n, p);
		r = digamma((p+1)) - digamma((p+1)/n) - log(n);
		return(r);
	}
	digamma(1/n) + Euler + log(n);
}

### I on [0, Inf]
# I( Polynomial(x as x^p) / (x^n + 1)^k )
int.FrUInf = function(n, p=0, pow=1, coeff=1) {
	k = 1/pow; # actually (x^n + 1)^(1/k);
	tmp = sapply(p, function(p) {
		gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
	});
	tmp = sum(coeff * tmp);
	return(tmp);
}

##########################
##########################

################
### Radicals ###
################

### [0, Inf]
p = 1 - sqrt(2)
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n
beta((p+1)/n, 1/k - (p+1)/n) / n


### "Divergent Radicals"

###
integrate(\(x) (x^3 + 1)^(1/3) - x, 0, Inf)
gamma(-1/3 - 1/3) / gamma(-1/3) * gamma(1/3) / 3


### I( (x^n + 1)^(1/n) ) on [0, Inf]
n = sqrt(5)
integrate(\(x) (x^n + 1)^(1/n) - x, 0, Inf)
gamma(-2/n) / gamma(-1/n) * gamma(1/n) / n


### I( (x^n + 1)^(2/n) ) on [0, Inf]
n = sqrt(11); # n > 3;
integrate(\(x) (x^n + 1)^(2/n) - x^2, 0, Inf)
gamma(-3/n) / gamma(-2/n) * gamma(1/n) / n


### I( 1 / (x^n + 1)^(1/n) ) on [0, Inf]
# OK for any n > 0;
n = 6 # n = sqrt(3);
integrate(\(x) 1 / (x^n + 1)^(1/n) - 1/(x+1), 0, Inf)
- (digamma(1/n) + Euler) / n


### Special Cases:

### n = 3
# Convergence: more complicated;
integrate(\(x) (x^3 + 1)^(2/3) - x^2 - 2/3 / (x+1), 0, Inf)
- (digamma(-2/3) + Euler - 5/2) * 2/9
pi*sqrt(3) / 27 + log(3) / 3 + 2/9;

# library(Rmpfr)
integrate(\(x) {
	x = mpfr(x, 240); f23 = mpfr(2, 240) / 3;
	y = (x^3 + 1)^f23 - x^2 - f23 / (x+1);
	as.numeric(y); }, 0, Inf, subdivisions=1024)
- (digamma(-2/3) + Euler - 1) * 2/9 + 1/3;

### n = 4
integrate(\(x) {
	x = mpfr(x, 240);
	y = (x^4 + 1)^(3/4) - x^3 - 3/4 / (x+1);
	as.numeric(y); }, 0, Inf)
- (digamma(-3/4) + Euler - 1) * 3/16 + 1/4;

### n = 5
n = 5; # n = sqrt(7); # OK
integrate(\(x, n) {
	x = mpfr(x, 240); ni = 1 - mpfr(1, 240)/n;
	y = (x^n + 1)^ni - x^(n - 1) - ni / (x+1);
	as.numeric(y); }, 0, Inf, n=n)
- (digamma(1/n - 1) + Euler - 1) * (n-1)/n^2 + 1/n;


### Mixed Fractions

### I( 1 / (x^2+1)^(1/2) * 1/(x+1) )
integrate(\(x) 1 / (x^2+1)^(1/2) * 1/(x+1), 0, Inf)
log(3 + 2*sqrt(2)) * sqrt(2) / 2;


### I( (x^2+1)^(1/2) / (x+1) )
integrate(\(x) (x^2+1)^(1/2) / (x+1) - x/(x+1), 0, Inf)
2*sqrt(2) * log(sqrt(2)+1) - log(2) - 1;

# Derivation:
integrate(\(x) (1 - sin(x)) / (sin(x) + cos(x)) / cos(x)^2, 0, pi/2)
integrate(\(x) (1 - cos(x)) / (sin(x) + cos(x)) / sin(x)^2, 0, pi/2)
integrate(\(x) (cos(x) - sin(x) - cos(x)^2 + sin(x)*cos(x)) / cos(2*x) / sin(x)^2, 0, pi/2)
integrate(\(x) (cos(x) - sin(x) + sin(x)*cos(x) - 1) / cos(2*x) / sin(x)^2 +
	+ 1/tan(2*x), 0, pi/2);
#
integrate(\(x) (cos(x) - sin(x) + sin(x)*cos(x)) / cos(2*x) / sin(x)^2 +
	+ 1/tan(2*x) - 1/sin(x)^2 - 1/2/sin(x) + 1/2/cos(x) - 1/cos(2*x), 0, pi/2);
integrate(\(x) cos(x) / cos(2*x) / sin(x)^2 + cos(x) / cos(2*x) / sin(x) +
	- 1/cos(2*x)/sin(x) - 1/x^2 - (2/pi)^2 +
	- 1/cos(2*x), 0, pi/2);
#
integrate(\(x) sqrt(2)/2 / (pi/4-x) + 1/x - 1/cos(2*x)/sin(x), 0, pi/4)$value +
integrate(\(x) sqrt(2)/2 / (pi/4-x) + 1/x - 1/cos(2*x)/sin(x), pi/4, pi/2)$value +
	- log(pi/2) + log(3 + 2*sqrt(2)) * sqrt(2) / 2 - 1;
#
2*sqrt(2) * log(sqrt(2)+1) - log(2) - 1;


# I( cos(x) / cos(2*x) )
integrate(\(x) cos(x) / (1 - 2*sin(x)^2) + sqrt(2)/4 / (x-pi/4), pi/4, pi/2)
sqrt(2)/4 * log(pi/4) - log(2)/2/sqrt(2) +
	+ (digamma(5/8) - digamma(1/8) - digamma(7/8) + digamma(3/8)) / 8;
sqrt(2)/4 * log(pi/4) - log(2)/2/sqrt(2) + log(3 + 2*sqrt(2)) * sqrt(2) / 4;

# I( cos(x) / cos(2*x) / sin(x) )
integrate(\(x) cos(x) / cos(2*x) / sin(x) - 1/x - 1/2/(pi/4-x), 0, pi/4)
-1/2*log(pi/4) - log(pi/2);
#
integrate(\(x) cos(x) / cos(2*x) / sin(x) - 1/x - 1/2/(pi/4-x), pi/4, pi/2)
log(pi/4)/2;

# I( 1/cos(2*x) / sin(x) )
integrate(\(x) sqrt(2)/2 / (pi/4-x) + 1/x - 1/cos(2*x)/sin(x), 0, pi/4)
log(pi/4) + sqrt(2)/2 * log(pi/4) - log(2)/2 +
	+ (log((1+sqrt(2)/2)/(1-sqrt(2)/2))/2 + log(1/2) / sqrt(2)) +
	- (log(2)/2 + log((sqrt(2)-1)/(1+sqrt(2))) / sqrt(2));
(1 + sqrt(2)/2) * log(pi) - 3*(1 + sqrt(2)/2) * log(2) - log(2)/2 +
	+ log(2+sqrt(2)) - log((sqrt(2)-1)) * sqrt(2);
#
integrate(\(x) sqrt(2)/2 / (pi/4-x) + 1/x - 1/cos(2*x)/sin(x), pi/4, pi/2)
- sqrt(2)/2 * log(pi/4) + log(2) + log(2) / sqrt(2) +
	- log((1+sqrt(2)/2)/(1-sqrt(2)/2)) / 2;
#
up = 1/5
integrate(\(x) 2/(2*x^2-1) + 1/ (1-x^2), 0, up)
log((1+up)/(1-up))/2 + log((1-sqrt(2)*up)/(1+sqrt(2)*up)) / sqrt(2);
# Limits:
x = 1E-4;
log(1-cos(x))/2 - log(x);
- log(2)/2;
#
x = pi/4 - 1E-4;
log(sqrt(2)*cos(x) - 1) / sqrt(2) - log(pi/4-x) / sqrt(2);
0;


# alternative:
integrate(\(x) 1 / (2*sin(x)*cos(x) + 2*cos(x)^2 - 1) / cos(x)^2, 0, pi/4)


#################
#################

##############
### [0, 1] ###

### Note:
# - the fraction decomposition of x^p / (x^n - 1)
#   (for n, p = integers) is described in file:
#   Integrals.Fractions.Unity.R;

### I( 1 / (x^3 + 1)^(1/3) )
integrate(\(x) 1 / (x^3 + 1)^(1/3), 0, 1)
- (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan((2^(1/3) + 1/2) * 2/sqrt(3)) +
	+ 1/sqrt(3) * atan((1 + 1/2) * 2/sqrt(3));
# Simplification:
- (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan(1/sqrt(3) * (2^(1/3) - 1) / (2^(1/3) + 1));

### I( (x^3 + 1)^(1/3) )
integrate(\(x) (x^3 + 1)^(1/3), 0, 1)
2^(-2/3) + gamma(1/3)^2/gamma(2/3) / 12;


### I( x / (x^3 + 1)^(1/3) )
# - based on simplified formula: see below;
integrate(\(x) x / (x^3 + 1)^(1/3), 0, 1)
1/2^(1/3) - gamma(2/3)^2/gamma(4/3) / 6


### Pow = 2/3

### I( 1 / (x^3 + 1)^(2/3) )
# - based on simplified formula: see below;
integrate(\(x) 1 / (x^3 + 1)^(2/3), 0, 1)
gamma(1/3)^2/gamma(2/3) / 6

### I( (x^3 + 1)^(2/3) )
integrate(\(x) (x^3 + 1)^(2/3), 0, 1)
(2^(-1/3) - (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan(1/sqrt(3) * (2^(1/3)-1) / (2^(1/3)+1)) ) * 2/3;

### I( x / (x^3 + 1)^(2/3) )
integrate(\(x) x / (x^3 + 1)^(2/3), 0, 1)
- (digamma(2/3) + Euler)/3 +
	+ integrate(\(x) (x^2 - 1) / (x^3 - 1), 1, 2^(1/3))$value;
- (digamma(2/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	+ 1/sqrt(3) * atan(1/sqrt(3) * (2^(1/3) - 1) / (2^(1/3) + 1));


###########
### Pow = 5

### I( 1 / (x^5 + 1)^(1/5) )
integrate(\(x) 1 / (x^5 + 1)^(1/5), 0, 1)
- (digamma(1/5) + Euler)/5 +
	+ integrate(\(x) x^3 * (x-1) / (x^5 - 1), 1, 2^(1/5))$value;
x = 2^(1/5); cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
- (digamma(1/5) + Euler)/5 + sum(
	+ 1/5 * (1 - cs) * log(x^2 - 2*cs*x + 1) +
	- 1/5 * (1 - cs) * log(2 - 2*cs) +
	- 2/5 * sn * atan((x - cs) / sn) +
	+ 2/5 * sn * atan((1 - cs) / sn) );


### I( x   / (x^5 + 1)^(1/5) ) &
### I( x^3 / (x^5 + 1)^(1/5) )
integrate(\(x) x / (x^5 + 1)^(1/5), 0, 1)
2^(4/5) - gamma(4/5) * gamma(2/5) / gamma(1/5) +
	- integrate(\(x) 3 * x^3 / (x^5 + 1)^(1/5), 0, 1)$value;
# TODO: ???

### I( x^2 / (x^5 + 1)^(1/5) )
integrate(\(x) x^2 / (x^5 + 1)^(1/5), 0, 1)
1/2^(6/5) - gamma(3/5)^2/gamma(1/5) / 4


### I( 1 / (x^5 + 1)^(2/5) )
integrate(\(x) 1 / (x^5 + 1)^(2/5), 0, 1)
gamma(1/5)^2/gamma(2/5) / 10

### I( x / (x^5 + 1)^(2/5) )
integrate(\(x) x / (x^5 + 1)^(2/5), 0, 1)
- (digamma(2/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - x^2) / (x^5 - 1), 1, 2^(1/5))$value;
# Solution:
x = 2^(1/5)
cs = cos(c(2,4)*pi/5); cs2 = cos(c(2,4)*3*pi/5);
sn = sin(c(2,4)*pi/5); sn2 = sin(c(2,4)*3*pi/5);
cs4 = cos(c(2,4)*pi); sn4 = sin(c(2,4)*pi);
- (digamma(2/5) + Euler)/5 +
	+ sum((cs4-cs2) * log((x^2 - 2*cs*x + 1) / (2-2*cs)) +
	- 2*(sn4-sn2) * (atan((x - cs)/sn) - atan((1-cs)/sn)) ) / 5;

# Alternative formula with Arbitrary Interval, see file:
# Integrals.Fractions.Unity.Radicals.Mixed.R
n = 5; p = 1; # Integer in [0, n-1]
lim = 4/5; # Arbitrary Interval: can be > 1;
integrate(\(x) x^p / (x^n + 1)^((p+1)/n), 0, lim)
id = seq(2, n-1, by=2); x = lim / (lim^n + 1)^(1/n); 
cs = cos(id*pi/n); csp = cos(id*(p+1)*pi/n);
sn = sin(id*pi/n); snp = sin(id*(p+1)*pi/n);
- 1/n*(log(1-x) +
	+ sum(csp*log(x^2 - 2*cs*x + 1) +
	- 2*snp * (atan((x - cs)/sn) + atan(cs/sn))));


### I( x^2 / (x^5 + 1)^(2/5) ) &
### I( x^3 / (x^5 + 1)^(2/5) )
integrate(\(x) x^2 / (x^5 + 1)^(2/5), 0, 1)
2^(3/5) - gamma(4/5) * gamma(3/5) / gamma(2/5) +
	- integrate(\(x) 2 * x^3 / (x^5 + 1)^(2/5), 0, 1)$value;
# TODO: ???


### I( x^2 / (x^5 + 1)^(3/5) )
integrate(\(x) x^2 / (x^5 + 1)^(3/5), 0, 1)
- (digamma(3/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - x) / (x^5 - 1), 1, 2^(1/5))$value;
# Solution:
# Note: for Arbitrary Interval, see a few paragraphs below;
cs = cos(c(2,4)*pi/5); csp = cos(3*c(2,4)*pi/5);
sn = sin(c(2,4)*pi/5); snp = sin(3*c(2,4)*pi/5); x = 2^(-1/5);
-1/5 * log(1-x) - sum(csp*log(x^2 - 2*cs*x + 1) +
	- 2*snp*(atan((x - cs)/sn) + atan(cs/sn))) / 5;

### I( x^3 / (x^5 + 1)^(3/5) )
integrate(\(x) x^3 / (x^5 + 1)^(3/5), 0, 1)
1/2^(3/5) - gamma(4/5)^2/gamma(3/5) / 2


### I( x / (x^5 + 1)^(4/5) )
integrate(\(x) x / (x^5 + 1)^(4/5), 0, 1)
gamma(2/5)^2/gamma(4/5) / 10

### I( x^3 / (x^5 + 1)^(4/5) )
integrate(\(x) x^3 / (x^5 + 1)^(4/5), 0, 1)
- (digamma(4/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - 1) / (x^5 - 1), 1, 2^(1/5))$value;
# TODO: compute;


### Pow = 7

### I( 1 / (x^7 + 1)^(1/7) )
integrate(\(x) 1 / (x^7 + 1)^(1/7), 0, 1)
- (digamma(1/7) + Euler)/7 +
	+ integrate(\(x) x^5 * (x - 1) / (x^7 - 1), 1, 2^(1/7))$value
x = 2^(1/7); cs = cos(c(2,4,6)*pi/7); sn = sin(c(2,4,6)*pi/7);
- (digamma(1/7) + Euler)/7 + sum(
	+ 1/7 * (1 - cs) * log(x^2 - 2*cs*x + 1) +
	- 1/7 * (1 - cs) * log(2 - 2*cs) +
	- 2/7 * sn * atan((x - cs) / sn) +
	+ 2/7 * sn * atan((1 - cs) / sn) );


### Higher Powers

### I( 1 / (x^8 + 1)^(1/8) )
n = 8
integrate(\(x) 1 / (x^n + 1)^(1/n), 0, 1)
- (digamma(1/n) + Euler)/n +
	+ integrate(\(x) x^(n-2) * (x - 1) / (x^n - 1), 1, 2^(1/n))$value

### I( 1 / (x^9 + 1)^(1/9) )
n = 9
# Note: n = ODD for explicit formula!
integrate(\(x) 1 / (x^n + 1)^(1/n), 0, 1)
- (digamma(1/n) + Euler)/n +
	+ integrate(\(x) x^(n-2) * (x - 1) / (x^n - 1), 1, 2^(1/n))$value
id = seq(2, n-1, by = 2)
x = 2^(1/n); cs = cos(id*pi/n); sn = sin(id*pi/n);
pi/(2*n) / tan(pi/n) + pi/(2*n) / sin(pi/n) +
+ sum(
	+ 1/n * (1 - cs) * log(x^2 - 2*cs*x + 1) +
	- 2/n * sn * atan((x - cs) / sn) );


# [old]
- (digamma(1/n) + Euler)/n + sum(
	+ 1/n * (1 - cs) * log(x^2 - 2*cs*x + 1) +
	- 1/n * (1 - cs) * log(2 - 2*cs) +
	- 2/n * sn * atan((x - cs) / sn) +
	+ 2/n * sn * atan((1 - cs) / sn) );

###
n = sqrt(5)
integrate(\(x) 1 / (x^n + 1)^(1/n), 0, 1)
- (digamma(1/n) + Euler)/n +
	+ integrate(\(x) x^(n-2) * (x - 1) / (x^n - 1), 1, 2^(1/n))$value


### Note:
# - the fraction decomposition of x^p / (x^n - 1)
#   (for n, p = integers) is described in file:
#   Integrals.Fractions.Unity.R;


### Derivation:
# Limit: (p+1) -> 0
n = 3; p = -1 + 1E-6; k = 3;
integrate(\(x) 1 / (x^3 + 1)^(1/3)/x - 1/x * exp(-x), 0, Inf)
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n - gamma(p+1)
- (digamma(1/3) - 2*Euler)/3

#
integrate(\(x) 1/x * exp(-x) - 1/(x*(x+1)), 0, Inf)
- Euler

#
integrate(\(x) x / (x^2 + x + 1), 1, 2^(1/3))
1/2*log((2^(2/3) + 2^(1/3) + 1)/3) - 1/sqrt(3)*atan((2^(1/3) + 1/2)*2/sqrt(3)) +
	+ 1/sqrt(3)*atan((1 + 1/2)*2/sqrt(3))

# simplifications:
- 1/sqrt(3) * atan((2^(1/3) + 1/2) * 2/sqrt(3)) +
	+ 1/sqrt(3) * atan((1 + 1/2) * 2/sqrt(3));
- 1/sqrt(3) * atan((2^(1/3)-1)*2/sqrt(3) / (1 + 4/3*(2^(1/3) + 1/2)*(1 + 1/2)))
- 1/sqrt(3) * atan((2^(1/3) - 1) / (2^(1/3) + 1) / sqrt(3));


### Derivation: Pow = 5
integrate(\(x) x / (x^5 + 1)^(1/5), 0, 1)
1 - gamma(4/5) * gamma(2/5) / gamma(1/5) +
	- integrate(\(x) 1/(x^5 + 1)^(1/5) / x^2 - 1/x^2, 0, 1)$value;
1 - gamma(4/5) * gamma(2/5) / gamma(1/5) +
	+ integrate(\(x) x^3 * (x - 1) / (x^5 - 1)^(6/5), 1, 2^(1/5))$value;
# TODO: ???

### Derivation: on Arbitrary Interval;
# - for technique, see:
#   Integrals.Fractions.Unity.Radicals.Mixed.R;
lim = 3/4
integrate(\(x) x^2 / (x^5 + 1)^(3/5), 0, lim)
integrate(\(x) 1/5 / (x^2 * (x + 1)^3)^(1/5), 0, lim^5)
x = 1 / (1/lim^5 + 1)^(1/5);
integrate(\(x) x / (x^5 - 1), 1/x, Inf)
integrate(\(x) x^2 / (1 - x^5), 0, x)
cs = cos(c(2,4)*pi/5); csp = cos(3*c(2,4)*pi/5);
sn = sin(c(2,4)*pi/5); snp = sin(3*c(2,4)*pi/5);
-1/5 * log(1-x) - sum(csp*log(x^2 - 2*cs*x + 1) +
	- 2*snp*(atan((x - cs)/sn) + atan(cs/sn))) / 5;


# Fraction decomposition:
integrate(\(x) x^3 * (x-1) / (x^5 - 1), 1, 2^(1/5))
integrate(\(x) x^3 * (x-1) * (1/(x-1) +
	+ 2*(cos(2*pi/5)*x - 1) / (x^2 - 2*cos(2*pi/5)*x + 1) +
	+ 2*(cos(4*pi/5)*x - 1) / (x^2 - 2*cos(4*pi/5)*x + 1) ) / 5, 1, 2^(1/5));
cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
integrate(\(x) x^3 / 5 +
	+ x^3 * (x-1) * (2*(cs[1]*x - 1) / (x^2 - 2*cs[1]*x + 1) +
	+ 2*(cs[2]*x - 1) / (x^2 - 2*cs[2]*x + 1) ) / 5, 1, 2^(1/5));
integrate(\(x) x^3 / 5 +
	+ (2 * (cs[1]*x^3 - x^2 + 2*cs[1]^2*x^2 - cs[1]*x^2 +
		+ (4*cs[1]^3 - 2*cs[1]^2 - 3*cs[1] + 1)*x) +
		+ 2*(8*cs[1]^4 - 4*cs[1]^3 - 8*cs[1]^2 + 3*cs[1] + 1)) / 5 +
	# LOG:
	# + ((16*cs[1]^5 - 8*cs[1]^4 - 20*cs[1]^3 + 8*cs[1]^2 + 5*cs[1] - 1) *
	#	2*(x - cs[1]) / (x^2 - 2*cs[1]*x + 1) +
	+ ((1 - cs[1]) * 2*(x - cs[1]) / (x^2 - 2*cs[1]*x + 1) +
	# + 2*(cs[1]^2 - 1) * (16*cs[1]^4 - 8*cs[1]^3 - 12*cs[1]^2 + 4*cs[1] + 1) /
	#   (x^2 - 2*cs[1]*x + 1) + # Parenth == 1;
	+ 2*(cs[1]^2 - 1) / (x^2 - 2*cs[1]*x + 1) +
	+ 2*(cs[2]*x - 1) * x^4 / (x^2 - 2*cs[2]*x + 1) +
	- 2*(cs[2]*x^4 - x^3) / (x^2 - 2*cs[2]*x + 1) ) / 5, 1, 2^(1/5));
x = 2^(1/5);
sum(
	+ 1/5 * (1 - cs) * log(x^2 - 2*cs*x + 1) +
	- 1/5 * (1 - cs) * log(2 - 2*cs) +
	- 2/5 * sn * atan((x - cs) / sn) +
	+ 2/5 * sn * atan((1 - cs) / sn) );

# Note: == 0
x^4 / 20 - 1/20 + 2*(- 2/5 *(x^3/3 - x^2/2 - x - 1/3 + 1/2 + 1)) +
	+ sum(
	+ 2/5 * (cs*x^4/4 + 2*cs^2*x^3/3 +
		- (cs*x^3/3 + 2*cs^2*x^2/2 + 4*cs^3*x - 3*cs*x) +
		- 3*cs*x^2/2 + 4*cs^3*x^2/2 + (8*cs^4 - 8*cs^2)*x) +
	- 2/5 * (cs/4 + 2*cs^2/3 - (cs/3 + 2*cs^2/2 + 4*cs^3 - 3*cs) +
		- 3*cs/2 + 4*cs^3/2 + (8*cs^4 - 8*cs^2)) );

#
x = 2^(1/5)
integrate(\(x) (x^4 - x^2) / (x^5 - 1), 1, x)
cs = cos(c(2,4)*pi/5); sn = sin(c(2,4)*pi/5);
cs2 = cos(c(2,4)*3*pi/5); sn2 = sin(c(2,4)*3*pi/5);
cs4 = cos(c(2,4)*pi); sn4 = sin(c(2,4)*pi);
integrate(\(x) 1/5 * (x^3 + x^2) +
	+ (x^4 - x^2) * 2*(cs[1]*x-1) / (x^2 - 2*cs[1]*x + 1) / 5 +
	+ (x^4 - x^2) * 2*(cs[2]*x-1) / (x^2 - 2*cs[2]*x + 1) / 5, 1, x)
integrate(\(x)
	+ (2*(cs4[1]-cs2[1])*(x-cs[1]) - 2*(sn4[1]-sn2[1])*sn[1]) /
		(x^2 - 2*cs[1]*x + 1) / 5 +
	+ (2*(cs4[2]-cs2[2])*(x-cs[2]) - 2*(sn4[2]-sn2[2])*sn[2]) /
		(x^2 - 2*cs[2]*x + 1) / 5, 1, x)
sum((cs4-cs2) * log((x^2 - 2*cs*x + 1) / (2-2*cs)) +
	- 2*(sn4-sn2) * (atan((x - cs)/sn) - atan((1-cs)/sn)) ) / 5


#
integrate(\(x) 1/(x^5 + 1)^(1/5) / x^2 - (1/x^2 + 1/x) * exp(-x), 0, Inf)
n = 5; p = -2 + 1E-6; k = 5;
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n - gamma(p+1) - gamma(p+2)
1 - gamma(4/5) * gamma(2/5) / gamma(1/5)

#
integrate(\(x) (1/x^2 + 1/x) * exp(-x) - 1/x^2, 0, Inf)
-1;

# Unsolved:
integrate(\(x) -1 / (x^5 + 1)^(1/5) / x^2 + 1/x^2, 0, 1)
integrate(\(x) (1 - 1/(x^5 + 1)^(1/5)) / x^2, 0, 1)
integrate(\(x) 1/5 * (1 - 1/(x + 1)^(1/5)) / x^(6/5), 0, 1)
integrate(\(x) x^3 * (x - 1) / (x^5 - 1)^(6/5), 1, 2^(1/5))
2^(4/5) - 2 + integrate(\(x) (4*x^4 - 3*x^3) / (x^5 - 1)^(1/5), 1, 2^(1/5))$value;
2^(4/5) - 1 - integrate(\(x) 3*x^3 / (x^5 - 1)^(1/5), 1, 2^(1/5))$value;
2^(4/5) - 1 - integrate(\(x) 3/5*(x+1)^(-1/5) / x^(1/5), 0, 1)$value;
2^(4/5) - 1 - integrate(\(x) 3 * x^3 / (x^5 + 1)^(1/5), 0, 1)$value;

# TODO: ???

#
integrate(\(x) (4*x^4 - 3*x^3) / (x^5 - 1)^(1/5), 1, 2^(1/5))
(2 - 2^(4/5)) + integrate(\(x) x^3 * (x-1) / (x^5 - 1)^(6/5), 1, 2^(1/5))$value


# cyclic redundancy: x^2 / (...)^(1/5)
integrate(\(x) (3*x^4 - 2*x^3) /(x^5-1)^(2/5), 1, 2^(1/5))
2 - 2^(4/5) + integrate(\(x) 2 * x^3 * (x-1) / (x^5 - 1)^(7/5), 1, 2^(1/5))$value
#
integrate(\(x) 2 * x^3 * (x-1) / (x^5 - 1)^(7/5), 1, 2^(1/5))
2^(4/5) - 1 - integrate(\(x) 2 * x^3 / (x^5-1)^(2/5), 1, 2^(1/5))$value
2^(4/5) - 1 - integrate(\(x) 2 * x^2 / (x^5 + 1)^(1/5), 0, 1)$value


##################

### Special Cases:

### I( 1 / (x^(2*n) + 1)^(1/n) )

###
n = sqrt(5)
integrate(\(x) 1 / (x^(2*n) + 1)^(1/n), 0, 1)
gamma(1/(2*n))^2 / gamma(1/n) / (4*n)


### Series:

### I( x / (x^3 + 1)^(1/3) )
integrate(\(x) 1 / (x^(3/2) + 1)^(4/3), 0, 1)
gamma(2/3)^2/gamma(4/3) / 3
#
integrate(\(x) x / (x^3 + 1)^(4/3), 0, 1)
gamma(2/3)^2/gamma(4/3) / 6
#
integrate(\(x) x / (x^3 + 1)^(1/3), 0, 1)
1/2^(1/3) - gamma(2/3)^2/gamma(4/3) / 6

### I( 1 / (x^3 + 1)^(2/3) )
integrate(\(x) 1 / (x^3 + 1)^(2/3), 0, 1)
gamma(1/3)^2/gamma(2/3) / 6


### Trivial:
# see section further below;
integrate(\(x) 1 / (x^3 + 1)^(4/3), 0, 1)
1/2^(1/3)

###
integrate(\(x) 1 / (x^(3/2) + 1)^(5/3), 0, 1)
integrate(\(x) 2 * x / (x^3 + 1)^(5/3), 0, 1)
1/2^(2/3)
#
integrate(\(x) x / (x^3 + 1)^(5/3), 0, 1)
1/2^(5/3)


### Series: Pow = 5

###
integrate(\(x) x^2 / (x^5 + 1)^(6/5), 0, 1)
gamma(3/5)^2/gamma(6/5) / 10
#
integrate(\(x) x^2 / (x^5 + 1)^(1/5), 0, 1)
1/2^(6/5) - gamma(3/5)^2/gamma(6/5) / 20


###
integrate(\(x) x^3 / (x^5 + 1)^(8/5), 0, 1)
gamma(4/5)^2/gamma(8/5) / 10
#
integrate(\(x) x^3 / (x^5 + 1)^(3/5), 0, 1)
1/2^(3/5) - 3*gamma(4/5)^2/gamma(8/5) / 10


### Other Series:

### Type: (x^n + 1)^(2/n)
integrate(\(x) 1 / (x^6 + 1)^(1/3), 0, 1)
gamma(1/6)^2 / gamma(1/3) / 12

###
integrate(\(x) 1 / (x^(6/5) + 1)^(5/3), 0, 1)
integrate(\(x) 5*x^4 / (x^6 + 1)^(5/3), 0, 1)
gamma(5/6)^2 / gamma(5/3) / (6/5)/2

###
integrate(\(x) 1 / (x^8 + 1)^(1/4), 0, 1)
gamma(1/8)^2 / gamma(1/4) / 16


### Other with [0,1] == [1,Inf]

### I( sqrt(x^4+1) / x^2 )
integrate(\(x) sqrt(x^4+1) / x^2 - 1/x^2, 0, 1)
gamma(-1/4)*gamma(1/2 + 1/4) / gamma(1/2) / 4 + 1;

# =>

### I( x^2 / sqrt(x^4+1) )
integrate(\(x) x^2 / sqrt(x^4+1), 0, 1)
gamma(-1/4)*gamma(1/2 + 1/4) / gamma(1/2) / 8 + sqrt(2)/2;
# beta(3/4, 1/2 - 3/4) / 8 + sqrt(2)/2;

### I( x^2 * sqrt(x^4+1) )
integrate(\(x) x^2 * sqrt(x^4+1), 0, 1)
integrate(\(x) 2/5 * x^2 / sqrt(x^4+1), 0, 1)$value + sqrt(2)/5;
gamma(-1/4)*gamma(1/2 + 1/4) / gamma(1/2) / 20 + 2*sqrt(2)/5;

### Regular:

### I( 1 / sqrt(x^4+1) )
# Type: 1 / (x^(2*n) + 1)^(1/n);
integrate(\(x) 1 / sqrt(x^4+1), 0, 1)
beta(1/4, 1/4) / 8;

### I( sqrt(x^4+1) )
integrate(\(x) sqrt(x^4+1), 0, 1)
(beta(1/4, 1/4) / 8 + 1/sqrt(2)) * 2/3;


### Pow = 8

### I( 1/(x^8+1)^(1/4) )
integrate(\(x) 1/(x^8+1)^(1/4), 0, 1)
beta(1/8, 1/4 - 1/8) / 16;


### I( x^4 / (x^8+1)^(1/4) )

### I( (x^8+1)^(3/4) / x^4 )
integrate(\(x) (x^8+1)^(3/4) / x^4 - 1/x^4, 0, 1)
# beta(-3/8, -3/4 + 3/8) / 16 + 1/3;
gamma(-3/8) * gamma(-3/4+3/8) / gamma(-3/4) / 16 + 1/3;
# =>
integrate(\(x) x^4 / (x^8+1)^(1/4), 0, 1)
gamma(-3/8) * gamma(-3/4+3/8) / gamma(-3/4) / 32 + 2^(3/4) / 6;


### I( x^6 / (x^8+1)^(3/4) )

### I( (x^8+1)^(1/4) / x^2 )
integrate(\(x) (x^8+1)^(1/4) / x^2 - 1/x^2, 0, 1)
# beta(-1/8, -1/4 + 1/8) / 16 + 1;
gamma(-1/8) * gamma(-1/8) / gamma(-1/4) / 16 + 1;
# =>
integrate(\(x) x^6 / (x^8+1)^(3/4), 0, 1)
gamma(-1/8) * gamma(-1/8) / gamma(-1/4) / 32 + 2^(1/4) / 2;


#######################
#######################

### Integer Powers:

# TODO:

###
integrate(\(x) 1/(x^3 + 1)^5, 0, 1, rel.tol=1E-8)
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3 + 2/3 * int.FrU01(3, 0))));
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3))) +
	+ 11/12 * 8/9 * 5/6 * 2/3 * int.FrU01(3, 0);
# TODO:
# (3!/(3*2^4) + 11 * 2!/(3^2*2^3) + 11*8 * 1!/(3^3*2^2) + 11*8*5 * 0!(3^4*2)) / gamma(5)
1/2^4 / 12 + 1/2^3 * 11/12 / 9 + 1/2^2 * 11/12*8/9 / 6 + 1/2 * 11/12*8/9*5/6 / 3 +
	+ int.FrU01(3, 0) * gamma(5 - 1/3)/gamma(1 - 1/3)/gamma(5);

###
# I[3,6] = 1/2^5 / (3*5) + 14/15 * I[3,5];
integrate(function(x) 1/(x^3 + 1)^6, 0, 1, rel.tol=1E-8)
1/2^5 / 15 + 1/2^4 * 14/15 / 12 + 1/2^3 * 14/15*11/12 / 9 +
	+ 1/2^2 * 14/15*11/12*8/9 / 6 + 1/2 * 14/15*11/12*8/9*5/6 / 3 +
	+ int.FrU01(3, 0) * gamma(6 - 1/3)/gamma(1 - 1/3)/gamma(6);


#############

### [old]

### Radicals:

###
n = 5
I0 = integrate(\(x) sqrt(x^n + 1), 0, 1)
Ii = integrate(\(x) 1/sqrt(x^n + 1), 0, 1)
(n+2)/2*I0$value - n/2*Ii$value - sqrt(2) # == 0


###
n = 5
I1 = integrate(\(x) (x^n + 1)^(1/3), 0, 1)
I2 = integrate(\(x) (x^n + 1)^(2/3), 0, 1)
In1 = integrate(\(x) 1/(x^n + 1)^(1/3), 0, 1)
In2 = integrate(\(x) 1/(x^n + 1)^(2/3), 0, 1)
(n+3)/3*I1$value - n/3*In2$value - 2^(1/3) # == 0
(2*n+3)/3*I2$value - 2*n/3*In1$value - 2^(2/3) # == 0


### [old]
p = sqrt(5) - 2
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
IInf = gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n
print(IInf)

# on [0, 1]
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=1)$value
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=1, upper=Inf)$value
# TODO: find formula for -0.3248496;


#####################
#####################

### Diff Type

### I( x^p * (1 - x^n)^r )
### on [0, 1]

###
p = sqrt(3); r = sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^r, lower=0, upper=1)
gamma((p+1)/n) * gamma(r+1) / gamma((p+1)/n + r + 1) / n

### r > -1; p > - 1;
p = sqrt(7); r = - 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^r, lower=0, upper=1)
gamma((p+1)/n) * gamma(r+1) / gamma((p+1)/n + r + 1) / n


### Special Cases:

### I( (1 - x^n)^r )
r = sqrt(7); n = sqrt(5)
integrate(function(x) (1 - x^n)^r, lower=0, upper=1)
gamma(1/n) * gamma(r+1) / gamma(1/n + r + 1) / n


### Fractions:
# works for r < 2; (e.g. up to r = 1.95)
p = sqrt(7); r = sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p / (1 - x)^r - n^r * x^p / (1 - x^n)^r, lower=0, upper=1)
gamma(p+1) * gamma(1 - r) / gamma(p - r + 2) +
	- n^r * gamma((p+1)/n) * gamma(1 - r) / gamma((p+1)/n - r + 1) / n;


ff = \(x) sapply(x, \(r) {
	integrate(function(x) x^p / (1 - x)^r - n^r * x^p / (1 - x^n)^r, 0, 1)$value;
})
fg = \(r) gamma(p+1) * gamma(1 - r) / gamma(p - r + 2) +
	- n^r * gamma((p+1)/n) * gamma(1 - r) / gamma((p+1)/n - r + 1) / n;
curve(ff(x), 1.5, 1.95, col="green")
curve(fg(x), add=TRUE, col="red", lty=2)


### Lim: r -> 1
p = sqrt(7); n = 1/sqrt(5)
integrate(function(x) x^p / (1 - x) - n * x^p / (1 - x^n), lower=0, upper=1)
digamma((p+1)/n) - digamma((p+1)) + log(n);


### Special Variant:
p = sqrt(7); r = 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^(1-r) / (1 + x^n)^r, 0, 1)
gamma(1-r)	/ (2*n) * (gamma((p+1)/(2*n)) / gamma((p+1)/(2*n) - r + 1) +
	- gamma((p+1)/(2*n) + 1/2) / gamma((p+1)/(2*n) - r + 3/2) );
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n)

###
p = sqrt(7); r = 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 + x^n)^(1-r) / (1 - x^n)^r, 0, 1)
gamma(1-r)	/ (2*n) * (gamma((p+1)/(2*n)) / gamma((p+1)/(2*n) - r + 1) +
	+ gamma((p+1)/(2*n) + 1/2) / gamma((p+1)/(2*n) - r + 3/2) );
(beta(1-r, (p+1)/(2*n)) + beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n)


# Note:
# - Variants derived from (1 - x^n)^p / (1 + x^n)^q
#   are presented in a section below;


#####################

#############
### Other ###
#############

### Special Cases:

n = 8
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 9
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 10
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = sqrt(11)
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
integrate(function(x) 2/(2 - x^n)^(1+1/n), lower=0, upper=1)
# == 1


### Gen: I( 1 / (a - x^n)^(1/n + 1) )
a = sqrt(7); n = sqrt(5);
integrate(function(x) 1 / (a - x^n)^(1/n + 1), lower=0, upper=1)
1/(a-1)^(1/n) / a


### Gen: I( 1 / (x^n + a)^(1/n+1) )
n = sqrt(7); a = sqrt(5);
integrate(function(x) 1 / (x^n + a)^(1/n + 1), lower=0, upper=1)
1/(a+1)^(1/n) / a


### Gen: I( x^p / (x^n + a)^((p+1)/n + 1) )
p = sqrt(3); n = sqrt(7); a = sqrt(5);
integrate(function(x) x^p / (x^n + a)^((p+1)/n + 1), lower=0, upper=1)
1/(a+1)^((p+1)/n) / (a*(p+1))


##########################
##########################

### Mixed Types

### I( (1 - x^8)^(1/4) / (x^4 + 1) )
# Maths 505: An awesome calculus result I cooked up
# https://www.youtube.com/watch?v=WoxHJyOOma4
# Note: substitution tan(x) = y^n;

### Tan Pow: 1

### Gen: Pow 1
p = 2/5; n = 5;
# p = 3/5; n = 5/2; # x^0
integrate(\(x) x^(n*(1-p)-1) * (1 - x^(2*n))^p / (x^(2*n) + 1), 0, 1)
integrate(\(x) (1/tan(x) - tan(x))^p / n, 0, pi/4)
2^p / (4*n) * beta((p+1)/2, (1-p)/2)


### I( (1 - x^4)^(1/2) / (x^4 + 1) )
integrate(\(x) sqrt(1 - x^4) / (x^4 + 1), 0, 1)
integrate(\(x) sqrt(1/tan(x) - tan(x)) / 2, 0, pi/4)
pi/4

### I( x^2 * (1 - x^8)^(1/4) / (x^8 + 1) )
integrate(\(x) x^2 * (1 - x^8)^(1/4) / (x^8 + 1), 0, 1)
integrate(\(x) (1/tan(x) - tan(x))^(1/4) * 1/4, 0, pi/4)
2^(1/4)/16 * beta(5/8, 3/8)

### I( (1 - x^3)^(1/3) / (x^3 + 1) )
integrate(\(x) 1/2 * (1 - x^3)^(1/3) / (x^3 + 1), 0, 1)
integrate(\(x) x * (1 - x^6)^(1/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (1/tan(x) - tan(x))^(1/3), 0, pi/4)
2^(1/3)/12 * beta(2/3, 1/3)

### I( (1 - x^6)^(1/3) / (x^6 + 1) )
integrate(\(x) (1 - x^6)^(2/3) / (x^6 + 1), 0, 1)
integrate(\(x) 1/3 * (1/tan(x) - tan(x))^(2/3), 0, pi/4)
2^(2/3)/12 * beta(5/6, 1/6)


### Tan Pow: 2

# - can be generalized based on the cos/sin-integral:
#   the TAN-subintegral can be actually skipped;


### I( (1 - x^8)^(1/4) / (x^4 + 1) )
integrate(\(x) (1 - x^8)^(1/4) / (x^4 + 1), 0, 1)
integrate(\(x) (1/tan(x)^2 - tan(x)^2)^(1/4) / 2, 0, pi/4)
integrate(\(x) sqrt(2)/4 * (cos(x) / sin(x)^2)^(1/4), 0, pi/2)
sqrt(2)/8 * beta(5/8, 1/4)

### Derived:

### I( x^4 / (x^8 + 1)^(3/4) )
integrate(\(x) x^4 / (x^8 + 1)^(3/4), 0, 1)
sqrt(2)/16i * beta(5/8, 1/4) * sinh(1i*pi/8)
sqrt(2)/16 * beta(5/8, 1/4) * sin(pi/8)

### I( x^2 / (x^8 + 1)^(3/4) )
integrate(\(x) x^2 / (x^8 + 1)^(3/4), 0, 1)
beta(3/8, 3/8) / 16

### I( 1 / (x^8 + 1)^(3/4) )
integrate(\(x) 1 / (x^8 + 1)^(3/4), 0, 1)
gamma(1/8)*gamma(3/4 - 1/8) / gamma(3/4) / 8 +
	- sqrt(2)/16i * beta(5/8, 1/4) * sinh(1i*pi/8);
beta(1/8, 5/8) / 8 - sqrt(2)/16 * beta(1/4, 5/8) * sin(pi/8);

### I( (x^8 + 1)^(1/4) )
integrate(\(x) (x^8 + 1)^(1/4), 0, 1)
2^(1/4)/3 + 1/12*(beta(1/8, 5/8) - sqrt(2)/2 * beta(1/4, 5/8) * sin(pi/8));

### I( x^8 / (x^8 + 1)^(3/4) )
integrate(\(x) x^8 / (x^8 + 1)^(3/4), 0, 1)
2^(1/4)/3 - 1/24*(beta(1/8, 5/8) - sqrt(2)/2 * beta(1/4, 5/8) * sin(pi/8));

### I( x^6 / (x^8 + 1)^(7/4) )
integrate(\(x) x^6 / (x^8 + 1)^(7/4), 0, 1)
beta(7/8, 7/8) / 16

### I( x^6 / (x^8 + 1)^(3/4) )
integrate(\(x) x^6 / (x^8 + 1)^(3/4), 0, 1)
2^(-3/4) - 6*beta(7/8, 7/8) / 16

### I( 1 / (x^8 + 1)^(3/4) / x^2 )
integrate(\(x) 1 / (x^8 + 1)^(3/4) / x^2 - 1/x^2 + 1, 0, 1)
gamma(7/8) * gamma(3/4-7/8) / gamma(3/4) / 8 +
	+ 6*beta(7/8, 7/8) / 16 - 2^(-3/4) + 2;
-3/8 * beta(7/8, 7/8) - 2^(-3/4) + 2;


### Helper:
integrate(\(x) x^6 / (x^8 + 1)^(3/4) - 1, 0, Inf)
gamma(7/8) * gamma(3/4-7/8) / gamma(3/4) / 8


### I( 1 / (x^8 + 1)^(1/4) )
# see section above:
integrate(\(x) 1 / (x^8 + 1)^(1/4), 0, 1)
gamma(1/8)^2 / gamma(1/4) / 16
beta(1/8, 1/8) / 16

# Note: there are 2 trivial variants for even power;
integrate(\(x) x^4 / (x^8 + 1)^(5/4), 0, 1)
beta(5/8, 5/4 - 5/8) / 16



### Pow-Series: 1/8

###  I( (1 - x^16)^(1/8) / (x^8 + 1) )
integrate(\(x) (1 - x^16)^(1/8) / (x^8 + 1), 0, 1)
p = 0; n = 8; r = 7/8;
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n);
(beta(1/8, 1/16) - beta(1/8, 9/16)) / 16

### I( x^2 * (1 - x^16)^(1/8) / (x^8 + 1) )
integrate(\(x) x^2 * (1 - x^8)^(1/8) / (x^8 + 1)^(7/8), 0, 1)
2^(1/4-2)/4 * beta(9/16, 3/8)
(beta(1/8, 3/16) - beta(1/8, 11/16)) / 16

### Div: I( x^4 / (1 - x^16)^(1/8) / (x^8 + 1) )
integrate(\(x) x^4 / (1 - x^16)^(1/8) / (x^8 + 1), 0, 1)
2^(-1/4-2)/4 * beta(7/16, 5/8)
### Prod: I ( x^4 * (1 - x^16)^(1/8) / (x^8 + 1) )
integrate(\(x) x^4 * (1 - x^16)^(1/8) / (x^8 + 1), 0, 1)
p = 4; n = 8; r = 7/8;
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n);
(beta(1/8, 5/16) - beta(1/8, 13/16)) / 16

### I( x^6 * (1 - x^16)^(1/8) / (x^8 + 1) )
integrate(\(x) x^6 * (1 - x^16)^(1/8) / (x^8 + 1), 0, 1)
p = 6; n = 8; r = 7/8;
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n);
(beta(1/8, 7/16) - beta(1/8, 15/16)) / 16


### Pow-Series: 3/8

### I( (1 - x^16)^(3/8) / (x^8 + 1) )
integrate(\(x) (1 - x^8)^(3/8) / (x^8 + 1)^(5/8), 0, 1)
2^(3/4-2)/4 * beta(11/16, 1/8)

### Gen: I( x^p * (1 - x^16)^(3/8) / (x^8 + 1) )
p = sqrt(3)
integrate(\(x) x^p * (1 - x^16)^(3/8) / (x^8 + 1), 0, 1)
n = 8; r = 5/8;
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n);


# Based on Pow 3 & Higher:
integrate(\(x) x^6 * (1 - x^16)^(1/8) / (x^16 + 1), 0, 1)
2^(1/8-5) * beta(9/16, 7/16)

integrate(\(x) x^4 * (1 - x^16)^(1/8) / (x^16 + 1)^(3/4), 0, 1)
2^(3/8-5) * beta(9/16, 5/16)

# I( 1 / (x^16 + 1)^(1/8) )
integrate(\(x) 1 / (x^16 + 1)^(1/8), 0, 1)
beta(1/16, 1/8 - 1/16) / 32


integrate(\(x) x * (1 - x^8)^(3/8) / (x^8 + 1)^(7/8), 0, 1)
2^(1/2-4) * beta(11/16, 1/4)

integrate(\(x) (1 - x^8)^(1/8) / (x^8 + 1)^(3/8), 0, 1)
2^(-1/4-3) * beta(9/16, 1/8)

integrate(\(x) (1 - x^16)^(1/8) / x^2 - 1/x^2, 0, 1)
# 1 + 2^(-3/4 - 2) * beta(9/16, -1/8)
1 + 2^(-3/4 - 2) * gamma(9/16) * gamma(- 1/8) / gamma(7/16)


########

### (1 - x^12)^p * (x^6 + 1)^q

### I( (1 - x^12)^(1/3) / (x^6 + 1) )
integrate(\(x) (1 - x^12)^(1/3) / (x^6 + 1), 0, 1)
integrate(\(x) (1/tan(x)^2 - tan(x)^2)^(1/3) / 3, 0, pi/4)
integrate(\(x) 2^(2/3)/6 * (cos(x) / sin(x)^2)^(1/3), 0, pi/2)
2^(2/3)/12 * beta(2/3, 1/6)

### I( x^6 * (1 - x^12)^(1/3) / ((x^6 + 1) * (1 - x^12)) )
integrate(\(x) x^6 * (1 - x^12)^(-2/3) / (x^6 + 1), 0, 1)
integrate(\(x) (1/tan(x)^2 - tan(x)^2)^(-2/3) / 3, 0, pi/4)
integrate(\(x) 2^(-4/3)/6 * (cos(x) / sin(x)^2)^(-2/3), 0, pi/2)
2^(-1/3)/24 * beta(7/6, 1/6)


### I( (1 - x^12)^(1/3) )
integrate(\(x) (1 - x^12)^(1/3), 0, 1)
2^(2/3)/15 * beta(2/3, 1/6) * sin(5*pi/12) / cos(3*pi/12)
beta(1/12, 4/3) / 12

### I( (1 - x^12)^(-2/3) )
integrate(\(x) 1 / (1 - x^12)^(2/3), 0, 1)
2^(2/3)/12 * beta(2/3, 1/6) * sin(5*pi/12) / cos(3*pi/12)
beta(1/12, 1/3) / 12

### I( 1 / ((x^6 + 1) * (1 - x^12)^(2/3)) )
integrate(\(x) 1 / ((x^6 + 1) * (1 - x^12)^(2/3)), 0, 1)
2^(2/3)/12 * beta(2/3, 1/6) * sin(5*pi/12) / cos(3*pi/12) +
	- 2^(-1/3)/24 * beta(7/6, 1/6);
beta(1/12, 1/3) / 12 - 2^(-1/3)/24 * beta(7/6, 1/6);

### I( x^6 / (1 - x^12)^(2/3) )
integrate(\(x) x^6 / (1 - x^12)^(2/3), 0, 1)
2^(2/3)/24 * beta(2/3, 1/6) * sin(3*pi/12) / cos(1*pi/12)
beta(7/12, 1/3) / 12

### I( (1 - x^12)^(1/3) / (1 - x^6) )
integrate(\(x) (1 - x^12)^(1/3) / (1 - x^6), 0, 1)
2^(2/3)/12 * beta(2/3, 1/6) * (2*sin(5*pi/12) / cos(3*pi/12) - 1)

### I( x^6 * (1 - x^12)^(1/3) / (1 - x^6) )
integrate(\(x) x^6 * (1 - x^12)^(1/3) / (1 - x^6), 0, 1)
integrate(\(x) (1 - x^12)^(1/3) / (1 - x^6) - (1 - x^12)^(1/3), 0, 1)
2^(2/3)/12 * beta(2/3, 1/6) * (2*sin(5*pi/12) / cos(3*pi/12) - 1) +
	- beta(1/12, 4/3) / 12;


### I( 1 / (x^12 + 1)^(2/3) )
integrate(\(x) (x^12 + 1)^(1/3), 0, 1)

integrate(\(x) 1 / (x^12 + 1)^(2/3), 0, 1)
integrate(\(x) 5/4 * (x^12 + 1)^(1/3) - 2^(1/3-2), 0, 1)
integrate(\(x) (1 - x)^(-5/12) * x^(-11/12) * 1/12, 0, 1/2)

# TODO


### I( x^4 * (1 - x^12)^(-5/12) )
integrate(\(x) x^4 * (1 - x^12)^(-5/12), 0, 1)
integrate(\(x) 1/6 * tan(x)^(-1/6), 0, pi/2)
integrate(\(x) 1/6 * x^(-1/6) / (1 + x^2), 0, Inf)
beta(5/12, 7/12) / 12;


### Various intervals:
integrate(\(x) x^4 / (1 - x^12)^(5/12), 0, 2^(-1/12))
integrate(\(x) 1/6 * tan(x)^(-1/6), 0, pi/4)
integrate(\(x) 1/6 * x^(-1/6) / (1 + x^2), 0, 1)
(digamma(17/24) - digamma(5/24)) / 24
# Derived =>
integrate(\(x) x^(-1) / (x^12 - 1)^(5/12), 2^(1/12), Inf)
# =>
integrate(\(x) x^(-1) / (x^12 - 1)^(5/12), 1, 2^(1/12))
beta(5/12, 7/12) / 12 +
	- (digamma(17/24) - digamma(5/24)) / 24;

### I( x^(-1) / (1 - x^12)^(5/12) )
integrate(\(x) x^(-1) / (1 - x^12)^(5/12) - 1/x, 0, 1)
integrate(\(x) x^4 / (x^12 - 1)^(5/12) - 1/x, 1, Inf)
integrate(\(x) 1/12 * x^(-1) * (1 - x)^(-5/12) - 1/12/x, 0, 1)
- (digamma(7/12) + Euler) / 12;
pi*(sqrt(3) - 2)/24 + log(3)/8 + log(2) / 4 - sqrt(3)/6 * acoth(sqrt(3))

# Limit:
eps = 1E-5;
gamma(7/12) * (gamma(eps) / gamma(7/12 + eps) +
	- 1/gamma(7/12) / eps) / 12
(beta(eps, 7/12) - 1/eps) / 12
# library(Rmpfr)
eps = mpfr("1E-12", 240)
gamma(7/12) * (gamma(eps) / gamma(mpfr(7, 240)/12 + eps) +
	- 1/gamma(mpfr(7, 240)/12) / eps) / 12;

#
ff = \(x) x^4 / abs(1 - x^12)^(5/12) - 1/x + exp(-x)/x;
integrate(ff, 0, 1, rel.tol=1E-8)$value +
integrate(ff, 1, Inf, rel.tol=1E-8)$value;
(beta(5/12, 7/12) - digamma(7/12) - Euler) / 12 - Euler;

#
integrate(\(x) x^4 / (1 - x^12)^(5/12) - 1/x + exp(-x)/x, 0, 1)
beta(5/12, 7/12) / 12 - pracma::expint(1) - Euler;


#############
### Gen Tan^2
p = 1/8; n = 4;
integrate(\(x) x^(n*(1-2*p)-1) * (1 - x^(2*n))^p / (x^(2*n) + 1)^(1-p), 0, 1)
integrate(\(x) 2^(2*p)/n * ((cos(x)^2 - sin(x)^2) / sin(2*x)^2)^p, 0, pi/4)
integrate(\(x) 2^(2*p-1)/n * (cos(x) / sin(x)^2)^p, 0, pi/2)
2^(2*p-2)/n * beta((p+1)/2, 1/2 - p)


#######
### Tan Pow: 3

### Gen:
p = 1/5; n = 5/2;
integrate(\(x) x^(n*(1-3*p)-1) * (1 - x^(2*n))^p / (x^(2*n) + 1)^(1-2*p), 0, 1)
integrate(\(x) 1/n * (1/tan(x)^3 + 1/tan(x) - tan(x) - tan(x)^3)^p, 0, pi/4)
integrate(\(x) 2^(3*p-1) / n * (cos(x) / sin(x)^3)^p, 0, pi/2)
2^(3*p-2) / n * beta((p+1)/2, (1-3*p)/2)

### I( (1 - x^5)^(1/5) / (x^5 + 1)^(3/5) )
integrate(\(x) (1 - x^5)^(1/5) / (x^5 + 1)^(3/5), 0, 1)
integrate(\(x) 2/5 * (1/tan(x)^3 + 1/tan(x) - tan(x) - tan(x)^3)^(1/5), 0, pi/4)
2^(-2/5) / 5 * beta(3/5, 1/5)


### Variant: cos^2
p = 1/5; n = 5/2;
integrate(\(x) x^(n*(1-3*p)-1) * (1 - x^(2*n))^(2*p) / (x^(2*n) + 1)^(1-p), 0, 1)
integrate(\(x) 1/n * (1/tan(x)^3 - 1/tan(x) - tan(x) + tan(x)^3)^p, 0, pi/4)
integrate(\(x) 2^(3*p-1) / n * (cos(x)^2 / sin(x)^3)^p, 0, pi/2)
2^(3*p-2) / n * beta(p+1/2, (1-3*p)/2)

### I( (1 - x^5)^(2/5) / (x^5 + 1)^(4/5) )
integrate(\(x) (1 - x^5)^(2/5) / (x^5 + 1)^(4/5), 0, 1)
integrate(\(x) 2/5 * (1/tan(x)^3 - 1/tan(x) - tan(x) + tan(x)^3)^(1/5), 0, pi/4)
2^(-2/5) / 5 * beta(7/10, 1/5)
# TODO: DONE?


#######
### Tan Pow: 4

### Gen:
p = 1/5; n = 5;
integrate(\(x) x^(n*(1-4*p)-1) * (1 - x^(2*n))^p / (x^(2*n) + 1)^(1-3*p), 0, 1)
integrate(\(x) 1/n * (1/tan(x)^4 + 2/tan(x)^2 +
	- 2*tan(x)^2 - tan(x)^4)^p, 0, pi/4)
integrate(\(x) ((cos(x)^2 - sin(x)^2) * (cos(x)^2 + sin(x)^2)^3)^p *
	2^(4*p) / n / sin(2*x)^(4*p), 0, pi/4)
integrate(\(x) 2^(4*p-1) / n * (cos(x) / sin(x)^4)^p, 0, pi/2)
2^(4*p-2) / n * beta((p+1)/2, (1-4*p)/2)

### I( (1 - x^10)^(1/5) / (x^10 + 1)^(2/5) )
# p = 1/5; n = 5;
integrate(\(x) (1 - x^10)^(1/5) / (x^10 + 1)^(2/5), 0, 1)
2^(-6/5) / 5 * beta(3/5, 1/10)


### Variant: cos^3
p = 1/5; n = 5;
integrate(\(x) x^(n*(1-4*p)-1) * (1 - x^(2*n))^(3*p) / (x^(2*n) + 1)^(1-p), 0, 1)
integrate(\(x) ((cos(x)^2 - sin(x)^2)^3 * (cos(x)^2 + sin(x)^2))^p *
	2^(4*p) / n / sin(2*x)^(4*p), 0, pi/4)
integrate(\(x) 2^(4*p-1) / n * (cos(x)^3 / sin(x)^4)^p, 0, pi/2)
2^(4*p-2) / n * beta((3*p+1)/2, (1-4*p)/2)

### I( (1 - x^10)^(3/5) / (x^10 + 1)^(4/5) )
# p = 1/5; n = 5;
integrate(\(x) (1 - x^10)^(3/5) / (x^10 + 1)^(4/5), 0, 1)
2^(-6/5) / 5 * beta(4/5, 1/10)


###################

### Generalization:
p = 1/8; q = 5/8; n = 3;
integrate(\(x) x^(n*(q-p)-1) * (1 - x^(2*n))^p / (x^(2*n) + 1)^q, 0, 1)
integrate(\(x) 2^(p-q+1) / n * (cos(x)^2 - sin(x)^2)^p / sin(2*x)^(1+p-q), 0, pi/4)
integrate(\(x) 2^(p-q) / n * cos(x)^p / sin(x)^(1+p-q), 0, pi/2)
2^(p-q-1) / n * beta((p+1)/2, (q-p)/2)


### Generalization: Simple Case (1-r) vs r;
# - see a section above for details;
p = sqrt(7); r = 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^(1-r) / (1 + x^n)^r, 0, 1)
gamma(1-r)	/ (2*n) * (gamma((p+1)/(2*n)) / gamma((p+1)/(2*n) - r + 1) +
	- gamma((p+1)/(2*n) + 1/2) / gamma((p+1)/(2*n) - r + 3/2) );
(beta(1-r, (p+1)/(2*n)) - beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n)

###
p = sqrt(7); r = 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 + x^n)^(1-r) / (1 - x^n)^r, 0, 1)
gamma(1-r)	/ (2*n) * (gamma((p+1)/(2*n)) / gamma((p+1)/(2*n) - r + 1) +
	+ gamma((p+1)/(2*n) + 1/2) / gamma((p+1)/(2*n) - r + 3/2) );
(beta(1-r, (p+1)/(2*n)) + beta(1-r, (p+1)/(2*n) + 1/2)) / (2*n)


####################
####################

# Note: moved "Mixed Radicals" to file:
# Integrals.Fractions.Unity.Radicals.Mixed.R;
