########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration
## Polynomial Fractions: Unity
## Definite Integrals: Difference-Type
##
## draft v.0.1b



####################

### Helper Functions

constEuler = 0.57721566490153286060651209008240243079;
Euler = constEuler;


### I on [0, 1]

# - code based on the Digamma function:
#   enables continuous values of n & p;
int.FrU01 = function(n, p=0) {
	# I( x^p / (x^n + 1))
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}

### Diff Type
int.FrDU01 = function(n, p=0) {
	# I( x^p/(1 - x) -  n*x^p / (1 - x^n) )
	if( p!= 0) {
		# r = int.FrDU01(1, p) - int.FrDU01(n, p);
		r = digamma((p+1)) - digamma((p+1)/n) - log(n);
		return(r);
	}
	digamma(1/n) + Euler + log(n);
}
int.FrDUp1 = function(n, p=0) {
	# integrate(\(x) (1 - x^p) / (1 - x^n), 0, 1)
	(digamma((p+1)/n) - digamma(1/n)) / n;
}


####################
####################

#################
### Diff Type ###
#################

### Simple Case:
# - see section below;
p = sqrt(3)
integrate(\(x) (1 - x^p) / (1-x), 0, 1)
(digamma(p) - digamma(1)) + 1/p;
(digamma(p+1) - digamma(1));

### Gen: I( (1 - x^p) / (1 - x^n) )
p = sqrt(3); n = sqrt(5)
integrate(\(x) (1 - x^p) / (1 - x^n), 0, 1)
(digamma((p+1)/n) - digamma(1/n)) / n;

### Gen: I( x^p * (1 - x^q) / (1 - x^n) )
p = sqrt(3); q = sqrt(2); n = sqrt(5)
integrate(\(x) x^p * (1 - x^q) / (1 - x^n), 0, 1)
(digamma((p+q+1)/n) - digamma((p+1)/n)) / n;


### Derived: log
p = sqrt(3); q = sqrt(2); n = sqrt(5)
integrate(\(x) x^p * (1 - x^q) * log(x) / (1 - x^n), 0, 1)
(pracma::psi(1, (p+q+1)/n) - pracma::psi(1, (p+1)/n)) / n^2;


### [old]

### 1 / (x^n - 1)

# I() on [0, 1]
n = 5
#
id = seq(floor((n - 1)/2))
cs = cos(2*pi*id/n); # ONLY 1 * cos()!
sn = sin(2*pi*id/n);
sn2 = sin(pi*id/n);
integrate(\(x) (x^3 + 2*x^2 + 3*x + 4) / (x^4 + x^3 + x^2 + x + 1), 0, 1)
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)

###
n = 7; # n = integer: odd & even!
#
id = seq(floor((n - 1)/2))
cs = cos(2*pi*id/n);
sn = sin(2*pi*id/n);
sn2 = sin(pi*id/n);
integrate(\(x) (x^5 + 2*x^4 + 3*x^3 + 4*x^2 + 5*x + 6) /
	(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1), 0, 1)
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)


# - cs * log(x^2 - 2*cs*x + 1) + 2*sn*atan((x - cs)/sn)
- sum(cs * log(2 - 2*cs)) + 2*sum(sn*atan((1 - cs)/sn) - sn*atan(- cs/sn))
log(2)/2 - sum(cs * log(1 - cs)) + 2*sum(sn*atan((1 - cs)/sn) - sn*atan(- cs/sn))
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)


### Fraction Decomposition:
x = sqrt(11)
#
n = 5
1/(x - 1) - n/(x^n - 1)
(x^3 + 2*x^2 + 3*x + 4) / (x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1))
#
n = 6
1/(x - 1) - n/(x^n - 1)
(x^4 + 2*x^3 + 3*x^2 + 4*x + 5) / (x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1)) + 1/(x+1)
#
n = 7
1/(x - 1) - n/(x^n - 1)
(x^5 + 2*x^4 + 3*x^3 + 4*x^2 + 5*x + 6) / (x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1))
#
n = 8
1/(x - 1) - n/(x^n - 1)
(x^6 + 2*x^5 + 3*x^4 + 4*x^3 + 5*x^2 + 6*x + 7) /
	(x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1)) + 1/(x+1)


##########################
##########################

Euler = 0.57721566490153286060651209008240243079;

### Digamma Function
# see: Lines That Connect: Extending the Harmonic Numbers to the Reals
# https://www.youtube.com/watch?v=9p_U_o1pMKo


###
n = 7
integrate(\(x) (1 - (1-x)^n)/x, 0, 1)
integrate(\(x) (1 - x^n) / (1 - x), 0, 1)
sum(1/seq(n))
digamma(n + 1) + Euler


###
n = sqrt(5)
integrate(\(x) (1 - (1-x)^n)/x, 0, 1)
integrate(\(x) (1 - x^n) / (1 - x), 0, 1)
digamma(n + 1) + Euler


####################
####################

### Diff Type

### Any n:
n = sqrt(7)
integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)
digamma(1/n) + Euler + log(n)


### I( x^p/(x-1) - n*x^p/(x^n - 1) )
p = sqrt(3); n = sqrt(11);
# can be: p > n;
integrate(\(x) x^p/(x-1) - n*x^p/(x^n - 1), 0, 1)
int.FrDU01(n, p)

### [old]

# [alternative] [n = 5]
pi/tan(pi/n)/2 + 1/4*log(n) + sqrt(5)/2*atanh(1/sqrt(5))

### ODD n:
n = 5
id = seq(1, floor(n/2))
cs = cos(2*pi*id/n)
sn = sin(2*pi*id/n);sh = sin(pi*id/n);
integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) - sn * atan(-cs/sn))

sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) - sn * atan(-cs/sn))
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) + sn * (pi/2 - 2*id*pi/n))
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn)) +
	+ 4*pi * sum(id*sn)/n - pi * sum(sn)
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn)) +
	+ pi/sin(pi/n) - pi/tan(pi/n/2)/2
2*sum(cs * log(sh)) - log(2) +
	+ pi/sin(pi/n)/2 - pi/tan(pi/n/2)/2
digamma(1/n) + Euler + log(n) +
	+ pi/tan(pi/n)/2 - pi/tan(pi/n/2)/2 + pi/sin(pi/n)/2; # sum(TRIG) == 0;

# Helper:
pi * sum(sn)
pi/tan(pi/n/2)/2
#
sum(id*sn)
n/sin(pi/n)/4


################

### Higher Power

### I( x^p / (x^n - 1)^2 )
p = sqrt(3); n = sqrt(7);
# can be: p > n;

# if: (p - n) != integer or > 0:
integrate(\(x) x^p / (x^n - 1)^2 - x^p * (1/(x-1)^2 - (n-1)/(x-1)) / n^2, 0, 1)
(digamma((p-n+1)/n) - digamma(p-n+1) + log(n)) * (p-n+1)/(n^2) +
	+ (3*n-1) / (2*n^2) + (0) / (2*p*n^3) # - (.../p + ...)/(2*n^2)

# TODO: compute the missing term;

# if: (p - n) == integer < 0
p = 0; n = sqrt(7);
integrate(\(x) x^p / (x^n - 1)^2 - x^p * (1/(x-1)^2 - (n-1)/(x-1)) / n^2, 0, 1)
(digamma((p-n+1)/n) + Euler + log(n)) * (p-n+1)/(n^2) +
	+ (3*n-1) / (2*n^2) + ifelse(p == 0, 0, (0) / (2*p*n^3))

# TODO: works only for p = 2!
ff = \() (integrate(\(x) x^p / (x^n - 1)^2 - x^p * (1/(x-1)^2 - (n-1)/(x-1)) / n^2, 0, 1)$value
	- (digamma((p-n+1)/n) + Euler + log(n)) * (p-n+1)/(n^2) +
	- (3*n-1)/2/n^2 - p*((5-p)*n - 9)/(4*n^2)) * (n^2)

p = 2; # p = 3;
n = 4; ff();
n = 5; ff();
n = 6; ff();
n = 7; ff();

x0 = c(1,2,3)
solve(cbind(x0^2, x0, 1), c(4,6,6 + 1 + 1/4 + 1/16 + 1/64)/ x0) * 96

###
p = 1/2; n = 3;
integrate(\(x) x^p / (x^n - 1)^2 - x^p * (1/(x-1)^2 - (n-1)/(x-1)) / n^2, 0, 1)
(8 - 3*log(3)) / 18
(2 - 3*log(3)) /(2*n^2) + 3/n^2


# Base:
p = sqrt(3); n = sqrt(7);
integrate(\(x) x^(p-n) / (x^n - 1) - x^(p-n) / (x-1) / n, 0, 1)
- (digamma(p-n+1) - digamma((p-n+1)/n) - log(n)) / n;


### Lim: x -> 1
n = sqrt(5)
x = 1 - 1/2^10
1/ (x^n - 1)^2 -  (1/(x-1)^2 - (n-1)/(x-1)) / n^2;
(5*n^2 - 6*n +1) / (12*n^2)

### Lim: x -> 1
x = 1 - 1/2^10
(1/(x^n - 1) - 1/n / (x-1))
1/(2*n) - 1/2;

