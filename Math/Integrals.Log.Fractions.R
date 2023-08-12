########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Fractions
###
### draft v.0.2d


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


### Helper Constants
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


### Solved:
# I( x^p * log(x) / (x^n + 1) ) on [0, Inf];
# I( x^p * log(x) / (x^n + 1)^k ) on [0, Inf]
# I( x^p * log(x^n + 1) / (x^n + 1)^k ) on [0, Inf]
# I( log(x^n + 1) / x^p ) on [0, Inf]
# I( log(x^n + 1) / x^p ) on [0, 1]
# I( log(x^n + a^n) / (x^n + b^n) ) on [0, Inf]
# I( (x^p - 1) / ((x^n + 1) * log(x)) ) on [0, 1]
### Diff-Type:
# I( x^p * log(x)^2 / (x^n - 1) ) on [0, Inf]
# - any n, p, k;


### Refactor:
# - I( log(x) / (x + 1)^p ) on [0, Inf]
#   moved to file Integrals.Log.Fractions.Simple.R;
# - I( log(x)^s / (x^n + 1) ) on [0, Inf]
#   moved to file: Integrals.Log.Fractions.Higher.R;
# - moved old / simplistic derivations to file:
#   Integrals.Log.Fractions.Old.R;


####################

#################
### Fractions ###
#################

### Feynman trick & Other tricks
# - usually (much) simpler than Contour integration;

### Generalization

### Base-Formulas:

### I( log(x) / (x^n + 1) )
n = 7
integrate(function(x) log(x) / (x^n + 1), 0, Inf)
- pi^2*cos(pi/n) / sin(pi/n)^2 / n^2

### I( x^p * log(x) / (x^n + 1) )
n = 8
p = sqrt(2)
integrate(function(x) x^p * log(x) / (x^n + 1), 0, Inf)
- pi^2*cos(pi*(p+1)/n) / sin(pi*(p+1)/n)^2 / n^2

### I( x^p * log(x) / (x^n + 1)^k )
n = 5; p = sqrt(2); k = sqrt(3);
integrate(function(x) x^p * log(x) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma((p+1)/n) - digamma(k - (p+1)/n)) / gamma(k) / n^2;


### on [0, 1]

### I( x^p * log(x) / (x^n + 1) )
n = 8
p = sqrt(2)
integrate(function(x) x^p * log(x) / (x^n + 1), 0, 1)
(pracma::psi(1, ((p+1)/n + 1)/2) - pracma::psi(1, (p+1)/n/2)) / (4*n^2);


### I( log(x^n + 1) / x^(p+1) )
n = sqrt(5); p = sqrt(3);
integrate(function(x) log(x^n + 1) / x^(p+1), 0, 1)
(digamma(((n-p)/n + 1)/2) - digamma((n-p)/(2*n))) / (2*p) - log(2)/p;


### I( x^p * (1 - x)^q * log(x) )
# - is equivalent to: I( x^p1 * log(x) / (x^n + 1)^k ) on [0, Inf];
p = sqrt(2); q = sqrt(3);
integrate(\(x) x^p * (1 - x)^q * log(x), 0, 1, rel.tol=1E-8)
gamma(p+1) * gamma(q+1) / gamma(p+q+2) * (digamma(p+1) - digamma(p+q+2))


####################
####################

### Basic Fractions:

Catalan = 0.915965594177219015054603514;


### I( x^p * log(x^n + 1) / (x^n + 1)^k )
p = sqrt(3); n = sqrt(5); k = sqrt(3) + sqrt(2);
integrate(function(x) x^p * log(x^n + 1) / (x^n + 1)^k, lower=0, upper=Inf)
- gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k - (p+1)/n) - digamma(k)) / gamma(k) / n;


### Lim: k -> 0
p = - sqrt(3); n = 5;
integrate(function(x) x^p * log(x^n + 1), lower=0, upper=Inf)
- gamma((p+1)/n) * gamma(-(p+1)/n) / n

#
k = 1E-6;
- gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k - (p+1)/n) - digamma(k)) / gamma(k) / n;


### Derivation:

### Basic / Helper
b = sqrt(3)
integrate(function(x) log(x) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b)/(2*b)
#
integrate(function(x) log(x) / (b^2*x^2 + 1), lower=0, upper=Inf)
- pi*log(b)/(2*b)

###
# - using Feynman's trick;
# - Note: Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)
integrate(function(x) log(x + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)/4 + Catalan;

### TODO:
integrate(function(x) log(x + 1) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b^2 + 1)/(4*b) + log(b)*atan(1/b)/b +
	- integrate(function(x) log(x) / (x^2 + b^2), 0, 1)$value;
pi*log(b^2 + 1)/(4*b) +
	- integrate(function(x) log(x) / (x^2 + 1), 0, 1/b)$value / b;

### TODO:
a = sqrt(5); b = sqrt(3)
integrate(function(x) log(x + a) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(a^2 + b^2)/(4*b) + log(b)*atan(a/b)/b +
	- integrate(function(x) log(x) / (x^2 + b^2), 0, a)$value;

### TODO:
integrate(function(x) x^(1/2) * log(x + 1) / (x^2 + 1), 0, 1)
integrate(function(x) 2 * x^2 * log(x^2 + 1) / (x^4 + 1), 0, 1)
integrate(\(b) 4*b^2 / (b^4 + 1) * atan(b), 0, 1)$value +
	+ (pracma::psi(1, 7/8) - pracma::psi(1, 3/8)) / 16 +
	- pi * (digamma(7/8) - digamma(3/8)) / 4 +
	+ (digamma(5/8) - digamma(1/8)) / 4 * pi/4 +
	+ (digamma(7/8) - digamma(3/8)) / 8 * log(2);


# Fraction Decomposition:
p = 1/2; b = sqrt(3);
# possibly solvable for p = rational?
# p = 1/2 => (x => x^2)
integrate(\(x) x^p * (1/(x + b) - (x - b)/(x^2 + 1)) / (b^2 + 1), 0, 1)
2/(b^2 + 1) * integrate(\(x) x^2 / (x^2 + b) - x^4 / (x^4 + 1) + b * x^2 / (x^4 + 1), 0, 1)$value
2/(b^2 + 1) * integrate(\(x) -b / (x^2 + b) + 1 / (x^4 + 1) + b * x^2 / (x^4 + 1), 0, 1)$value;
- 2*b / (b^2 + 1) * integrate(\(x) 1 / (x^2 + b), 0, 1)$value +
	+ (digamma(5/8) - digamma(1/8)) / (4*(b^2 + 1)) +
	+ (digamma(7/8) - digamma(3/8)) * b / (4*(b^2 + 1));
2*sqrt(b) / (b^2 + 1) * atan(sqrt(b)) - pi*sqrt(b) / (b^2 + 1) +
	+ (digamma(5/8) - digamma(1/8)) / (4*(b^2 + 1)) +
	+ (digamma(7/8) - digamma(3/8)) * b / (4*(b^2 + 1));


#######################

### I( log(x^n + a^n) / (x^n + b^n) )

#############
### log(P[2])
# - using Feynman's trick;
b = sqrt(3); a = sqrt(2)
integrate(function(x) log(x^2 + a^2) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(a + b)/b

### I( log(x^2 + a^2) / (x^2 + b^2)^2 )
b = sqrt(3); a = sqrt(2)
integrate(function(x) log(x^2 + a^2) / (x^2 + b^2)^2, lower=0, upper=Inf)
pi*(log(a + b) - b/(a+b)) / (2*b^3)

### I( log(x^2 + a^2) / (x^2 + b^2)^3 )
b = sqrt(3); a = sqrt(2)
integrate(function(x) log(x^2 + a^2) / (x^2 + b^2)^3, lower=0, upper=Inf)
pi*(3*log(a + b) - 4*b/(a+b) + a*b/(a+b)^2) / (8*b^5);


### Helper
integrate(function(x) log(x^2 + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)

###
b = sqrt(3)
integrate(function(x) log(x^2 + 1) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b+1)/b

### [Special Case]
b = sqrt(3)
integrate(function(x) log(x^2 + b^2) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(2*b)/b

# "cyclic redundancy" / Extra Cases
b = sqrt(3)
integrate(function(x) log(x^2 + 1) / (x^2 + b^2), lower=0, upper=Inf)
integrate(function(x) log(x^2 + 1) / (1 + b^2*x^2), lower=0, upper=Inf)$value +
	+ pi*log(b)/b;
integrate(function(x) log(x^2 + b^2) / (x^2 + 1), lower=0, upper=Inf)$value / b;
pi*log(b+1)/b


### Generalization:
# - using Feynman's trick;
# - using formulas for Fractions with roots of unity, see:
#   Integrals.Fractions.Unity.R;
n = 7
a = sqrt(3); b = sqrt(5)
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)
#
Int   = integrate(function(x) 1/(1 - x^n), lower=0, upper=a/b)$value
Const = - pi/sin(pi/n) / b^(n-1) * (pi*cos(pi/n)/sin(pi/n) / n - log(b));
(log(abs(a^n - b^n)) - n*log(b) + n * Int) * pi/sin(pi/n) / n / b^(n-1) + Const;


###################

### Explicit Cases:

### Case: n = 3
n = 3; # n = hard-coded;
a = sqrt(3); b = sqrt(5)
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)
sqrt(3)*pi/3 * log((a^n - b^n)/(a - b)) / b^2 - pi^2/3/b^2 +
	+ 2*pi/3 * atan((2*a/b + 1)/sqrt(3)) / b^2;


# Helper:
integrate(function(x) log(x) / (x^3 + 1), lower=0, upper=Inf)
- 2*pi^2/27

#
b = sqrt(5)
integrate(function(x) log(x) / (x^3 + b^3), lower=0, upper=Inf)
- 2*pi^2/27/b^2 + 2*pi/(3*sqrt(3))*log(b)/b^2

# Derivation:
- 2*pi^2/27/b^2 + log(b)/b^2 * integrate(function(x) 1 / (x^3 + 1), lower=0, upper=Inf)$value
# - 2*pi^2/27/b^2 + log(b)/b^2/3 * (log(x+1) - 1/2*log(x^2-x+1) + 3/sqrt(3)*atan((2*x-1)/sqrt(3)))
- 2*pi^2/27/b^2 + log(b)/b^2*(pi/2 - atan(-1/sqrt(3)))/sqrt(3)

# dI: evaluated at Inf & at 0;
n = 3
integrate(function(x) n*a^(n-1)/(b^n-a^n)* (1/(x^n + a^n) -  1/(x^n + b^n)), lower=0, upper=Inf)
# [not run]
1/2/(b^n - a^n) * log((x+a)^2/(x^2-a*x+a^2)) +
	- 1/2*a^2/b^2/(b^n - a^n)*log((x+b)^2/(x^2-b*x+b^2)) +
	+ 3/sqrt(3)/(b^n - a^n)*(atan(2*x/(a*sqrt(3)) - 1/sqrt(3)) - a^2/b^2*atan(2*x/(b*sqrt(3)) - 1/sqrt(3)))
3/sqrt(3)/(b^n - a^n)*(atan(2*x/(a*sqrt(3)) - 1/sqrt(3)) - a^2/b^2*atan(2*x/(b*sqrt(3)) - 1/sqrt(3)))
# =>
sqrt(3)*(1 - a^2/b^2)/(b^n - a^n)*(pi/2 + atan(1/sqrt(3)))
2*sqrt(3)*pi/3 * (1 - a^2/b^2)/(b^n - a^n)
pi/sin(pi/3) * (1 - (a/b)^2) / (b^n - a^n)

# back-Integration:
2*sqrt(3)*pi/3 * integrate(function(x) (x^2/b^2 - 1)/(x^n - b^n), lower=0, upper=a)$value +
	- 2*pi^2/9/b^2 + 2*pi/sqrt(3)*log(b)/b^2;
2*sqrt(3)*pi/9 * log((a^n - b^n)/(a - b)) / b^2 +
	- 2*pi^2/9/b^2 + 2*sqrt(3)*pi/9*log(b) / b^2 +
	+ 2*sqrt(3)*pi/9/b^2 * integrate(function(x) (x + 2*b)/(x^2+b*x+b^2), lower=0, upper=a)$value
sqrt(3)*pi/3 * log((a^n - b^n)/(a - b)) / b^2 - 2*pi^2/9/b^2 +
	+ 2*pi/3 *(atan((2*a/b + 1)/sqrt(3)) - atan(1/sqrt(3))) / b^2;


### Case: n = 5
n = 5
a = sqrt(3); b = sqrt(5);
#
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)
Int = integrate(function(x) 1/(1 - x^n), lower=0, upper=a/b)$value
Const = - pi/sin(pi/n) / b^(n-1) * (pi*cos(pi/n)/sin(pi/n) / n - log(b));
(log(abs(a^n - b^n)) - n*log(b) + n * Int) * pi/sin(pi/n) / n / b^(n-1) + Const;

# Int = can be computed exactly for n = integer;
# if(a > b): Int = Int[0, 1 - eps] + Int[1 + eps, a/b];

### Special Case: a == b
integrate(function(x) log(x^n + b^n) / (x^n + b^n), lower=0, upper=Inf)
# Closed formula;
# Int = integrate( \(x) (x^3 + 2*x^2 + 3*x + 4) / (x^4+x^3+x^2+x+1), lower=0, upper=1)
Int = - (digamma(1/n) + Euler + log(n));
Const = - pi/sin(pi/n) / b^(n-1) * (pi*cos(pi/n)/sin(pi/n) / n - log(b));
(log(n) + Int) * pi/sin(pi/n) / n / b^(n-1) + Const;
- (digamma(1/n) + Euler) * pi/sin(pi/n) / n / b^(n-1) + Const;


# Derivation:
# dI/da:
integrate(\(x) n*a^(n-1)/(b^n - a^n)* (1/(x^n + a^n) -  1/(x^n + b^n)), lower=0, upper=Inf)
n*a^(n-1)/(b^n - a^n) * (1/a^(n-1) - 1/b^(n-1)) * pi/sin(pi/n) / n;
(1 - (a/b)^(n-1)) / (1 - (a/b)^n) / b^n * pi/sin(pi/n);
# Back-integration: I( dI(a) da );
# a = 0 =>
Const = - pi/sin(pi/n) / b^(n-1) * (pi*cos(pi/n)/sin(pi/n) / n - log(b));


# Helper
integrate(function(x) log(x)/(x^5 + 1), 0, Inf)
- pi^2*cos(pi/5)/sin(pi/5)^2 / 25
#
p = sqrt(3)
integrate(function(x) x^p * log(x) / (x^5 + 1), 0, Inf)
- pi^2*cos(pi*(p+1)/5)/sin(pi*(p+1)/5)^2 / 25

#
b = 3^(1/4);
integrate(function(x) log(x)/(x^5 + b^5), 0, Inf)
- pi^2*cos(pi/5)/sin(pi/5)^2 / 5^2 / b^4 +
	+ pi*log(b)/sin(pi/5) / 5 / b^4;


# Derivation:
m = cos(pi/5) + 1i*sin(pi/5)
v = pi/sin(pi/5)/5
integrate(function(x) log(x)/(x^5 + 1), 0, Inf)
2i*pi*(log(m)/((m+1)*prod(m - m^c(3,7,9))) + (1/5)*v*exp(2i*pi/5)) / (1 - exp(2i*pi/5))
2i/25*pi^2*(1i/m^4 + exp(2i*pi/5)/sin(pi/5)) / (1 - exp(2i*pi/5))
2i/25*pi^2*(1i*exp(-4i*pi/5) + exp(2i*pi/5)/sin(pi/5)) / (1 - exp(2i*pi/5))
- pi^2*(exp(1i*pi/5)/sin(pi/5) - 1i) / sin(pi/5) / 25



### Generalization
### log(x) / (x^n + 1)
n = 7
integrate(function(x) log(x)/(x^n + 1), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2


#
n = 7; b = 3^(1/4);
integrate(function(x) log(x)/(x^n + b^n), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2 / b^(n-1) +
	+ pi*log(b)/sin(pi/n) / n / b^(n-1);


###  log(x) * x^p / (x^n + 1)
n = 7
p = sqrt(2)
integrate(function(x) x^p * log(x)/(x^n + 1), 0, Inf)
- pi^2*cos(pi*(p+1)/n)/sin(pi*(p+1)/n)^2 / n^2

###
n = 7
p = sqrt(2)
b = 3^(1/4)
integrate(function(x) x^p * log(x)/(x^n + b^n), 0, Inf)
- pi^2*cos(pi*(p+1)/n)/sin(pi*(p+1)/n)^2 / n^2 * b^(p + 1 - n) +
	+ pi*log(b)/(n*sin(pi*(p+1)/n)) * b^(p + 1 - n)


### I( x^p * log(x^n + b^n) / (x^n + b^n) )
n = sqrt(11)
b = sqrt(5);
integrate(function(x) log(x^n + b^n) / (x^n + b^n), lower=0, upper=Inf)
Const = - pi/sin(pi/n) / b^(n-1) * (pi*cos(pi/n)/sin(pi/n) / n - log(b));
- (digamma(1/n) + Euler) * pi/sin(pi/n) / n / b^(n-1) + Const;

### Full:
p = sqrt(3); n = sqrt(11);
integrate(function(x) x^p * log(x^n + b^n) / (x^n + b^n), lower=0, upper=Inf)
Const = - pi/sin(pi*(p+1)/n) / b^(n - (p+1)) *
	(pi*cos(pi*(p+1)/n)/sin(pi*(p+1)/n) / n - log(b));
- (digamma((p+1)/n) + Euler) * pi/sin(pi*(p+1)/n) / n / b^(n - (p+1)) + Const;


##################
### Transforms ###
##################

###
Catalan = 0.915965594177219015054603514;

###
integrate(function(x) log(x+1) / (x^2 + 1), lower=-1, upper=1)
pi*log(2)/4 - Catalan;

# Base:
integrate(function(x) log(x + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)/4 + Catalan;
# x => (1-x)/(1+x)
integrate(function(x) (log(2) - log(x+1)) / (x^2 + 1), lower=-1, upper=1)
# =>
pi*log(2)/2 - integrate(function(x) log(x+1) / (x^2 + 1), lower=-1, upper=1)$value


### x => (a-x)/(a+x)
a = 7/5
integrate(function(x) log(x+a) / (x^2 + a^2), lower=-a, upper=a)
pi*log(a)/(2*a) + pi*log(2)/(4*a) - Catalan/a;


### x => (1-x)/(a+x)
a = 7/5
integrate(function(x) log(x+a) / (2*x^2 + 2*(a-1)*x + a^2 + 1), lower=-a, upper=1)
pi*log(a+1)/(2*(a+1)) - pi*log(2)/(4*(a+1)) - Catalan/(a+1);


#######################

### Composite Fractions

# qncubed3: Complex Analysis: An Integral from @MichaelPennMath
# https://www.youtube.com/watch?v=LH4i9XJsz_I
# - Generalization of the 2nd (sub-) Integral;

p = 1/5
k = 3
integrate(function(x) log(x) / x^p / (x+k)^2, lower=0, upper=Inf)
pi*((p*log(k) - 1)*sin(pi*p) + pi*p*cos(pi*p)) / k^(p+1) / sin(pi*p)^2


# Residue: x = - k;
# x^(p-1)*(1 - p*log(x))
2i*pi*k^(-p-1)*exp(-1i*pi*(p+1))*(1 - p*log(k) - 1i*pi*p)

# Div:
(1 - exp(-2i*pi*p))

# Helper
integrate(function(x) x^(-p)/(x+k)^2, lower=0, upper=Inf)
pi*p/sin(pi*p)/k^(p+1)


###
p = sqrt(1/5)
k = sqrt(7)
integrate(function(x) log(x)/x^p/(x+k)^2, lower=0, upper=Inf)
full = pi*((p*log(k) - 1)*sin(pi*p) + pi*p*cos(pi*p)) / k^(p+1) / sin(pi*p)^2
print(full)
#
- full/2 + integrate(function(x) log(x)/x^p/(x+k)^2, lower=0, upper=1)$value
- full/2 + integrate(function(x) log(x)/x^p/(x+k)^2, lower=1, upper=Inf)$value
# TODO: formula to compute 0.3017509;


##########################
##########################

### Composite Log 
### w. Simple Fraction

# qncubed3: Complex Analysis: Logarithms and Branch Cuts
# https://www.youtube.com/watch?v=8fYC8ldo__Y

k = sqrt(3)
p = sqrt(2) - 1
integrate(function(x) log(k*x + 1) / x^(p+1), 0, Inf)
pi*k^p / (p*sin(pi*p))


### Gen 1:
# - continuous power;

### n = 3
p = sqrt(2)
integrate(function(x) log(x^3 + 1) / x^(p+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/3))

# Int by parts:
integrate(function(x) 3/p * x^(2-p) / (x^3+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/3))


### n = 5
p = sqrt(2)
integrate(function(x) log(x^5 + 1)/x^(p+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/5))


### n = 6
p = sqrt(2)
n = 6;
integrate(function(x) log(x^n + 1)/x^(p+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/n))
1/p*pi/sin(pi*p/n)


### n = 7
p = sqrt(2)
n = 7;
integrate(function(x) log(x^n + 1)/x^(p+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/n))
1/p*pi/sin(pi*p/n)


### n = pi
p = sqrt(2)
n = pi;
integrate(function(x) log(x^n + 1)/x^(p+1), 0, Inf)
1/p*pi/sin(pi*p/n)


##############
##############

### on [0, 1]

# Note: formula based on Digamma available!
logcos.sh = function(n, p) {
	# shifted: p => n - p - 1;
	if(p == 0) {
		# r = - (digamma(1) - digamma(1/2))*n/2 + log(2);
		r = - (n-1)*log(2)
		return(r);
	}
	# shift: to match old variant;
	p = n - p - 1;
	# Note: NO division by n!
	r = (digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / 2;
	r - pi/sin((p+1)*pi/n)/2;
}
logcos.old = function(n, p) {
	nint = trunc(n);
	if(n != nint) {
		warning("n and p must be integers! Result will be inaccurate!");
	}
	if(nint %% 2 == 1) {
		pn = seq(2, n, by=2) * pi/n;
	} else {
		pn = seq(1, n, by=2) * pi/n;
	}
	cs = cos(pn*(n-p)); csH = cos(pn/2);
	sign = if((n-p-1) %% 2 == 0) 1 else -1;
	# Note: omitted division by n;
	r = 2*sign*sum(cs*log(csH));
	return(r);
}

### n = ODD or EVEN
n = 7;
p = sqrt(3); # INTEGER between [1, n-2]!
# [low precision for p = n - 2]
integrate(function(x) log(x^n + 1)/x^(p+1), 0, 1)
pracma::integral(function(x) log(x^n + 1)/x^(p+1), 0, 1)
#
pi/(2*p)/sin(pi*(n-p)/n) - log(2)/p + logcos.sh(n, p)/p;
(digamma(((n-p)/n + 1)/2) - digamma((n-p)/n/2)) / (2*p) - log(2)/p;


###
n = 10;
p = 1; # INTEGER between [1, n-2]!
# [low precision for p = n - 2]
integrate(function(x) log(x^n + 1)/x^(p+1), 0, 1)
pracma::integral(function(x) log(x^n + 1)/x^(p+1), 0, 1)
#
pi/(2*p)/sin(pi*(n-p)/n) - log(2)/p + logcos.sh(n, p)/p;


# Derivation:
- log(2)/p + integrate(function(x) n/p * x^(n-p-1)/(x^n + 1), 0, 1)$value
# [high precision even with p = n-2]


###
n = sqrt(11);
p = sqrt(3);
integrate(function(x) log(x^n + 1)/x^(p+1), 0, 1)
pracma::integral(function(x) log(x^n + 1)/x^(p+1), 0, 1)
#
pi/(2*p)/sin(pi*(n-p)/n) - log(2)/p + logcos.sh(n, p)/p;


########################
########################

# qncubed3: Complex Analysis: Fancy Branch Cuts
# https://www.youtube.com/watch?v=2EnE78LKY3Y

integrate(function(x) log(x^4 + 1)/(x^2 + 1), 0, Inf)
pi*log(2 + sqrt(2))

### Gen 1:
integrate(function(x) log(x^(4*2) + 1)/(x^2 + 1), 0, Inf)
2*2*pi*log(2) + 2*pi*log(prod(cos(pi/4 - pi/16), cos(pi/4 - 3*pi/16)))

# Derivation:
pi*log(2)/2 + pi*log((1 + sin(pi/8))*(1 + sin(3*pi/8))/(1 - sin(pi/8))/(1 - sin(3*pi/8)))/2
pi*log(2)/2 + pi*log((1 + sin(pi/8))*(1 + sin(3*pi/8))/cos(pi/8)/cos(3*pi/8))
pi*log(2)/2 + pi*log((1 + sin(pi/8))*(1 + sin(3*pi/8)) * 2^2/sqrt(2))
2*pi*log(2) + pi*log((1 + sin(pi/8))*(1 + sin(3*pi/8)))
# alternative:
2*pi*log(2) + pi*log(1 + 1/sin(pi/8)/2 + sin(2*pi/8)/2)
pi*log(2) + pi*log(2 + 1/sin(pi/8) + sin(2*pi/8))


###
k = 3
integrate(function(x) log(x^(4*k) + 1)/(x^2 + 1), 0, Inf)
2*k*pi*log(2) + 2*pi*log(prod(cos(pi/4 - pi * seq(1, 2*k-1, by=2)/(8*k))))

# Derivation:
k*pi*log(2) + pi*log(prod(1 + sin(pi * seq(1, 2*k-1, by=2)/(4*k))))
k*pi*log(2) + pi*log(prod(1 + cos(pi/2 - pi * seq(1, 2*k-1, by=2)/(4*k))))
2*k*pi*log(2) + 2*pi*log(prod(cos(pi/4 - pi * seq(1, 2*k-1, by=2)/(8*k))))


#################
#################

### Michael Penn: Can you guess the trick for this integral?
# https://www.youtube.com/watch?v=8R0MiRYmjbk

integrate(function(x) log(1 - x) / (x^2 + 1), 0, 1)
pi*log(2)/8 - Catalan

###
integrate(function(x) log(1 + x) / (x^2 + 1), 0, 1)
pi*log(2)/8

###
integrate(function(x) log(1 - x) / (x + 1) / sqrt(x), 0, 1)
log(2)*pi/2 - 2*Catalan


### Gen 1:
k = 3
integrate(function(x) log(k - x)/(x^2 + k^2), 0, k)
(pi*log(2)/8 - Catalan)/k + pi*log(k)/(4*k)

###
k = 1/3
integrate(function(x) log(x + k) / (x^2 + k^2), k, Inf)
pi*log(2)/(8*k) + pi*log(k)/(4*k) + Catalan/k
#
integrate(function(x) log(x + k) / (x^2 + k^2), 0, k)
pi*log(2)/(8*k) + pi*log(k)/(4*k)

###
k = sqrt(5) - sqrt(3)
integrate(function(x) log(k^2 - x) / (x + k^2) / sqrt(x), 0, k^2)
pi*log(2)/(2*k) + pi*log(k)/k - 2*Catalan/k


### Transforms

# x => 1 - 1/x
integrate(function(x) - log(1 - x) / (x^2 + 1), 0, 1)
integrate(function(x) log(x) / (2*x^2 - 2*x + 1), 1, Inf)
- pi*log(2)/8 + Catalan

###
integrate(function(x) log(x) / (2*x^2 + 2*x + 1), 0, Inf)
- pi*log(2)/8


###
# x => b - (b-a)/(x+1)
a = sqrt(2); b = sqrt(3);
p = 1/2
integrate(function(x) log(x) * ((x - a)*(b - x))^p, a, b)
integrate(function(x) log(b - (b-a)/(x+1)) * x^p / (x + 1)^(2*p + 2) * (b-a)^(2*p+1), 0, Inf)

### p = 1
a = sqrt(2); b = sqrt(3);
integrate(function(x) log(b*x + a) * x / (x + 1)^4, 0, Inf)
d = (b-a)^3;
b^2*(b - 3*a)*log(b)/(6*d) - a^2*(a - 3*b)*log(a)/(6*d) + 5/36 +
	+ (4*b^3 - 9*(a+b)*b^2 + 36*a*b^2)/(36*d) - (4*a^3 - 9*(a+b)*a^2 + 36*a^2*b)/(36*d)

# Derivation:
integrate(function(x) log(b - (b-a)/(x+1)) * x / (x + 1)^4, 0, Inf)
integrate(function(x) log((b*x + a)/(x+1)) * x / (x + 1)^4, 0, Inf)
integrate(function(x) - log(x) * (x^2 - (a+b)*x + a*b) / (b-a)^3, a, b)
#
integrate(function(x) log(x+1) * x / (x + 1)^4, 0, Inf)
5/36

### I( log(b*x + a) / (x + 1)^3 )
a = sqrt(2); b = sqrt(3);
d = (b-a)^3; # actually d3;
integrate(function(x) log(b*x + a) / (x + 1)^3, 0, Inf)
b^2*(b - 3*a)*log(b)/(6*d) - a^2*(a - 3*b)*log(a)/(6*d) + log(a)/3 +
	+ log(a/b) / (a/b-1)^3 / 3 +
	+ 1/3 * (1-a/b)/(a/b-1)^3 + 1/3 * a*b/(b-a)^2 + 1/6/(a/b-1);

# integrate(function(x) b/3/(b*x + a) / (x + 1)^3, 0, Inf)
# 1/3/(a/b-1)^3 * integrate(function(x) (x^2 + (3-a/b)*x + 3-3*a/b+a^2/b^2) / (x + 1)^3 +
#	- 1/(x + a/b), 0, Inf)$value
# - [above] more simplifications may be possible;
# - may have been simpler to integrate by parts directly;


### log(x+1) * x / (x + 1)^n
integrate(function(x) log(x+1) * x / (x + 1)^4, 0, Inf)
5/36

###
integrate(function(x) log(x+1) * x / (x + 1)^5, 0, Inf)
35 / factorial(6)

###
integrate(function(x) log(x+1) * x / (x + 1)^6, 0, Inf)
81/5 / factorial(6)


###########
### p = 1/2
a = sqrt(2); b = sqrt(3);
p = 1/2; d = (b - a);
integrate(function(x) log(x) * ((x - a)*(b - x))^p, a, b)
# TODO: ???

### [old]
f = function(x) {
	tmp = (x - (a+b)/2)*2/d;
	tmp = round(tmp, 12); # keep between [-1, 1];
	(d/2)^(2*p + 1)*(tmp*sqrt(1-tmp^2) + asin(tmp))/2;
}
log(b)*f(b) - log(a)*f(a) - integrate(function(x) f(x)/x, a, b)$value
pi/4*(d/2)^(2*p + 1)*log(a*b) - integrate(function(x) f(x)/x, a, b)$value
# TODO


# Helper:
p = 1/2
d = (b - a);
integrate(function(x) ((x - a)*(b - x))^p, a, b)
integrate(function(x) ((x + d/2)*(d/2 - x))^p, -d/2, d/2)
integrate(function(x) (d/2)^(2*p + 1)*((1 + x)*(1 - x))^p, -1, 1)
(d/2)^(2*p + 1) * pi/2

# Derivation:
(d/2)^(2*p + 1)*(sin(2*asin((x - (a+b)/2)*2/d))/2 + asin((x - (a+b)/2)*2/d))/2


####################
####################

### log(x) / sqrt(1 - x)
### on [0, 1]

integrate(\(x) x*log(x) / sqrt(1 - x^2), 0, 1)
log(2) - 1

###
integrate(\(x) log(x) / sqrt(1 - x), 0, 1)
4*log(2) - 4

### Subst: x => sin(x)^2
# or simpler: Integration by parts;
integrate(\(x) sin(x) * log(sin(x)), 0, pi/2)
log(2) - 1

###
integrate(\(x) sin(2*x) * log(sin(x)), 0, pi/2)
-1/2

###
integrate(\(x) sin(3*x) * log(sin(x)), 0, pi/2, rel.tol=1E-8)
log(2)/3 - 7/9


########################
########################

### I( (x^p - 1) / ((x^n + 1) * log(x)) )
# Maths 505: AN ABSOLUTE BEAST!!! And we're solving it using Feynman's technique
# https://www.youtube.com/watch?v=wnRZmd1vaxQ

###
integrate(\(x) (x - 1) / (log(x) * (x^2 + 1)), 0, 1)
log(gamma(1/4)^2 / gamma(1/2)^3) - log(2)/2

### Gen 1:
n = sqrt(5)
integrate(\(x) (x - 1) / (log(x) * (x^n + 1)), 0, 1)
log(gamma(1/(2*n))^2 * gamma(2/n) / gamma(1/n)^3) - log(2)/n

### Gen 2:
n = sqrt(5)
p = sqrt(3)
integrate(\(x) (x^p - 1) / (log(x) * (x^n + 1)), 0, 1)
log(gamma(1/(2*n))^2 * gamma((p+1)/n) / gamma((p+1)/(2*n))^2 / gamma(1/n)) - p*log(2)/n

### Gen 3: I( x^p1 * (x^p2 - 1) / (log(x) * (x^n + 1)) )
n = sqrt(5)
p1 = sqrt(3); p2 = sqrt(2)
integrate(\(x) x^p1 * (x^p2 - 1) / (log(x) * (x^n + 1)), 0, 1);
log(gamma((p1+p2+1)/n) / gamma((p1+p2+1)/(2*n))^2) +
	- log(gamma((p1+1)/n) / gamma((p1+1)/(2*n))^2) +
	- (p1+p2)*log(2)/n + p1*log(2)/n;


### Special Case:
p = 2; n = 1;
integrate(\(x) (x - 1) / log(x), 0, 1)
log(2)


### Gen 3: Pow = 2
n = sqrt(5)
p1 = sqrt(3); p2 = sqrt(2)
integrate(\(x) (x^p1 - 1) * (x^p2 - 1) / (log(x)^2 * (x^n + 1)), 0, 1);
# TODO:
# log(gamma((p1+p2+1)/n) / gamma((p1+p2+1)/(2*n))^2) +
#	- log(gamma((p1+1)/n) / gamma((p1+1)/(2*n))^2) +
	- p1*p2*log(2)/n;


########################
########################

#################
### Diff Type ###
#################

### I( x^p * log(x) / (x^n - 1) )
n = 7; p = sqrt(5)
tol = 1E-7;
integrate(\(x) x^p * log(x) / (x^n - 1), 0, 1 - tol)$value +
	+ integrate(\(x) x^p * log(x) / (x^n - 1), 1 + tol, Inf)$value
pi^2 / sin(pi*(p+1)/n)^2 / n^2;


### I( x^p * log(x)^2 / (x^n - 1) )
n = 7; p = sqrt(5)
tol = 1E-7;
integrate(\(x) x^p * log(x)^2 / (x^n - 1), 0, 1 - tol)$value +
	+ integrate(\(x) x^p * log(x)^2 / (x^n - 1), 1 + tol, Inf)$value
- 2*pi^3 * cos(pi*(p+1)/n) / sin(pi*(p+1)/n)^3 / n^3;


### Base:
n = 7; p = sqrt(5)
tol = 1E-7;
integrate(\(x) x^p / (x^n - 1), 0, 1 - tol)$value +
	+ integrate(\(x) x^p / (x^n - 1), 1 + tol, Inf)$value
- pi / tan(pi*(p+1)/n) / n;

