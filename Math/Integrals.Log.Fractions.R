########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### log-Fractions
###
### draft v.0.1a



##################
### Logarithms ###
##################

#################
### Fractions ###
#################

### Feynman trick & Other tricks
# - usually (much) simpler than Contour integration;

### Basic / Helper
b = sqrt(3)
integrate(function(x) log(x) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(b)/(2*b)
#
integrate(function(x) log(x) / (b^2*x^2 + 1), lower=0, upper=Inf)
- pi*log(b)/(2*b)

### TODO:
# - using Feynman's trick;
# - Note: Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)
Catalan = 0.915965594177219015054603514;
integrate(function(x) log(x + 1) / (x^2 + 1), lower=0, upper=Inf)
pi*log(2)/4 + Catalan;

### TODO:
Catalan = 0.915965594177219015054603514;
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

#############
### log(P[2])
# - using Feynman's trick;
b = sqrt(3); a = sqrt(2)
integrate(function(x) log(x^2 + a^2) / (x^2 + b^2), lower=0, upper=Inf)
pi*log(a + b)/b

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
n = 3
a = sqrt(3); b = sqrt(5)
integrate(function(x) log(x^n + a^n) / (x^n + b^n), lower=0, upper=Inf)

###############
### Case: n = 3
n = 3
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
# dI:
integrate(function(x) n*a^(n-1)/(b^n-a^n)* (1/(x^n + a^n) -  1/(x^n + b^n)), lower=0, upper=Inf)
# [more complicated] x^2 - (m+m^4)*a*x + a^2 =>
# TODO:
sqrt(n)*(1 - (a/b)^(n-1))/(b^n - a^n) * (pi/2 + ...)

# Helper
integrate(function(x) log(x)/(x^5 + 1), 0, Inf)
- pi^2*cos(pi/5)/sin(pi/5)^2 / 25

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
n = 7
integrate(function(x) log(x)/(x^n + 1), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2


#
n = 7; b = 3^(1/4);
integrate(function(x) log(x)/(x^n + b^n), 0, Inf)
- pi^2*cos(pi/n)/sin(pi/n)^2 / n^2 / b^(n-1) +
	+ pi*log(b)/sin(pi/n) / n / b^(n-1);


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

