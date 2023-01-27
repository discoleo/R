########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### log-Fractions
###
### draft v.0.1g


##################
### Logarithms ###
##################

# - definite integrals;
# - various types of Logarithms combined with fractions;


####################

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

###
n = 7
p = sqrt(2)
integrate(function(x) x^p*log(x)/(x^n + 1), 0, Inf)
- pi^2*cos(pi*(p+1)/n)/sin(pi*(p+1)/n)^2 / n^2

###
n = 7
p = sqrt(2)
b = 3^(1/4)
integrate(function(x) x^p*log(x)/(x^n + b^n), 0, Inf)
- pi^2*cos(pi*(p+1)/n)/sin(pi*(p+1)/n)^2 / n^2 * b^(p + 1 - n) +
	+ pi*log(b)/(n*sin(pi*(p+1)/n)) * b^(p + 1 - n)


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
#######################

### Simple Fractions

# qncubed3: Complex Analysis: Integral of log(x)/(x+1)^2
# https://www.youtube.com/watch?v=tPveHNdBWR8

# Note:
# - easy pole with higher multiplicity;

integrate(function(x) log(x)/(x+1)^2, 0, Inf)
# == 0

###
integrate(function(x) log(x)/(x+1)^3, 0, Inf)
-1/2

###
integrate(function(x) log(x)/(x+1)^4, 0, Inf)
-1/2

### n = 5
integrate(function(x) log(x)/(x+1)^5, 0, Inf)
- 11 / gamma(5)

### n = 6
integrate(function(x) log(x)/(x+1)^6, 0, Inf)
- 50 / gamma(6)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(- 100/x^5 + 48*log(x)/x^5)/gamma(6)


### n = 7
integrate(function(x) log(x)/(x+1)^7, 0, Inf)
- 274 / gamma(7)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(2*274 - 240*log(x)/x^6)/gamma(7)
# - 2i*I + 2*pi/(n-1) = 1i*(2*274 - 240*log(x)/x^6)/gamma(7)


#######################

### Composite Fractions

# qncubed3: Complex Analysis: An Integral from @MichaelPennMath
# https://www.youtube.com/watch?v=LH4i9XJsz_I
# - Generalization of the 2nd (sub-) Integral;

p = 1/5
k = 3
integrate(function(x) log(x)/x^p/(x+k)^2, lower=0, upper=Inf)
pi*((p*log(k) - 1)*sin(pi*p) + pi*p*cos(pi*p)) / k^(p+1) / sin(pi*p)^2


# Residue: x = - k;
# x^(p-1)*(1 - p*log(x))
2i*pi*k^(-p-1)*exp(-1i*pi*(p+1))*(1 - p*log(k) - 1i*pi*p)

# Div:
(1 - exp(-2i*pi*p))

# Helper
integrate(function(x) x^(-p)/(x+k)^2, lower=0, upper=Inf)
pi*p/sin(pi*p)/k^(p+1)


#####################
#####################

#####################
### Log-Fractions ###
### Higher Powers ###
#####################


# qncubed3: Complex Analysis: Is this solution way too overkill?
# https://www.youtube.com/watch?v=vfa5wGE7fRI

###
integrate(function(x) log(x)^2/(x^2 + 1), lower=0, upper=Inf)
pi^3/8

### Gen 1:
b = sqrt(5)
integrate(function(x) log(x)^2/(x^2 + b^2), lower=0, upper=Inf)
pi^3/(2^3*b) + pi/2 * log(b)^2/b


### Helper

# Derivation:
source("Polynomials.Helper.R")

# Generate Auxiliary Polynomial:
solve.polyInt = function(n, verbose=TRUE, subst=BULL) {
	px = toPoly.pm(paste0("x^", seq(n-1), "*b", seq(n-1), collapse="+"));
	px = sum.pm(px, toPoly.pm("x^n"));
	pd = diff.pm(replace.pm(px, toPoly.pm("x+k"), "x"), px);
	b = paste0("b", seq(n-1));
	r = list();
	for(i in seq(n-1, 1)) {
		tmp = pd[pd$x == (i-1), , drop=FALSE];
		tmp$x = NULL;
		tmp = drop.pm(tmp);
		tmp = simplify.pm.pow(tmp, do.gcd=TRUE);
		if(verbose) print.pm(tmp);
		if(nrow(tmp) == 1) {
			r = c(r, list(0, 1));
			pd = replace.pm(pd, 0, b[i]);
			px = replace.pm(px, 0, b[i]);
		} else {
			tmp = solve.pm(tmp, pd, b[i], verbose=FALSE);
			r = c(r, list(tmp$x0, tmp$div));
			div = tmp$div;
			if(is.null(div)) div = data.frame(coeff=1);
			pd = replace.fr.pm(pd, tmp$x0, div, b[i], verbose=FALSE);
			px = replace.fr.pm(px, tmp$x0, div, b[i], verbose=FALSE);
		}
	}
	if(px$coeff[px$x == n] < 0) px$coeff = -px$coeff;
	if( ! is.null(subst)) {
		px = replace.pm(px, subst, "k");
	}
	px = reduce.coef.pm(px);
	return(list(poly=px, Coeff=r));
}


### Gen: Power of log

### "Miss Piggy" Polynomials
# n = 3
x^4 - 4i*pi*x^3 - 4*pi^2*x^2
# n = 4
3*x^5 - 15i*pi*x^4 - 20*pi^2*x^3 - 8*pi^4*x
# n = 5; Factor = 30i*pi;
x^6 - 6i*pi*x^5 - 10*pi^2*x^4 - 8*pi^4*x^2
# n = 6; Factor = 42i*pi;
3*x^7 - 21i*pi*x^6 - 42*pi^2*x^5 - 56*pi^4*x^3 - 32*pi^6*x
# n = 7; Factor = 48i*pi;
3*x^8 - 24i*pi*x^7 - 56*pi^2*x^6 - 112*pi^4*x^4 - 128*pi^6*x^2
# n = 8; Factor = 90i*pi;
5*x^9 - 45i*pi*x^8 - 120*pi^2*x^7 - 336*pi^4*x^5 - 640*pi^6*x^3 - 384*pi^8*x


n = 9
px = solve.polyInt(n + 1, subst=as.pm("2i*pi"))
print.pm(px$poly, "x")


###
b = sqrt(5)
# Note: for b = 1 => 0;
integrate(function(x) log(x)^3/(x^2 + b^2), lower=0, upper=Inf)
pi*log(b)^3 / (2*b) + 3*pi^3*log(b)/(8*b)


###
integrate(function(x) log(x)^4/(x^2 + 1), lower=0, upper=Inf)
5*pi^5/32

#
b = sqrt(5)
integrate(function(x) log(x)^4/(x^2 + b^2), lower=0, upper=Inf)
5*pi^5/(2^5*b) + pi/2 * log(b)^4/b + 6*pi^3*log(b)^2/(8*b)


###
integrate(function(x) log(x)^6/(x^2 + 1), lower=0, upper=Inf)
61*pi^7/128

#
b = sqrt(5)
integrate(function(x) log(x)^6/(x^2 + b^2), lower=0, upper=Inf)
61*pi^7/(2^7*b) + pi/2 * log(b)^6/b + 75*pi^5*log(b)^2/(2^5*b) + 15*pi^3*log(b)^4/(8*b)


###
integrate(function(x) log(x)^8/(x^2 + 1), lower=0, upper=Inf)
1385*pi^9/512

#
b = sqrt(5)
integrate(function(x) log(x)^8/(x^2 + b^2), lower=0, upper=Inf)
1385*pi^9/(2^9*b) + pi/2 * log(b)^8/b + 28*61*pi^7*log(b)^2/(2^7*b) +
	+ 70*5*pi^5*log(b)^4/(2^5*b) + 28*pi^3*log(b)^6/(8*b)


###
integrate(function(x) log(x)^10/(x^2 + 1), lower=0, upper=Inf,
	abs.tol=1E-8, rel.tol=1E-8, subdivisions=4096)
50521*pi^11/2^11


### Gen 2: Both Powers

### Pow = 3
integrate(function(x) log(x)^3/(x^3 + 1), lower=0, upper=Inf)
pi^4*(25/3^4 - 1) / 12


# Aux. Polynomial:
x^4 - 4i*pi*x^3 - 4*pi^2*x^2

# Derivation:
- pi^4*(1 +
	- (1/3^4 - 4/3^3 + 4/9)*m +
	- (5^4/3^4 - 4*5^3/3^3 + 4*5^2/9)*m^5) / 12
- pi^4*(1 - 25*m/3^4 - 25*m^5/3^4) / 12
- pi^4*(1 - 25/3^4) / 12


###########
### Pow = 4
integrate(function(x) log(x)^4/(x^4 + 1), lower=0, upper=Inf)
2*pi^5*(545*sin(pi/4) - 595*sin(5*pi/4)) / 5 / 4^6
# alternative:
2*228*pi^5*sin(pi/4) / 4^6

# Aux. Polynomial:
3*x^5 - 15i*pi*x^4 - 20*pi^2*x^3 - 8*pi^4*x

#
px = toPoly.pm("3*x^5 - 15i*pi*x^4 - 20*pi^2*x^3 - 8*pi^4*x")
# x = c(1,3,5,7)*pi/4;
eval.pm(px, list(x=1i*pi/4, pi=pi)) / pi^5 * 4^5;
# =>
1i*pi^5*(595/exp(3i*pi/4) + 545/exp(9i*pi/4) +
	- 545/exp(15i*pi/4) - 595/exp(21i*pi/4)) / 5 / 4^6;


###########
### Pow = 5
integrate(function(x) log(x)^5/(x^5 + 1), lower=0, upper=Inf)
pi^6*(2*4779*cos(pi/5) - 4*31311*cos(pi/5)^2 + 2*31311 - 3*5^6) / (6*5*5^6)

# Aux. Polynomial:
x^6 - 6i*pi*x^5 - 10*pi^2*x^4 - 8*pi^4*x^2
b5 = - 3*k; b4 = 5/2*k^2; b3 = 0; b2 = -1/2*k^4;

# Derivation:
n = 6
px = solve.polyInt(n, subst=as.pm("2i*pi"))
print.pm(px$poly, "x")


px = toPoly.pm("x^6 - 6i*pi*x^5 - 10*pi^2*x^4 - 8*pi^4*x^2")
eval.pm(px, list(x=1i*pi/5, pi=pi)) / pi^6 * 5^6
# =>
- pi^6*(3*5^6/5 + 4779/(5*m^4) + 31311/(5*m^12) + 31311/(5*m^28) + 4779/(5*m^36)) / (6*5^6)
- pi^6*(3*5^6 - 4779*m - 31311*m^3 + 31311*m^2 + 4779*m^4) / (6*5*5^6)
- pi^6*(3*5^6 - 2*4779*cos(pi/5) + 2*31311*cos(2*pi/5)) / (6*5*5^6)
pi^6*(2*4779*cos(pi/5) - 4*31311*cos(pi/5)^2 + 2*31311 - 3*5^6) / (6*5*5^6)


###########
### Pow = 6
integrate(function(x) log(x)^6/(x^6 + 1), lower=0, upper=Inf)
2*pi^7*(461195*sin(pi/6) + 933849*sin(3*pi/6) - 473935*sin(7*pi/6)) / 7 / 6^8
# alternative:
2*pi^7*(935130*sin(pi/6) + 933849) / 7 / 6^8


# Aux. Polynomial:
3*x^7 - 21i*pi*x^6 - 42*pi^2*x^5 - 56*pi^4*x^3 - 32*pi^6*x

#
px = toPoly.pm("3*x^7 - 21i*pi*x^6 - 42*pi^2*x^5 - 56*pi^4*x^3 - 32*pi^6*x")
# x = c(1,3,5,7,9,11)*pi/6;
eval.pm(px, list(x=1i*pi/6, pi=pi)) / pi^7 * 6^7 / 3;
# =>
1i*pi^7*(473935/exp(5i*pi/6) - 473935/exp(55i*pi/6) +
	+ 933849/exp(15i*pi/6) - 933849/exp(45i*pi/6) +
	+ 461195/exp(25i*pi/6) - 461195/exp(35i*pi/6)) / 7 / 6^8


##########################
##########################

### Composite Log 
### w. Simple Fraction

# qncubed3: Complex Analysis: Logarithms and Branch Cuts
# https://www.youtube.com/watch?v=8fYC8ldo__Y

k = sqrt(3)
p = sqrt(2) - 1
integrate(function(x) log(k*x + 1)/x^(p+1), 0, Inf)
pi*k^p / (p*sin(pi*p))


### n = 3
p = sqrt(2)
integrate(function(x) log(x^3 + 1)/x^(p+1), 0, Inf)
1/p*pi/sin(pi*(1 - p/3))

# Int by parts:
integrate(function(x) 3/p * x^(2-p)/(x^3+1), 0, Inf)
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

