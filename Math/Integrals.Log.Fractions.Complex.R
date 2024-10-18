#########################
##
## Log: Complex Fractions


### Values

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


################
################

### Base:
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)) / gamma(k) / n

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n) + n*log(b)) /
	(gamma(k) * n * b^(n*k - p - 1))

# Special Case: k = 1
p = 1/sqrt(2); n = sqrt(5);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n), 0, Inf)
- gamma((p+1)/n) * gamma(1 - (p+1)/n) *
	(digamma(1 - (p+1)/n) + Euler - n*log(b)) /
	(n * b^(n - p - 1))


###########

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
a = sqrt(3) - sqrt(2); b = sqrt(2)
integrate(\(x) x^p * log(x^n + a^n) / (x^n + b^n)^k, 0, Inf)
# TODO


#################
#################

### Simple:

### I( log(x + 1i)/(x - 1i)) )
integrate(\(x) Re(log(x + 1i)/(x - 1i)), 0, 1)
integrate(\(x) Im(log(x + 1i)/(x - 1i)), 0, 1)
-3/32*pi^2 + log(2)^2/8 + 1i*(5/8*pi*log(2) - Catalan);

# Helper:

# Li2((1−i)/2) # =
5*pi^2/96 - log(2)^2/8 + 1i*(pi/8*log(2) - Catalan);

# Li2((1+i)/2) # =
5*pi^2/96 - log(2)^2/8 - 1i*(pi/8*log(2) - Catalan);


####################

### Generalised:

### I( x^p * atan(x) / (x^2 + 1) )
p = - 1/sqrt(3);
integrate(\(x) x^p * atan(x) / (x^2 + 1), 0, Inf)
pi/4 / sin(pi*p/2) * (digamma(-p/2) - 2*digamma(-p) - Euler);

# Derivation:
pi/4 * gamma((p+1)/2) * gamma(1 - (p+1)/2) +
	- gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) / 4 +
	+ gamma(p+1) * gamma(-p) / 2 *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );
# from below:

# Derivation:
p = - 1/sqrt(3);
integrate(\(x) x^p * (pi/2 - atan(x)) / (x^2 + 1), 0, Inf)
integrate(\(x) x^p * Re(log((x+1i)/(x-1i)) / 2i) / (x^2 + 1), 0, Inf)
gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) / 4 +
	- gamma(p+1) * gamma(-p) / 2 *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );

###
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
(gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) +
	- 1i*(gamma((p+1)/2) * gamma(1 - (p+1)/2) *
	(digamma(1 - (p+1)/2) + Euler))) / 4i +
	- gamma(p+1) * gamma(-p) / 2i *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n) *
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n +
	+ gamma(p+1) * gamma(-p) *
	(digamma(-p) + Euler - log(1i)) / 1i^(- p);


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n)*
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n

# TODO: log(x - 1i) / (x + 1i)


###
n = 4; k = 1; # !!!
p = 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
(gamma((p+1 + 2*k)/n) * gamma(k - (p+1 + 2*k)/n)*
	(digamma(k) - digamma(k - (p+1 + 2*k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n


#############

###
integrate(function(x) Im(log(x^2 + 1i)) / (x^2 + 1), 0, Inf)
pi^2 / 8
integrate(function(x) Re(log(x^2 + 1i)) / (x^2 + 1), 0, Inf)
pi*log(2*cos(pi/8))

# Derivation:
integrate(function(x) atan(x^2)/(x^2 + 1), 0, Inf)
pi^2 / 8
integrate(function(x) Im(log((x^2 + 1i)/(x^2 - 1i))) / (x^2 + 1), 0, Inf)
pi^2 / 4

integrate(function(x) log(x^4 + 1)/(x^2 + 1), 0, Inf)
2*pi*log(2*cos(pi/8))

integrate(function(x) atan(x^4)/(x^2 + 1), 0, Inf)
pi^2 / 8

#################
#################

### Base:
# see file: Integrals.Log.Fractions.R;
b = sqrt(3)
integrate(\(x) log(x^2 + b^2) / (x^4 + 1), 0, Inf)
sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4));


### I( log(x^4 - x^2 + 1) / (x^4 + 1) )
integrate(\(x) log(x^4 - x^2 + 1) / (x^4 + 1), 0, Inf)
b = cos(pi/3) + c(1i,-1i)*sin(pi/3);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );


### Gen: I( log(x^4 + 2*cos(a)*x^2 + 1) / (x^4 + 1) )
a = 1/pi;
integrate(\(x) log(x^4 + 2*cos(a)*x^2 + 1) / (x^4 + 1), 0, Inf)
b = cos(a/2) + c(1i,-1i)*sin(a/2);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );


### I( log(Poly(x^2)) / (x^4 + 1) )
integrate(\(x) log(x^10 + 9*x^8 + 28*x^6 + 35*x^4 + 15*x^2 + 1) / (x^4 + 1), 0, Inf)
cs = 2*cos(2*seq(5)*pi/11);
b  = abs(cs);
sum(sqrt(2)*pi/8 * (2*atan(1/b^2) + log(b^4 + 1)) +
	- pi/2 * (atan(1/b*exp(1i*pi/4))*exp(1i*pi/4) +
		+ atan(1/b*exp(-1i*pi/4))*exp(-1i*pi/4)) );
# Test:
x = - cs^2;
x^5 + 9*x^4 + 28*x^3 + 35*x^2 + 15*x + 1 # == 0

