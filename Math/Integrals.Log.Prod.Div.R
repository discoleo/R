########################
##
## Leonard Mada
## [the one and only]
##
## Integrals: Logarithms
## Log-Products: LOG((1-x)/(1+x))^k
##
## draft v.0.1a


### Derived from:
# LOG((1-x)/(1+x))^k

### Base-Formulas
# in file: Integrals.Log.Prod.R;

################

### Constants

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


#######################

### I( log((1-x)/(1+x))^p )
# Maths 505: This integral is actually one of your favorite constants
# https://www.youtube.com/watch?v=83mUOaF7G9A
# Subst: (1-x)/(1+x) = y;

###
integrate(\(x) log((1-x)/(1+x))^3, 0, 1)
- 12 * (1 - 2^(1-3)) * pracma::zeta(3)

### Gen: I( log((1-x)/(1+x))^p )
p = sqrt(3)
integrate(\(x) abs(log((1-x)/(1+x)))^p, 0, 1)
gamma(p+1) * (2 - 2^(2-p)) * pracma::zeta(p)

### I( x * abs(log(1-x))^p )
p = sqrt(3)
integrate(\(x) x * abs(log(1-x))^p, 0, 1)
gamma(p+1)*(1 - 2^(-p-1))

### I( x^2 * abs(log(1-x))^p )
p = sqrt(3)
integrate(\(x) x^2 * abs(log(1-x))^p, 0, 1)
gamma(p+1) * (1 - 2^(-p) + 3^(-p-1))

### I( x^3 * abs(log(1-x))^p )
p = sqrt(3)
integrate(\(x) x^3 * abs(log(1-x))^p, 0,1)
gamma(p+1) * (1 - 3*2^(-p-1) + 3*3^(-p-1) - 4^(-p-1))


### Derived:

### I( log(1-x) * log(1+x)^2 )
integrate(\(x) log(1-x) * log(1+x)^2, 0, 1)
integrate(\(x) log(1-x^2)^3 / 6, 0, 1)$value +
	- 2 * (1 - 2^(1-3)) * pracma::zeta(3) + 2;
- gamma(1/2) * ((digamma(3/2) + Euler)^3 +
	+ 3*(digamma(3/2) + Euler) * (pracma::psi(1, 1) - pracma::psi(1, 3/2)) +
	+ pracma::psi(2, 3/2) - pracma::psi(2, 1)
	) / gamma(3/2) / 12 - 2 * (1 - 2^(1-3)) * pracma::zeta(3) + 2;

### I( log(1-x)^2 * log(1+x) )
integrate(\(x) log(1-x)^2 * log(1+x), 0, 1)
- gamma(1/2) * ((digamma(3/2) + Euler)^3 +
	+ 3*(digamma(3/2) + Euler) * (pracma::psi(1, 1) - pracma::psi(1, 3/2)) +
	+ pracma::psi(2, 3/2) - pracma::psi(2, 1)
	) / gamma(3/2) / 12 + 3/2 * pracma::zeta(3) +
	- 2/3*log(2)^3 + 2*log(2)^2 - 4*log(2) + 2;


### I( log(1-x) * log(1+x) * (log(1-x) - log(1+x)) )
integrate(\(x) log(1-x) * log(1+x) * (log(1-x) - log(1+x)), 0, 1)
3 * pracma::zeta(3) - 2/3*log(2)^3 + 2*log(2)^2 - 4*log(2);


### Helper:

### Base (D1):
p = 1/3; q = 1/5; n = 1/sqrt(3);
integrate(\(x) x^p * (1 - x^n)^q * log(1-x^n), 0, 1)
gamma((p+1)/n) * gamma(q+1) * (digamma(q+1) - digamma((p+1)/n+q+1)) / gamma((p+1)/n+q+1) / n;

### D2:
integrate(\(x) x^p * (1 - x^n)^q * log(1-x^n)^2, 0, 1)
gamma((p+1)/n) * gamma(q+1) *
	((digamma(q+1) - digamma((p+1)/n+q+1))^2 +
		+ pracma::psi(1, q+1) - pracma::psi(1, (p+1)/n+q+1)) / gamma((p+1)/n+q+1) / n;

### D3:
integrate(\(x) x^p * (1 - x^n)^q * log(1-x^n)^3, 0, 1)
gamma((p+1)/n) * gamma(q+1) *
	((digamma(q+1) - digamma((p+1)/n+q+1))^3 +
	+ 3*(digamma(q+1) - digamma((p+1)/n+q+1)) *
		(pracma::psi(1, q+1) - pracma::psi(1, (p+1)/n+q+1)) +
	+ pracma::psi(2, q+1) - pracma::psi(2, (p+1)/n+q+1)
	) / gamma((p+1)/n+q+1) / n;

# q = 0 =>
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * log(1-x^n)^3, 0, 1)
- gamma((p+1)/n) *
	((digamma((p+1)/n+1) + Euler)^3 +
	+ 3*(digamma((p+1)/n+1) + Euler) *
		(pracma::psi(1, 1) - pracma::psi(1, (p+1)/n+1)) +
	+ pracma::psi(2, (p+1)/n+1) - pracma::psi(2, 1)
	) / gamma((p+1)/n+1) / n;


####################
####################
