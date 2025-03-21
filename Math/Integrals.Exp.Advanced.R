

### Examples:
# I( atan(x) / (exp(k*x) - 1) )
# I( atan(x) / (exp(k*x) + 1) )
# I( atan(x) / sinh(k*x) )
# I( atan(x^2) / cosh(k*x) )
# I( atan(x^4) / cosh(k*x) )
# I( atan(x^6) / cosh(k*x) )
# I( atan(x^2 + 1) / cosh(k*x) )
# I( atan(exp(-x)) / x )


####################

### Helper Functions

Euler   = 0.57721566490153286060651209008240243079;
constEuler = Euler;


#######################
#######################

### I( atan(x) / (exp(k*x) - 1) )
# Flammable Maths: A ""-relaxing Integral Experience
# https://www.youtube.com/watch?v=QVgfL8Le0I0
# Flammable Maths:: One Spicy Class of Integrals.
# https://www.youtube.com/watch?v=HDaQAHhwtLo


###
integrate(\(x) atan(x) / (exp(2*pi*x) - 1), 0, Inf)
1/2 - log(2*pi)/4

###
integrate(\(x) atan(x) / (exp(x) - 1), 0, Inf)
pi*log(gamma(1/(2*pi))) - (pi - 1/2)*log(2*pi) + 1/2

### I( atan(x) / (exp(n*x) - 1) )
n = sqrt(3)
integrate(\(x) atan(x) / (exp(n*x) - 1), 0, Inf)
log(gamma(n/(2*pi))) * pi/n + (1 - pi/n)/2*log(2*pi/n) - log(2*pi) * pi/(2*n) + 1/2

### I( atan(x) / (exp(n*x) + 1) )
n = sqrt(3)
integrate(\(x) atan(x) / (exp(n*x) + 1), 0, Inf)
log( gamma(n/(2*pi)) / gamma(n/pi) ) * pi/n +
	- 1/2*log(pi/n) + (1 - pi/n)/2*log(2) - 1/2;

### I( atan(x) / sinh(k*x) )
k = sqrt(3)
integrate(\(x) atan(x) / sinh(k*x), 0, Inf)
log( gamma(k/(2*pi))^2 / gamma(k/pi) / pi ) * pi/k +
	+ pi/(2*k) * log(k/2) + (1 - pi/k)*log(2);


### Transformations
integrate(\(x) - atan(x) / (exp(x) - 1), 0, Inf)
#
integrate(\(x) log(1 - exp(-x)) / (x^2 + 1), 0, Inf)
integrate(\(x) log(1 - exp(-1/x)) / (x^2 + 1), 0, Inf)

### I( log((exp(n*x) - 1)/(exp(n*x) + 1)) / (x^2 + 1) )
n = sqrt(5) - sqrt(3); # numerical problems for n > 1.5
integrate(\(x) - log((exp(n*x) - 1)/(exp(n*x) + 1)) / (x^2 + 1), 0, Inf)
log( gamma(n/(2*pi))^2 / gamma(n/pi) ) * pi +
	- pi*log(pi/n) - pi/2*log(n) + (n - 3/2*pi) * log(2);

# - see also: I( x^p * log((exp(x) - 1) / (exp(x) + 1)) )
#   in file: Integrals.Exp.R;


### Perverse Transformations
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - x/(x^2+2), 0, Inf)
log(2)/2 - pi*log(gamma(1/(2*pi))) + (pi - 1/2)*log(2*pi) - 1/2

### including also Euler's constant:
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - exp(-1/x)/x, 0, Inf)
- pi*log(gamma(1/(2*pi))) + (pi - 1/2)*log(2*pi) - 1/2 + constEuler
# numerically problematic:
integrate(\(x) log(exp(x) - 1) / (x^2 + 1) - exp(-x)/x, 0, Inf)

# library(Rmpfr)
ff = \(x) {
	x = mpfr(x, 240);
	y = log(exp(x) - 1) / (x^2 + 1) - exp(-x)/x;
	as.numeric(y);
}
# very crude:
integrate(ff, 1/500, 1)$value + integrate(ff, 1, 500)$value


###
# x => log(y)
integrate(\(x) atan(x) / (exp(x) - 1), 0, Inf)
integrate(\(x) atan(log(x)) / (x*(x - 1)), 1, Inf)
integrate(\(x) atan(log(x)) / (x - 1), 0, 1)
integrate(\(x) - log(1-x) / (x * (log(x)^2 + 1)), 0, 1)
integrate(\(x) - log(1 - exp(-x)) / (x^2 + 1), 0, Inf) # redundant
pi*log(gamma(1/(2*pi))) - (pi - 1/2)*log(2*pi) + 1/2

### I( atan(x) / (exp(x) + 1) )
integrate(\(x) atan(x) / (exp(2*x) - 1), 0, Inf)
log(gamma(1/pi)) * pi/2 + (1 - pi/2)/2*log(pi) - log(2*pi) * pi/4 + 1/2
# =>
integrate(\(x) atan(x) / (exp(x) + 1), 0, Inf)
integrate(\(x) atan(log(x)) / (x*(x + 1)), 1, Inf)
integrate(\(x) - atan(log(x)) / (x + 1), 0, 1)
integrate(\(x) log(x + 1) / (x*(log(x)^2 + 1)), 0, 1)
pi*log( gamma(1/(2*pi)) / gamma(1/pi) ) +
	- 1/2*log(pi) - pi/2*log(2) + 1/2*log(2) - 1/2;

# TODO:
# - find any utility?


#####################

### ATAN(x^2)

### I( atan(x^2) / cosh(k*x) )
# - for Base-formula, see file:
#   Integrals.Exp.Log.R;
k = sqrt(3)
integrate(\(x) atan(x^2) / cosh(k*x), 0, Inf)
b = exp(c(1i,-1i) * pi / 4)
pi^2/(4*k) - diff( log(pracma::gammaz(k*b/(2*pi) + 3/4) /
	pracma::gammaz(k*b/(2*pi) + 1/4)) ) * pi/k * 1i;


### I( atan(x^4) / cosh(k*x) )
k = sqrt(3)
integrate(\(x) atan(x^4) / cosh(k*x), 0, Inf)
b = exp(c(1i,-1i,3i,-3i) * pi / 8)
pi^2/(4*k) - sum( c(1,-1,-1,1) * log(pracma::gammaz(k*b/(2*pi) + 3/4) /
	pracma::gammaz(k*b/(2*pi) + 1/4)) ) * pi/k * 1i;


### I( atan(x^6) / cosh(k*x) )
k = sqrt(3)
integrate(\(x) atan(x^6) / cosh(k*x), 0, Inf)
b  = exp(c(1,3,5,-5,-3,-1) * 1i * pi / 12);
sg = c(-1,1,-1,1,-1,1);
pi^2/(4*k) - sum( sg * log(pracma::gammaz(k*b/(2*pi) + 3/4) /
	pracma::gammaz(k*b/(2*pi) + 1/4)) ) * pi/k * 1i;


### Exercises:

# Tough Exercise for students:
#   "Prove the following results,
#    where k is a parameter:"
k = sqrt(3) - sqrt(2);
integrate(\(x) atan(x^2 + 1) / cosh(k*x), 0, Inf)
b = sqrt(1 + c(1i,-1i));
pi^2/(4*k) - diff( log(pracma::gammaz(k*b/(2*pi) + 3/4) /
	pracma::gammaz(k*b/(2*pi) + 1/4)) ) * pi/k * 1i;


#####################
#####################

### I( atan(x/k) / (exp(n*x) - 1) )
n = sqrt(3); k = sqrt(5);
integrate(\(x) atan(x/k) / (exp(n*x) - 1), 0, Inf)
log(gamma(k*n/(2*pi))) * pi/n + (k - pi/n)/2*log(2*pi/(k*n)) - log(2*pi) * pi/(2*n) + k/2;

### Derived:

### I( x / (x^2 + b^2) * ... )

### I( x / (x^2 + k^2) * 1/(exp(n*x) - 1) )
n = sqrt(3); k = sqrt(5);
integrate(\(x) x / (x^2 + k^2) *  1/(exp(n*x) - 1), 0, Inf)
- 1/2 * digamma(k*n/(2*pi)) +
	- 1/2*log(2*pi/(k*n)) - pi/(2*k*n);


### I( x / (x^2 + k^2)^2 / (exp(n*x) - 1) )
n = 1/sqrt(3); k = sqrt(5);
integrate(\(x) x / (x^2+k^2)^2 / (exp(n*x) - 1), 0, Inf)
pracma::psi(1, k*n/(2*pi)) * n / (8*pi*k) +
	- 1/(4*k^2) - pi/(4*k^3*n);


### Series: exp(n*x) + 1

### I( x / (x^2 + k^2) * 1/(exp(n*x) + 1) )
n = sqrt(3); k = 1 / sqrt(5);
integrate(\(x) x / (x^2 + k^2) *  1/(exp(n*x) + 1), 0, Inf)
digamma(k*n/pi) - 1/2 * digamma(k*n/(2*pi)) +
	+ 1/2*log(pi/(2*k*n));


### I( x / (x^2 + k^2)^2 / (exp(n*x) + 1) )
n = 1/sqrt(3); k = 1/sqrt(5);
integrate(\(x) x / (x^2+k^2)^2 / (exp(n*x) + 1), 0, Inf)
(pracma::psi(1, k*n/(2*pi)) - 4*pracma::psi(1, k*n/pi)) * n / (8*pi*k) + 1/(4*k^2);


### Series: sinh()

### I( x / (x^2 + b^2) * 1/sinh(n*x) )
n = sqrt(3); b = 1 / sqrt(5);
integrate(\(x) x / (x^2 + b^2) *  1/sinh(n*x), 0, Inf)
digamma(b*n/pi) - digamma(b*n/(2*pi)) - log(2) - pi/(2*b*n);


### I( 1 / (x^2 + b^2) * ... )

### I( 1 / ((x^2 + b^2) * cosh(x)) ) on [0, Inf]
# - see Integrals.Exp.Log.R;
b = sqrt(5);
integrate(\(x) 1 / ((x^2 + b^2) * cosh(x)), 0, Inf)
(digamma(b/(2*pi) + 3/4) - digamma(b/(2*pi) + 1/4)) / (2*b)

### Gen: I( 1 / ((x^2 + b^2) * cosh(n*x)) )
b = sqrt(5); n = sqrt(3);
integrate(\(x) 1 / ((x^2 + b^2) * cosh(n*x)), 0, Inf)
(digamma(b*n/(2*pi) + 3/4) - digamma(b*n/(2*pi) + 1/4)) / (2*b)

###
integrate(\(x) 1 / ((x^2 + b^2) * sinh(n*x)) - exp(-x)/x/(n*b^2), 0, Inf)
# TODO
# Note: + Euler;


### Other:

### I( atan(k/x) / (exp(n*x) - 1) )
n = sqrt(3); k = sqrt(5);
integrate(\(x) (atan(k/x) - atan(1/x)) / (exp(n*x) - 1), 0, Inf)
k = c(1, k)
diff( -1/2 * log(gamma(k*n/(2*pi))) * 2*pi/n +
	- k/2*log(2*pi/n) + 1/2*(k*log(k) - k) - pi/(2*n)*log(k) );
# TODO

# Note:
# - redundant: atan(k/x) = pi/2 - atan(x/k);
# - however, still useful as proof of concept of back-integration;


###################
###################

### I( atan(exp(-x)) / x ) on [0, Inf]
integrate(\(x) atan(exp(-x))/x - pi/4*exp(-x)/x, 0, Inf)
pi/4 * log(4*pi^3) - pi*log(gamma(1/4)) + pi*Euler/4

###
pracma::integral(\(x) log(x) / cosh(x), 0, Inf)
pi/2 * log(4*pi^3) - 2*pi*log(gamma(1/4))


####################
####################

### I( exp(cos(x)) * cos(sin(x)) / (x^2 + 1) )
# Maths 505: MESMERIZING integral result! int 0 to infty e^(cos(x))cos(sin(x))/(1+x^2)
# https://www.youtube.com/watch?v=hxmCDKuZLsU

Int = function(FUN, lower, upper, subdiv = 1024, rel.tol=1E-6) {
	integrate(FUN, lower, upper,
		subdivisions = subdiv, rel.tol=rel.tol);
}
Int(\(x) exp(cos(x)) * cos(sin(x)) / (x^2 + 1), 0, Inf, 4400)
pi/2 * exp(exp(-1))

