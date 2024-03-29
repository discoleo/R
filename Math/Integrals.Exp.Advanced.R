
#######################
#######################

### I( atan(x) / (exp(x) - 1) )
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
#####################

### I( atan(x/k) / (exp(n*x) - 1) )
n = sqrt(3); k = sqrt(5);
integrate(\(x) atan(x/k) / (exp(n*x) - 1), 0, Inf)
log(gamma(k*n/(2*pi))) * pi/n + (k - pi/n)/2*log(2*pi/(k*n)) - log(2*pi) * pi/(2*n) + k/2;


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

