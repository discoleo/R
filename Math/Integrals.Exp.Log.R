

### Integrals: Log & Exp


### Helper

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;
gStjelt1 = - 0.0728158454836767248605863758749013191377363383;
# https://en.wikipedia.org/wiki/Stieltjes_constants
dzeta2   = -0.937548254316;
### Glaisher–Kinkelin Constant:
# https://en.wikipedia.org/wiki/Glaisher%E2%80%93Kinkelin_constant
A = exp((log(2*pi) + Euler - 6*dzeta2/pi^2)/12);

# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)

###
dzeta = function(x, dx = 1E-6) {
	(pracma::zeta(x + dx) - pracma::zeta(x)) / dx;
}


#################
#################

### Trig * Log

### Gen: Sin
p = sqrt(3); a = sqrt(5) - sqrt(2); k = sqrt(3);
integrate(\(x) x^p * sin(k*x) * log(x) / exp(a*x), 0, Inf)
gamma(p + 1) * digamma(p + 1) * Im(1 / (a - k*1i)^(p + 1)) +
	- gamma(p + 1) * Im(1 / (a - k*1i)^(p + 1) * log(a - k*1i))


### Gen: Cos
p = sqrt(3); a = sqrt(5) - sqrt(2); k = sqrt(3);
integrate(\(x) x^p * cos(k*x) * log(x) / exp(a*x), 0, Inf)
gamma(p + 1) * digamma(p + 1) * Re(1 / (a - k*1i)^(p + 1)) +
	- gamma(p + 1) * Re(1 / (a - k*1i)^(p + 1) * log(a - k*1i))


###########################
###########################

### I( exp(-x) * log(x)^2 )
# 1.) Maths 505: Another awesome integral with a beautiful result
# https://www.youtube.com/watch?v=wUY6TJxfFTk
# 2.) Dr Peyam: A mathematically stunning formula
# https://www.youtube.com/watch?v=hmYYYhQF1RM
# 3.) Michael Penn: Euler's other constant
# https://www.youtube.com/watch?v=UEqI9GKYozU
# [the basic proof]


###
integrate(\(x) exp(-x) * log(x), 0, Inf)
- Euler

###
n = sqrt(3)
integrate(\(x) exp(-n*x) * log(x), 0, Inf)
- Euler/n - log(n)/n


### Pow = 2
integrate(\(x) exp(-x) * log(x)^2, 0, Inf)
pi^2/6 + Euler^2

### Pow = 3
integrate(\(x) exp(-x) * log(x)^3, 0, Inf)
- 2 * pracma::zeta(3) - Euler*(pi^2/2 + Euler^2)

# Derivation:
# psi2(1)*gamma(1) + 2*psi1(1)*d(gamma)(1) + psi(1)*d2(gamma)(1)
- 2* pracma::zeta(3) - 2*pi^2/6 * Euler - Euler*(pi^2/6 + Euler^2)


### I( x^n * log(x) * exp( - x^n) )
n = sqrt(5)
integrate(\(x) x^n * log(x) * exp( - x^n), 0, Inf)
gamma(1/n) * digamma(1/n + 1) / n^3


### Gen: I( x^p * log(x) * exp( - x^n) )
p = sqrt(5); n = sqrt(3);
integrate(\(x) x^p * log(x) * exp( - x^n), 0, Inf)
gamma((p+1)/n) * digamma((p+1)/n) / n^2;


#####################

### I( x^p * log(x) * (exp(-a1*x) - exp(-a2*x)) )
p = -1/3; # p != -1 && p > - 2 !!!
a = sqrt(c(5, 2))
integrate(\(x) x^p * log(x) * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(p + 1) * digamma(p + 1) * (a[1]^(-p-1) - a[2]^(-p-1)) +
	- gamma(p + 1) * (a[1]^(-p-1)*log(a[1]) - a[2]^(-p-1)*log(a[2]))

###
p = - sqrt(3)
a = sqrt(c(5, 2))
integrate(\(x) x^p * log(x) * (exp(-a[1]*x) - exp(-a[2]*x)), 0, Inf)
gamma(p + 1) * digamma(p + 1) * (a[1]^(-p-1) - a[2]^(-p-1)) +
	- gamma(p + 1) * (a[1]^(-p-1)*log(a[1]) - a[2]^(-p-1)*log(a[2]))


### I( log(x) * (exp(-a[1]*x) - exp(-a[2]*x)) / x )
# Lim: p -> -1
a = sqrt(c(5, 2))
integrate(\(x) log(x) * (exp(-a[1]*x) - exp(-a[2]*x)) / x, 0, Inf)
log(a[1]/a[2]) * (log(a[1]) + log(a[2]) + 2*Euler)/2


#####################
#####################

### I( log(x) / (exp(x) + 1) )
# Maths 505: ONE OF THE COOLEST INTEGRALS EVER!!! int ln(x)/(1+e^x) from 0 to infty
# https://www.youtube.com/watch?v=qY_sLn8yYLM
# 2.) Michael Penn: a stylized integral [alternative method]
# https://www.youtube.com/watch?v=_orPdt5r1Yg


integrate(\(x) log(x) / (exp(x) + 1), 0, Inf)
- log(2)^2/2


### Gen 1:
integrate(\(x) log(x)^2 / (exp(x) + 1), 0, Inf, rel.tol=1E-8)
log(2)*(pi^2/6 - Euler^2) - 2*gStjelt1*log(2) + log(2)^3 / 3


id = seq(120000)
log(2)*(pi^2/6 - Euler^2) + Euler*log(2)^2 +
	- sum((-1)^id * log(id)^2 / id)
# TODO: d2 eta(s) (1);
# after 8:00 in the presentation;
# d eta(s) = (1 - 2^(1-s)) * d zeta(s) + log(2)*2^(1-s) * zeta(s);
# d2 eta(s) = (1 - 2^(1-s)) * d2 zeta(s) +
#   + 2*log(2)*2^(1-s) * d zeta(s) - log(2)^2 * zeta(s);


### Gen 2: I( log(x) / (exp(k*x) + 1) )
integrate(\(x) log(x)/(exp(2*x) + 1), 0, Inf)
- log(2)^2 * 3/4

###
k = 3
integrate(\(x) log(x)/(exp(k*x) + 1), 0, Inf)
- log(2)^2/(2*k) - log(2)*log(k)/k

###
# up = Inf; numerical issue!
integrate(\(x) log(x) * (exp(x) - 2) / (exp(2*x) - exp(x) + 1), 0, 100)
log(2)*log(3)


### I( log(x) * exp(x) / (exp(2*x) + 1) )
# up = Inf; numerical issue!
integrate(\(x) log(x) / (exp(x) + exp(-x)), 0, 20)
integrate(\(x) log(x) / cosh(x) / 2, 0, 100)
pracma::integral(\(x) log(x) / cosh(x) / 2, 0, Inf)
id = seq(160000)
- pi*Euler / 4 - sum((-1)^id * log(2*id + 1) / (2*id + 1))
# TODO: how?


### I( x * log(x) * exp(k*x) / (exp(k*x) + 1)^2 )
k = sqrt(5)
# up = Inf; numerical issue!
integrate(\(x) x * log(x) * exp(k*x) / (exp(k*x) + 1)^2, 0, 30)
- log(2)^2/(2*k^2) + log(2)*(1 - log(k))/k^2

### I( x * log(x) / (exp(k*x) + 1)^2 )
k = 1/sqrt(5)
# up = Inf; numerical issue!
integrate(\(x) x * log(x) / (exp(k*x) + 1)^2, 0, 30)
(1 - Euler - log(k/2)) * pi^2 / (12*k^2) +
	+ (log(2)^2 - 2*log(2) + 2*log(2)*log(k) + dzeta2) / (2*k^2);


### I( x^2 * log(x) * exp(k*x) / (exp(k*x) + 1)^3 )
k = 1/sqrt(5)
# up = Inf; numerical issue!
integrate(\(x) x^2 * log(x) * exp(k*x) / (exp(k*x) + 1)^3, 0, 30)
(log(2)^2 - 3*log(2) + dzeta2 + pi^2*(log(2) + 3/2 - Euler)/6) / (2*k^3) +
	+ (12*log(2) - pi^2) * log(k) / (12*k^3);


### I( x * log(x) / (exp(k*x) - 1) )

###
integrate(\(x) x * log(x) / (exp(x) - 1), 0, Inf)
(1 - Euler)*pi^2 / 6 + dzeta2;
# TODO: dzeta2 ???

###
k = sqrt(5)
integrate(\(x) x * log(k*x) / (exp(k*x) - 1), 0, Inf)
((1 - Euler)*pi^2 / 6 + dzeta2) / k^2;

# =>
integrate(\(x) x * log(x) / (exp(k*x) - 1), 0, Inf)
((1 - Euler)*pi^2 / 6 + dzeta2) / k^2 - log(k) * gamma(2) * pracma::zeta(2)/k^2;
((1 - Euler)*pi^2 / 6 + dzeta2) / k^2 - pracma::zeta(2) * log(k) / k^2;


### Glaisher–Kinkelin Constant:
k = 2*pi;
integrate(\(x) x * log(x) / (exp(k*x) - 1), 0, Inf)
((1 - Euler)*pi^2 / 6 + dzeta2) / k^2 - log(k) * gamma(2) * pracma::zeta(2)/k^2;
#
A = exp(1/12 - ((1 - Euler)*pi^2 / 6 + dzeta2 - log(2*pi) * pracma::zeta(2))/(2*pi^2));
A = exp((log(2*pi) + Euler - 6*dzeta2/pi^2)/12);


### I( x * log(x) / (exp(k*x) + 1) )
# =>
k = 1/sqrt(3)
integrate(\(x) x * log(x) / (exp(k*x) + 1), 0, Inf)
(1 - Euler + dzeta2 * 6/pi^2 - log(k/2)) * pi^2 / (12*k^2);

### Special Case:
integrate(\(x) x * log(x) / (exp(2*x) + 1), 0, Inf)
(1 - Euler + dzeta2 * 6/pi^2) * pi^2 / 48;


### Gen: I( x^p * log(x) / (exp(k*x) + 1) )
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p * log(x) / (exp(k*x) + 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * digamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) / k^(p + 1) +
	+ gamma(p + 1) * dzeta(p + 1) * (1 - 1/2^p) / k^(p + 1) +
	+ gamma(p + 1) * pracma::zeta(p + 1) * log(2) / (2^p * k^(p + 1)) +
	- gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) * log(k) / k^(p + 1);


### Special Case:
# I( x^2 * log(x) / (exp(k*x) + 1) )
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^2 * log(x) / (exp(k*x) + 1), 0, Inf, rel.tol=1E-8)
(3*digamma(3) + log(2) - 3*log(k)) * pracma::zeta(3) / (2*k^3) +
	+ 3/2 * dzeta(3) / k^3;


### I( x^p * log(x) / (exp(k*x) + 1)^2 )
k = sqrt(3) - sqrt(2); p = sqrt(5);
integrate(\(x) x^p * log(x) / (exp(k*x) + 1)^2, 0, Inf)
gamma(p + 1) * digamma(p + 1) *
	( pracma::zeta(p + 1) * (2^p - 1) +
	- pracma::zeta(p) * (2^p - 2) ) / (2^p * k^(p+1)) +
	+ gamma(p + 1) *
	( dzeta(p + 1) * (2^p - 1) - dzeta(p) * (2^p - 2) +
	+ (pracma::zeta(p + 1) - pracma::zeta(p)) * 2^p * log(2) ) / (2^p * k^(p+1)) +
	- gamma(p + 1) *
	( pracma::zeta(p + 1) * (2^p - 1) +
	- pracma::zeta(p) * (2^p - 2) ) * log(2*k) / (2^p * k^(p+1));


### Pow: x^2

### I( x^2 * log(x) * exp(k*x) / (exp(k*x) - 1)^2 )
k = sqrt(5)
# Up = Inf; numerical issue!
pracma::integral(\(x) x^2 * log(x) * exp(k*x) / (exp(k*x) - 1)^2, 0, 100)
((3/2 - Euler)*pi^2 / 3 + 2*dzeta2) / k^3 +
	- 2 * pracma::zeta(2) * log(k) / k^3;


### I( x^2 * log(x) * exp(k*x) / (exp(k*x) + 1)^2 )
k = sqrt(5)
# Up = Inf; numerical issue!
pracma::integral(\(x) x^2 * log(x) * exp(k*x) / (exp(k*x) + 1)^2, 0, 100)
((3/2 - Euler)*pi^2 + 6*dzeta2) / (6*k^3) +
	- pi^2 * log(k/2) / (6*k^3);


### I( x^2 * log(x) / (exp(k*x) + 1)^2 )
# =>
k = 1/sqrt(5)
# Up = Inf; numerical issue!
pracma::integral(\(x) x^2 * log(x) / (exp(k*x) + 1)^2, 0, 100)
(3*digamma(3) + log(2) - 3*log(k)) * pracma::zeta(3) / (2*k^3) +
	- ((3/2 - Euler)*pi^2/6 + dzeta2 - 3/2*dzeta(3)) / k^3 +
	+ pi^2 * log(k/2) / (6*k^3);
# Note: digamma(3) = 3/2 - Euler;

