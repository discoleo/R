

### Integrals: Log & Exp


### Examples:

### on [0, Inf]
# I( x^p * log(x) / (exp(k*x) - 1) )
# I( x^p * log(x) / (exp(k*x) + 1) )
# I( x^p * log(x) * (exp(-a1*x) - exp(-a2*x)) )
# I( log(x^2 + b^2) / cosh(x) ) & Variants:
#   I( log(x^4 + b^2*x^2 + b^4) / cosh(x) )
#   I( log(x^4 - b^2*x^2 + b^4) / cosh(x) )
# I( log(x^2 + b^2) / (cosh(k*x) + cos(phi)) )
# I( log(x^2 + 1) / (sinh(3*k*x) * sinh(k*x)) )


####################

### Helper Functions

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
# Note: generalizations are a few sections below;


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
integrate(\(x) log(x) / (exp(2*x) + 1), 0, Inf)
- log(2)^2 * 3/4

###
k = 3
integrate(\(x) log(x) / (exp(k*x) + 1), 0, Inf)
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


#######################
### Generalizations ###

### I( x^p * log(x) / (exp(k*x) - 1) )
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p * log(x) / (exp(k*x) - 1), 0, Inf, rel.tol=1E-8)
((digamma(p+1) - log(k)) * pracma::zeta(p + 1) + dzeta(p+1)) *
	gamma(p + 1) / k^(p+1)


### I( x^p * log(x) / (exp(k*x) + 1) )
p = sqrt(5); k = sqrt(3);
integrate(\(x) x^p * log(x) / (exp(k*x) + 1), 0, Inf, rel.tol=1E-8)
gamma(p + 1) * digamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) / k^(p + 1) +
	+ gamma(p + 1) * dzeta(p + 1) * (1 - 1/2^p) / k^(p + 1) +
	+ gamma(p + 1) * pracma::zeta(p + 1) * log(2) / (2^p * k^(p + 1)) +
	- gamma(p + 1) * pracma::zeta(p + 1) * (1 - 1/2^p) * log(k) / k^(p + 1);


### Special Case:
# I( x^2 * log(x) / (exp(k*x) + 1) )
k = sqrt(3);
integrate(\(x) x^2 * log(x) / (exp(k*x) + 1), 0, Inf, rel.tol=1E-8)
(3*digamma(3) + log(2) - 3*log(k)) * pracma::zeta(3) / (2*k^3) +
	+ 3/2 * dzeta(3) / k^3;


### Pow: ()^2

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


#######################
#######################

### I( log(x^2 + b^2) / HYP )
# Blagouchine IV. Rediscovery of Malmsten’s integrals, their evaluation
# by contour integration methods and some related results.
# Ramanujan J (2014) 35:21–110; DOI: I 10.1007/s11139-013-9528-5

### I( log(x^2 + b^2) / cosh(x) ) on [0, Inf]
# see page 12 of article;
b = sqrt(5);
integrate(\(x) log(x^2 + b^2) / cosh(x), 0, Inf)
2*pi * log(gamma(b/(2*pi) + 3/4) / gamma(b/(2*pi) + 1/4) * sqrt(2*pi))

###
b = sqrt(5); k = sqrt(5);
integrate(\(x) log(x^2 + b^2) / cosh(k*x), 0, Inf)
2*pi/k * log(gamma(k*b/(2*pi) + 3/4) / gamma(k*b/(2*pi) + 1/4) * sqrt(2*pi/k));


### Variants:

### I( log(x^4 + b^4) / cosh(x) )
b = sqrt(5)
integrate(\(x) log(x^4 + b^4) / cosh(x), 0, Inf)
4*pi * Re(log(pracma::gammaz(b * exp(pi/4*1i)/(2*pi) + 3/4) /
	pracma::gammaz(b * exp(pi/4*1i)/(2*pi) + 1/4) * sqrt(2*pi)) );

### I( log(x^4 + b^2*x^2 + b^4) / cosh(x) )
b = sqrt(5)
integrate(\(x) log(x^4 + b^2*x^2 + b^4) / cosh(x), 0, Inf)
4*pi * Re(log(pracma::gammaz(b * exp(pi/6*1i)/(2*pi) + 3/4) /
	pracma::gammaz(b * exp(pi/6*1i)/(2*pi) + 1/4) * sqrt(2*pi)) );

### I( log(x^4 - b^2*x^2 + b^4) / cosh(x) )
b = sqrt(5)
integrate(\(x) log(x^4 - b^2*x^2 + b^4) / cosh(x), 0, Inf)
4*pi * Re(log(pracma::gammaz(b * exp(pi/3*1i)/(2*pi) + 3/4) /
	pracma::gammaz(b * exp(pi/3*1i)/(2*pi) + 1/4) * sqrt(2*pi)) );


### I( 1 / ((x^2 + b^2) * cosh(x)) ) on [0, Inf]
# - Derived I();
b = sqrt(5);
integrate(\(x) 1 / ((x^2 + b^2) * cosh(x)), 0, Inf)
(digamma(b/(2*pi) + 3/4) - digamma(b/(2*pi) + 1/4)) / (2*b)

#
integrate(\(x) 1 / ((x^2 + 1) * (log(x)^2 + b^2)), 1, Inf)
(digamma(b/(2*pi) + 3/4) - digamma(b/(2*pi) + 1/4)) / (4*b)


### I( log(x^2 + b^2) / (exp(x) + exp(-x) - 1) ) on [0, Inf]
# - see on page 48 of article (28 in pdf);
b = sqrt(5);
integrate(\(x) log(x^2 + b^2) / (exp(x) + exp(-x) - 1), 0, Inf)
4/3 * pi * sin(pi/3) * log(gamma(b/(2*pi) + 5/6) / gamma(b/(2*pi) + 1/6)) +
	+ 4*pi/3 / tan(pi/3) * log(2*pi);


### I( log(x^2 + b^2) / (exp(x) + exp(-x) + 1) ) on [0, Inf]
b = sqrt(5);
integrate(\(x) log(x^2 + b^2) / (exp(x) + exp(-x) + 1), 0, Inf)
4/3 * pi * sin(pi/3) * log(gamma(b/(2*pi) + 1 - 1/3) / gamma(b/(2*pi) + 1/3)) +
	+ 2*pi/3 / tan(pi/3) * log(2*pi);


# Derivation:
integrate(\(x) log(x^2 + b^2) / (exp(x*2*pi/3) + exp(-x*2*pi/3) + 1), 0, Inf)
integrate(\(x) 3/pi * (log(x^2 + (b*pi/3)^2) - 2*log(pi/3)) / (exp(2*x) + exp(-2*x) + 1), 0, Inf)
integrate(\(x) 3/pi * log(x^2 + (b*pi/3)^2) / (exp(2*x) + exp(-2*x) + 1), 0, Inf)$value +
	- log(pi/3) / tan(pi/3);
tan(pi/6)*log(3) + 2*sin(pi/3) * log(gamma(1 - 1/3 + b/3) / gamma((b+1)/3))


###
integrate(\(x) log(x^2 + b^2) / (exp(2*x) + exp(x) + exp(-x) + exp(-2*x) + 1), 0, Inf)
4/5 * pi * sin(pi/5) * log(gamma(1 - 1/5 + b/(2*pi)) / gamma(b/(2*pi) + 1/5)) +
	- 4/5 * pi * sin(2*pi/5) * log(gamma(1 - 2/5 + b/(2*pi)) / gamma(b/(2*pi) + 2/5)) +
	+ 2*pi/5 / tan(2*pi/5) * log(2*pi);

#
integrate(\(x) 1 / (exp(2*x)+exp(x) + exp(-x) + exp(-2*x) + 1), 0, Inf)
pi/5 / tan(2*pi/5)


### I( log(x^2 + b^2) / (exp(x) + exp(-x) + 2) ) on [0, Inf]
# - see on page 49 of article (page 29 in pdf);
b = sqrt(5);
integrate(\(x) log(x^2 + b^2) / (exp(x) + exp(-x) + 2), 0, Inf)
integrate(\(x) 1/2 * log(x^2 + b^2) / (cosh(x) + 1), 0, Inf)
digamma(b/(2*pi) + 1/2) + log(2*pi);


### Gen: I( log(x^2 + b^2) / (cosh(k*x) + cos(phi)) ) on [0, Inf]
# - see Exercise 2 on page 50 of article (page 30 in pdf);
b = sqrt(5); k = sqrt(3); phi = 1/3;
integrate(\(x) log(x^2 + b^2) / (cosh(k*x) + cos(phi)), 0, Inf)
2*pi/(k*sin(phi)) * log(gamma((k*b + phi)/(2*pi) + 1/2) /
		gamma((k*b - phi)/(2*pi) + 1/2)) +
	+ 2*phi/(k*sin(phi)) * log(2*pi/k);

# [Practice] Variants:
b = sqrt(5); k = sqrt(3); phi = 1/3;
#
integrate(\(x) log(x^2 + b^2) / (cosh(2*k*x) + cosh(k*x) + cos(phi)), 0, Inf)
integrate(\(x) log(x^2 + b^2) / (cosh(3/2*k*x) * cosh(k/2*x) * 2 + cos(phi)), 0, Inf)
dd = sqrt(1/4 + 4/2 - 4/2*cos(phi) + 0i); phi0 = acos(1/4 + c(1,-1)*dd / 2);
diff( 2*pi/(k*sin(phi0)) * log(pracma::gammaz((k*b + phi0)/(2*pi) + 1/2) /
		pracma::gammaz((k*b - phi0)/(2*pi) + 1/2)) +
	+ 2*phi0/(k*sin(phi0)) * log(2*pi/k) ) / (2*dd);
# Upper = Inf; # Numeric instability!
integrate(\(x) log(x^2 + b^2) / (cosh(2*k*x) - cosh(k*x) + cos(phi)), 0, 120)
integrate(\(x) log(x^2 + b^2) / (sinh(3/2*k*x) * sinh(k/2*x) * 2 + cos(phi)), 0, 120)
dd = sqrt(1/4 + 4/2 - 4/2*cos(phi) + 0i); phi0 = acos(-1/4 + c(1,-1)*dd / 2);
diff( 2*pi/(k*sin(phi0)) * log(pracma::gammaz((k*b + phi0)/(2*pi) + 1/2) /
		pracma::gammaz((k*b - phi0)/(2*pi) + 1/2)) +
	+ 2*phi0/(k*sin(phi0)) * log(2*pi/k) ) / (2*dd);


### I( log(x^2 + b^2) / (cosh(3*k*x) * cosh(k*x)) )
# Special Case: Limit phi = pi/2
integrate(\(x) log(x^2 + b^2) / (cosh(2*k*x) + cosh(k*x)), 0, Inf)
integrate(\(x) log(x^2 + b^2) / (cosh(3/2*k*x) * cosh(k/2*x)) / 2, 0, Inf)
phi0 = 2/3 * pi; # acos(-1/2); dd = 3/2;
2/3 * pi/(k*sin(phi0)) * log(pracma::gammaz((k*b + phi0)/(2*pi) + 1/2) /
		pracma::gammaz((k*b - phi0)/(2*pi) + 1/2)) +
	+ 2/3 * phi0/(k*sin(phi0)) * log(2*pi/k) +
	- 1/3 * (2*digamma(k*b/(2*pi) + 1/2) / k + 2/k * log(2*pi/k));


### I( log(x^2 + 1) / (sinh(3*k*x) * sinh(k*x)) )
# Special Case: Limit phi = pi/2 & b = 1;
k = sqrt(3) - sqrt(2);
integrate(\(x) log(x^2 + 1) / (sinh(3*k*x) * sinh(k*x)), 0, Inf)
integrate(\(x) log(x^2 + 1) / (cosh(4*k*x) - cosh(2*k*x)) * 2, 0, Inf)
# dd = 3/2; phi0 = c(pi, pi/3);
(2*pi * log(pracma::gammaz(k/pi + 1/3) / pracma::gammaz(k/pi + 2/3)) +
	+ 2*pi/3 * log(k/pi)) / (3*k*sin(pi/3)) +
	- 2*(digamma(k/pi) + pi/(2*k) + log(pi/k)) / (3*k);


#########
### SINH:

### Gen: I( log(x^2 + b^2) / (sinh(k*x) + sinh(phi)) ) on [0, Inf]
# - see Exercise 4a on page 50 of article (page 30 in pdf);
# - formula for principal value;

# TODO:

# library(Rmpfr)
b = sqrt(mpfr(3, 240)); phi = 4/mpfr(7, 240);
# Upper = Inf; # Numerical issues!
integrate(\(x) as.numeric(log(x^2 + b^2) / (sinh(x) + sinh(phi))),
	0, 150, rel.tol=1E-8)$value +
integrate(\(x) {
	x = mpfr(x, 240);
	y = log(x^2 + b^2) * 2*exp(x) / (1 - exp(x-phi)) / (1 + exp(x+phi)) +
		+ log(b^2 + phi^2) / cosh(phi) / (exp(abs(x - phi)) - 1);
	as.numeric(y); },
	phi, 150, rel.tol=1E-8)$value +
integrate(\(x) {
	x = mpfr(x, 240);
	y = log(x^2 + b^2) * 2*exp(x)/(1 - exp(x-phi)) / (1 + exp(x+phi)) +
		- log(b^2 + phi^2) / cosh(phi) / (exp(abs(x - phi)) - 1);
	as.numeric(y); }, 0, phi)$value

# OK
b = as.numeric(b); phi = as.numeric(phi);
4*pi/cosh(phi) * Im( log(pracma::gammaz((b + phi*1i)/(2*pi)) /
		pracma::gammaz((b - phi*1i)/(2*pi) + 1/2)) ) +
	+ 2*pi/cosh(phi) * (pi/2 - atan(b/phi)) + 4*phi*log(2*pi)/cosh(phi) +
	# Extra Limit-Terms:
	- log(b^2 + phi^2) / cosh(phi) * (log(exp(phi) - 1) - phi);


####################
####################

### I( log(x/b) / (sinh(x)^2 - sinh(b)^2) )
# Exercise 5 in Blagouchine

b = sqrt(5)
integrate(\(x) log(x/b) / (sinh(x)^2 - sinh(b)^2), 0, Inf)
- 2*pi / sinh(2*b) * Im(log(pracma::gammaz(1i*b/(2*pi)) /
		pracma::gammaz(1/2 - 1i*b/(2*pi)))) +
	+ (2*b*log(b / (2*pi)) - pi^2/2) / sinh(2*b);


####################

### Higher Powers:

### I( log(x^2 + b^2) / cosh(k*x)^2 )
# Exercise 6 in Blagouchine

b = sqrt(5); k = sqrt(3);
integrate(\(x) log(x^2 + b^2) / cosh(k*x)^2, 0, Inf)
2/k * (log(pi/k) + digamma(b*k/pi + 1/2))

