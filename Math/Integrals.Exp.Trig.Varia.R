#########################
##
## Integrals: Trig( EXP )
##
## Leonard Mada


### Unusual Types


#####################

### I( exp(k * cos(x)) ) on [0, pi]
# Maths 505: A surprisingly interesting integral
# https://www.youtube.com/watch?v=81qExKEYzo0
# Note: Series expansion of exp(exp(-1i*x));
# + Contour on half of Unit Circle;
# Pole of order k: Res = 1/k!;

### I( exp(2*cos(x)) )
integrate(\(x) exp(2*cos(x)), 0, pi)
integrate(\(x) exp(2*cos(x)) / 2, -pi, pi)
pi * sum(1 / factorial(seq(0, 15))^2)
pi * besselI(2,0);

### Gen: I( exp(k*cos(x)) )
k = 1/pi;
integrate(\(x) exp(k*cos(x)), 0, pi)
pi * besselI(k,0);

# TODO: closed formula for Bessel?


### Complex:

### I( exp(2i*sin(x)) )
integrate(\(x) Re(exp(2i*sin(x))), 0, pi)
pi * besselJ(2,0);

### Im(...)
integrate(\(x) Im(exp(2i*sin(x))), 0, pi, rel.tol=1E-13)
# TODO


### Gen: I( exp(k*1i * sin(x)) )
k = sqrt(pi);
integrate(\(x) Re(exp(k*1i*sin(x))), 0, pi)
pi * besselJ(k,0)


# Note:
sum(1 / factorial(seq(0, 15))^2);
besselI(2,0)
#
sum((-1)^seq(0, 15) / factorial(seq(0, 15))^2);
besselJ(2,0)


### I( exp(-k * tan(x)) )
# Maths 505: This is WILD!
# https://www.youtube.com/watch?v=LeR2ELsT_bI
# Note: Feynman trick (twice)
# => d2y + dy = I( exp(-k*tan(x)) * (1+tan(x)^2) ) = 1/k;

integrate(\(x) exp(-tan(x)), 0, pi/2, rel.tol=1E-13)
pracma::Ci(1) * sin(1) - pracma::Si(1) * cos(1) + pi/2 * cos(1);

### Gen: I( exp(-k * tan(x)) )
k = 1/3;
integrate(\(x) exp(-k * tan(x)), 0, pi/2, rel.tol=1E-13)
pracma::Ci(k) * sin(k) - pracma::Si(k) * cos(k) + pi/2 * cos(k);


###################

### Bessel Function

### Properties:
besselY(1, 2) - 2*besselY(1, 1) + besselY(1, 0) # == 0
besselY(1, 3) - 4*besselY(1, 2) + besselY(1, 1) # == 0
besselY(1, 4) - 6*besselY(1, 3) + besselY(1, 2) # == 0
#
n = 3/7;
besselY(1, n+1) - 2*n*besselY(1, n) + besselY(1, n-1) # == 0

###
besselY(2, 2) -   besselY(2, 1) + besselY(2, 0) # == 0;
besselY(2, 3) - 2*besselY(2, 2) + besselY(2, 1) # == 0;
besselY(2, 4) - 3*besselY(2, 3) + besselY(2, 2) # == 0;
#
n = 3/7
besselY(2, n+1) - n*besselY(2, n) + besselY(2, n-1) # == 0;

###
besselY(3, 2) - 2/3 * besselY(3, 1) + besselY(3, 0) # == 0
besselY(3, 3) - 4/3 * besselY(3, 2) + besselY(3, 1) # == 0
n = 5/7;
besselY(3, n+1) - 2*n/3 * besselY(3, n) + besselY(3, n-1) # == 0

### Gen:
n = 5/7; x = 4/3;
besselY(x, n+1) - 2*n/x * besselY(x, n) + besselY(x, n-1) # == 0
besselJ(x, n+1) - 2*n/x * besselJ(x, n) + besselJ(x, n-1) # == 0


# using Rmpfr:
# prec = 190; x2 = mpfr(2, prec);
# jn(2, 2) - jn(1, 2) + jn(0, 2)
# jn(3, x2) - jn(2, x2) + jn(0, x2)


### Values

### BesselJ:
besselJ(1, -1/2)
cos(1) / sqrt(pi/2);
#
besselJ(1, 1/2)
sin(1) / sqrt(pi/2);

###
besselJ(2, -1/2)
cos(2) / sqrt(pi);
#
besselJ(2, 1/2)
sin(2) / sqrt(pi);

### Gen: BesselJ(x, +/- 1/2)
x = sqrt(pi);
besselJ(x, -1/2)
cos(x) / sqrt(x*pi/2);
#
besselJ(x, 1/2)
sin(x) / sqrt(x*pi/2);


### BesselI:
x = sqrt(pi-2);
besselI(x, 1/2)
Im(sin(x*1i)) / sqrt(x*pi/2);
sinh(x) / sqrt(x*pi/2);
#
besselI(x, -1/2)
cosh(x) / sqrt(x*pi/2);

# Recurrence:
x = pi; n = 4/3;
besselI(x, n+1) # ==
-2*n/x * besselI(x, n) + besselI(x, n-1);


### Extension:
x = sqrt(5);
besselJ(x, 3/2)
(sin(x) - x*cos(x)) / sqrt(x^3*pi/2);
#
besselJ(x, -3/2)
- (x*sin(x) + cos(x)) / sqrt(x^3*pi/2);


### BesselY:
x = sqrt(5)
besselY(x, 1/2)
- besselJ(x, -1/2)
#
besselY(x, -1/2)
besselJ(x, 1/2)


### Extension: BesselY
x = sqrt(5);
besselY(x, 3/2)
- (x*sin(x) + cos(x)) / sqrt(x^3*pi/2);
#
besselY(x, -3/2)
- (sin(x) - x*cos(x)) / sqrt(x^3*pi/2);


### Bessel(k * 1i)
k = sqrt(3)
integrate(\(x) Re(cos(k * 1i * sin(x))), 0, pi/2)
pi/2 * besselI(k, 0);


### Other: Bessel Y vs Bessel J
k = sqrt(3);
n = sqrt(5) - sqrt(3);
besselJ(k, n+1) * besselY(k, n) - besselJ(k, n) * besselY(k, n+1) - 2/(k*pi) # == 0


### Errors

###
curve(besselJ(x, pi/2-x), 0, 2*pi, ylim = c(-2,2), n = 101)
abline(v = pi/2, col = "red", lty = 2)

besselJ(pi/2, 0)
# 0.4720012
besselJ(pi/2, 1E-14)
# 0.4720012
besselJ(pi/2, 1E-15)
# 4.720012e+14 # ERROR!

curve(besselJ(x, cos(x)), 0, 2*pi, ylim = c(-2,2), n=201)
abline(v = c(1,3) * pi/2, col = "red", lty = 2)

# Almost all values are affected:
median(besselJ(seq(0,1, by = 1/1024), 1E-15))


##################

### Bessel-Type

### I( cos(sin(x)) )
# Hmath: Интеграл с функцией Бесселя и трюками с рядами
# https://www.youtube.com/watch?v=JF1ikXax5rk
# [in Russian]

###
integrate(\(x) cos(sin(x)), 0, pi/2)
pi * besselJ(1, 0) / 2

### I( cos(k * sin(x)) )
k = sqrt(3)
integrate(\(x) cos(k * sin(x)), 0, pi/2)
pi * besselJ(k, 0) / 2


### I( Bessel( TRIG ) )

### Gen: I( BesselJn(k * sin(x)) )
k = exp(2/3); n = 3^(1/3);
integrate(\(x) besselJ(k*sin(x), n), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(k/2, n/2)^2;

### Gen: I( BesselIn(k * sin(x)) )
k = exp(2/3); n = 3^(1/3); # n != negative odd Integer;
integrate(\(x) besselI(k*sin(x), n), 0, pi/2, rel.tol=1E-13)
pi/2 * besselI(k/2, n/2)^2;

# Special Cases:
k = exp(-1/3);
n = -1; # n = -3;
# n = negative Integers; (can be even & odd)
integrate(\(x) besselI(k*sin(x), n), 0, pi/2, rel.tol=1E-13)
pi/2 * besselI(k/2, -n/2)^2;


### Explicit Cases:

### Gen: I( BesselJ0(k * sin(x)) )
k = exp(1/sqrt(2));
integrate(\(x) besselJ(k*sin(x), 0), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(k/2, 0)^2;


### Gen: I( BesselJ1(k * sin(x)) )
k = 3^(1/3);
integrate(\(x) besselJ(k*sin(x), 1), 0, pi/2, rel.tol=1E-13)
(1 - cos(k)) / k; # 2/k * sin(k/2)^2;
pi/2 * besselJ(k/2, 1/2)^2;


### I( BesselJ2(sin(x)) )
integrate(\(x) besselJ(sin(x), 2), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(1/2, 1)^2;

### Gen: I( BesselJ2(k * sin(x)) )
k = 3^(1/3);
integrate(\(x) besselJ(k*sin(x), 2), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(k/2, 1)^2;


### I( Bessel / Trig )

### Gen: I( BesselJn(k * sin(x)) / sin(x) )
k = exp(2/3); n = 3^(1/pi);
integrate(\(x) besselJ(k * sin(x), n) / sin(x), 0, pi/2, rel.tol=1E-13)
k/(4*n) * pi * (besselJ(k/2, (n+1)/2)^2 + besselJ(k/2, (n-1)/2)^2);

### Gen: I( BesselIn(k * sin(x)) / sin(x) )
k = exp(2/3); n = 3^(1/pi); # n > 0;
integrate(\(x) besselI(k * sin(x), n) / sin(x), 0, pi/2, rel.tol=1E-13)
k/(4*n) * pi * (besselI(k/2, (n-1)/2)^2 - besselI(k/2, (n+1)/2)^2);


### I( Trig * Bessel )

### Gen: I( sin(x) * BesselJn(k * sin(x)) )
k = sqrt(pi); n = exp(-1/3);
integrate(\(x) sin(x) * besselJ(k*sin(x), n), 0, pi/2)
pi/2 * besselJ(k/2, (n-1)/2) * besselJ(k/2, (n+1)/2);

### Gen: I( sin(x) * BesselIn(k * sin(x)) )
k = sqrt(pi); n = exp(-1/3); # n != negative even integer;
integrate(\(x) sin(x) * besselI(k*sin(x), n), 0, pi/2)
pi/2 * besselI(k/2, (n-1)/2) * besselI(k/2, (n+1)/2);

# Special Cases:
n = -2; # n = -4; # n = actually any negative integer;
k = exp(sqrt(2));
integrate(\(x) sin(x) * besselI(k*sin(x), n), 0, pi/2)
pi/2 * besselI(k/2, -(n-1)/2) * besselI(k/2, -(n+1)/2);


### Explicit Cases:

### I( sin(x) * BesselJ0(sin(x)) )
integrate(\(x) sin(x) * besselJ(sin(x), 0), 0, pi/2)
sin(1);

### I( sin(x) * BesselJ0(2*sin(x)) )
integrate(\(x) sin(x) * besselJ(2*sin(x), 0), 0, pi/2)
sin(2) / 2;

### Gen: I( sin(x) * BesselJ0(k*sin(x)) )
k = sqrt(pi);
integrate(\(x) sin(x) * besselJ(k*sin(x), 0), 0, pi/2)
sin(k) / k;
pi/2 * besselJ(k/2, -1/2) * besselJ(k/2, 1/2);


### I( sin(x) * BesselJ1(sin(x)) 
integrate(\(x) sin(x) * besselJ(sin(x), 1), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(1/2, 0) * besselJ(1/2, 1);

### Gen: I( sin(x) * BesselJ1(k*sin(x)) )
k = log(pi);
integrate(\(x) sin(x) * besselJ(k*sin(x), 1), 0, pi/2, rel.tol=1E-13)
pi/2 * besselJ(k/2, 0) * besselJ(k/2, 1);


### I( sin(x) * BesselJ2(sin(x)) )
integrate(\(x) sin(x) * besselJ(sin(x), 2), 0, pi/2, rel.tol=1E-13)
2*(1-cos(1)) - sin(1);

### Gen: I( sin(x) * BesselJ2(k * sin(x)) )
k = log(pi+1);
integrate(\(x) sin(x) * besselJ(k*sin(x), 2), 0, pi/2, rel.tol=1E-13)
(2 - 2*cos(k) - k*sin(k)) / k^2;


### I( sin(x) * BesselJ3(sin(x)) )
integrate(\(x) sin(x) * besselJ(sin(x), 3), 0, pi/2, rel.tol=1E-13)
4*pi/2 * besselJ(1/2, 1)^2 - pi/2 * besselJ(1/2, 0) * besselJ(1/2, 1);

### Gen: I( sin(x) * BesselJ3(k * sin(x)) )
k = exp(3/pi);
integrate(\(x) sin(x) * besselJ(k*sin(x), 3), 0, pi/2, rel.tol=1E-13)
2*pi/k * besselJ(k/2, 1)^2 - pi/2 * besselJ(k/2, 0) * besselJ(k/2, 1);
pi/2 * besselJ(k/2, 1) * (4/k * besselJ(k/2, 1) - besselJ(k/2, 0));


### I( sin(x)^2 * Bessel )

### Gen: I( sin(x)^2 * BesselJ0(k*sin(x)) )
k = exp(-2/3);
integrate(\(x) sin(x)^2 * besselJ(k*sin(x), 0), 0, pi/2, rel.tol=1E-13)
pi/4 * (besselJ(k/2, 0)^2 - besselJ(k/2, 1)^2);

### Gen: I( sin(x)^2 * BesselJ1(k*sin(x)) )
k = exp(-1);
integrate(\(x) sin(x)^2 * besselJ(k*sin(x), 1), 0, pi/2, rel.tol=1E-13)
besselJ(k, 3/2) * sqrt(pi/(2*k));


### BesselI

### Gen: I( sin(x) * BesselI0(k * sin(x)) )
k = sqrt(3)
integrate(\(x) sin(x) * besselI(k*sin(x), 0), 0, pi/2, rel.tol=1E-13)
sinh(k) / k;
besselI(k, 1/2) * sqrt(pi/(2*k));
pi/2 * besselI(k/2, -1/2) * besselI(k/2, 1/2);

# - see above for full generalisations;


### BesselY

### Gen: I( BesselY0(k * sin(x)) )
k = exp(1/2);
integrate(\(x) besselY(k*sin(x), 0), 0, pi/2, rel.tol=1E-13)
pi/2 * besselY(k/2, 0) * besselJ(k/2, 0);

### Gen: I( BesselY[1/2](k * sin(x)) )
k = exp(1/2);
integrate(\(x) besselY(k*sin(x), 1/2), 0, pi/2, rel.tol=1E-13)
- pi/2 * besselJ(k/2, -1/4)^2;

###
k = 1;
integrate(\(x) besselY(k*sin(x), 1) + 2/(k*pi) / x, 0, pi/2, rel.tol=1E-13)
# TODO


### I( sin(x) * BesselY0(sin(x)) )
integrate(\(x) sin(x) * besselY(sin(x), 0), 0, pi/2, rel.tol=1E-13)
(pracma::Ci(1) * sin(1) - pracma::Si(1) * cos(1)) * 2/pi;

### Gen: I( sin(x) * BesselY0(k*sin(x)) )
k = 3^(1/pi);
integrate(\(x) sin(x) * besselY(k*sin(x), 0), 0, pi/2, rel.tol=1E-13)
(pracma::Ci(k) * sin(k) - pracma::Si(k) * cos(k)) * 2 / (k*pi);

### Gen: I( sin(x) * BesselY[1/2](k*sin(x)) )
k = 3^(1/pi);
integrate(\(x) sin(x) * besselY(k*sin(x), 1/2), 0, pi/2, rel.tol=1E-13)
- pi/2 * besselJ(k/2, -3/4) * besselJ(k/2, 1/4);


### Bessel( cos(x) )

### I( sin(x) * besselJ0(cos(x)) )
integrate(\(x) sin(x) * besselJ(cos(x), 0), 0, pi/2)
integrate(\(x) besselJ(x, 0), 0, 1)
# TODO

### I( sin(x) * besselJ(cos(x), 1) )
integrate(\(x) sin(x) * besselJ(cos(x), 1), 0, pi/2)
1 - besselJ(1, 0);

### I( sin(x) * besselJ(cos(x), 2) )
integrate(\(x) sin(x) * besselJ(cos(x), 2), 0, pi/2)
integrate(\(x) sin(x) * besselJ(cos(x), 0), 0, pi/2)$value +
	- 2*besselJ(1, 1);


### I( x * Bessel )

### I( x * BesselJ0(x) )
integrate(\(x) x * besselJ(x, 0), 0, 1)
besselJ(1, 1);

### I( x * BesselJ1(x) )
integrate(\(x) x * besselJ(x, 1), 0, 1)
integrate(\(x) besselJ(x, 0), 0, 1)$value - besselJ(1, 0);


#########################

### I( exp(2*exp(x*1i)) )
integrate(\(x) Re(exp(2*exp(x*1i))), 0, pi)
pi

integrate(\(x) Im(exp(2*exp(x*1i))), 0, pi)
# TODO: ?

integrate(\(x) {
	x = mpfr(x, 240);
	y = exp(2*cos(x)) * sin(2*sin(x));
	as.numeric(y); }, 0, pi, rel.tol=1E-8)


######################

### I( exp(sin(2*x)/2) * cos(cos(x)^2) )
# Maths 505: A (literally) complex integral
# https://www.youtube.com/watch?v=4I6UgSwWkIA
# (intermediate transformation)

integrate(\(x) exp(sin(2*x)/2) * cos(cos(x)^2), 0, pi)
pi * cos(1/2)


############

### I( cos(x + sin(x)) * exp(cos(x)) / (x^2 + 1) )
# Hmath: Интегральный этюд с комплексным сюжетом
# [in Russian]
# => Re(exp(exp(1i*x)) * exp(1i*x))
# => sum( cos((n+1)*x) / (x^2+1) / n! );

integrate(\(x) cos(x + sin(x)) * exp(cos(x)) / (x^2 + 1), 0, Inf,
	subdivisions = 8000, rel.tol=1E-6)
exp(exp(-1) - 1) * pi/2;

# Derivation:
sapply(0:30, \(n) integrate(\(x, n) cos((n+1)*x) / (x^2 + 1) / gamma(n+1),
	0, Inf, n=n, subdivisions = 8000, rel.tol=1E-5)$value) |> sum();

# Var: sin
sapply(0:30, \(n) integrate(\(x, n) sin((n+1)*x) / (x^2 + 1) / gamma(n+1),
	0, Inf, n=n, subdivisions = 8000, rel.tol=1E-5)$value) |> sum();
# TODO


### I( Trig( x + sin(2*x) ) * exp(cos(2*x)) / Trig(x) )
# Maths 505: This wacky integral has a beautiful result
# https://www.youtube.com/watch?v=a2ZPqB2Syfo
# Note: series expansion + Dirichlet kernel;

### sin(...)
integrate(\(x) sin(x + sin(2*x)) * exp(cos(2*x)) / sin(x), 0, pi/2)
integrate(\(x) Im(exp(1i*x + exp(2i*x))) / sin(x), 0, pi/2)
pi/2 * exp(1)


### cos(...)
integrate(\(x) cos(x + sin(2*x)) * exp(cos(2*x)) / cos(x), 0, pi/2)
pi/2 * exp(-1)


### Dirichlet kernel:
k = 7; x = sqrt(5)
sin((2*k+1)*x) / sin(x) # ==
2*sum(cos(2*seq(k)*x)) + 1;

