#########################
##
## Integrals: Trig( EXP )
##
## Leonard Mada


### Unusual Types


#####################

### I( exp(2*cos(x)) ) on [0, pi]
# Maths 505: A surprisingly interesting integral
# https://www.youtube.com/watch?v=81qExKEYzo0
# Note: Series expansion of exp(exp(-1i*x));
# + Contour on half of Unit Circle;
# Pole of order k: Res = 1/k!;

### I( exp(2*cos(x)) )
integrate(\(x) exp(2*cos(x)), 0, pi)
integrate(\(x) exp(2*cos(x)) / 2, -pi, pi)
pi * sum(1 / factorial(seq(0, 15))^2)
pi * besselI(2,0)

# TODO: closed formula for Bessel?

### I( exp(2i*sin(x)) )
integrate(\(x) Re(exp(2i*sin(x))), 0, pi)
pi * besselJ(2,0)

###
integrate(\(x) Im(exp(2i*sin(x))), 0, pi)
# TODO

# Note:
sum(1 / factorial(seq(0, 15))^2);
besselI(2,0)
#
sum((-1)^seq(0, 15) / factorial(seq(0, 15))^2);
besselJ(2,0)


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


# using Rmpfr:
# prec = 190; x2 = mpfr(2, prec);
# jn(2, 2) - jn(1, 2) + jn(0, 2)
# jn(3, x2) - jn(2, x2) + jn(0, x2)


### Values
besselJ(1, -1/2)
cos(1) / sqrt(pi/2);


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

### I( Trig( x + sin(2*x) ) * exp(cos(2*x)) / Trig(x) )
# Maths 505: This wacky integral has a beautiful result
# https://www.youtube.com/watch?v=a2ZPqB2Syfo
# Note: series expansion + Dirichlet kernel;

###
integrate(\(x) sin(x + sin(2*x)) * exp(cos(2*x)) / sin(x), 0, pi/2)
integrate(\(x) Im(exp(1i*x + exp(2i*x))) / sin(x), 0, pi/2)
pi/2 * exp(1)


###
integrate(\(x) cos(x + sin(2*x)) * exp(cos(2*x)) / cos(x), 0, pi/2)
pi/2 * exp(-1)


### Dirichlet kernel:
k = 7; x = sqrt(5)
sin((2*k+1)*x) / sin(x) # ==
2*sum(cos(2*seq(k)*x)) + 1;

