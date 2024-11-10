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
# Note: Series expansion + Contour on half of Unit Circle;
# Pole of order k: Res = 1/k!;

integrate(\(x) exp(2*cos(x)), 0, pi)
integrate(\(x) exp(2*cos(x)) / 2, -pi, pi)
pi * sum(1 / factorial(seq(0, 15))^2)

# TODO: closed formula?


###
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

###
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

