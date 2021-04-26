########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Trigonometric Solutions
###
### draft v.0.1a


### Trigonometric Polynomials


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


##########################

###############
### Order 4 ###
###############

x^4 - 2*x^3 + 2*x^2 - x - R # = 0

### Solution

### sin(a)^10 + cos(a)^10 = 5*R + 1
(sin(a)^2 + cos(a)^2)^5 - 5*sin(a)^2 * cos(a)^2 + 5*sin(a)^4 * cos(a)^4 - 1 - 5*R # = 0
sin(a)^4 * cos(a)^4 - sin(a)^2 * cos(a)^2 - R # = 0
y = sin(2*a)^2 # =>
y^2 - 4*y - 16*R # = 0
xsqrt = sin(asin(sqrt(y)) / 2)^2;
x = xsqrt^2;

### Solver:
solve.Trig.P4 = function(R) {
	y = roots(c(1, -4, -16*R));
	y = sqrt(y + 0i); y = c(y, -y);
	print(y)
	a = asin(y) / 2;
	x = sin(a)^2;
	# TODO: all roots; HOW ???
	return(x);
}

### Examples:
R = 1
x = solve.Trig.P4(R);

### Test
err = x^4 - 2*x^3 + 2*x^2 - x - R
round0(err)


