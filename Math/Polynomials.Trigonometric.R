########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Trigonometric Solutions
###
### draft v.0.1b


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
	y = sqrt(y + 0i); # y = c(y, -y);
	print(y)
	a = asin(y);
	a = c(a, pi-a) / 2;
	x = sin(a)^2;
	# alternative: x^2 - x + y^2/4 = 0;
	# duplicate roots: filtered!
	return(x);
}

### Examples:
R = 1
x = solve.Trig.P4(R);

### Test
err = x^4 - 2*x^3 + 2*x^2 - x - R
round0(err)


##########################
##########################

###############
### Order 6 ###
###############

2*x^6 - 6*x^5 + 15*x^4 - 20*x^3 + 15*x^2 - 6*x - R # = 0

### Solution

### sin(a)^12 + cos(a)^12 = R + 1
(sin(a)^2 + cos(a)^2)^6 - 6*sin(a)^2 * cos(a)^2 + 9*sin(a)^4 * cos(a)^4 - 2*sin(a)^6 * cos(a)^6 - R - 1 # = 0
2*sin(a)^6 * cos(a)^6 - 9*sin(a)^4 * cos(a)^4 + 6*sin(a)^2 * cos(a)^2 + R # = 0
y = sin(2*a)^2 # =>
y^3 - 18*y^2 + 48*y + 32*R # = 0
xsqrt = sin(asin(sqrt(y)) / 2)^2;
x = xsqrt^2;

### Solver:
solve.Trig.P6 = function(R) {
	y = roots(c(1, -18, 48, 32*R));
	y = sqrt(y + 0i); # y = c(y, -y);
	### asin()-variant
	# a = asin(y);
	# a = c(a, pi-a) / 2;
	# x = sin(a)^2;
	# alternative: x^2 - x + y^2/4 = 0;
	x = sapply(y, function(y) roots(c(1, -1, y^2/4)));
	# duplicate roots: filtered!
	return(as.vector(x));
}

### Examples:
R = 2
x = solve.Trig.P6(R);

### Test
err = 2*x^6 - 6*x^5 + 15*x^4 - 20*x^3 + 15*x^2 - 6*x - R
round0(err)

