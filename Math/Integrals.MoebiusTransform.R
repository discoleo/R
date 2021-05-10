########################
###
### Leonard Mada
### [the one and only]
###
### Integrals:
### Moebius Transforms
###
### draft v.0.1a


### Integrals: Moebius Transforms

### Tr = Moebius Transform
# Tr((lower, upper)) = (upper, lower);
# f(Tr(x)) = k * f(x), where k = constant;

### Usage:
# - examples: P(Tr(x)) = - P(x);
# - f(x) = sin(P[odd](x))
#    => f(P(Tr(x)))  = - f(P(x));
#    => I( f(P(x)) ) = - I( f(P(x)) ) => I = 0;
# - f(x) = acos(P[odd](x))
#    => f(P(Tr(x))) = pi - f(P(x));
#    => I( f(P(x)) ) = pi/2;


##########################

###############
### History ###
###############


### draft v.0.1a:
# - moved to this file;
### [old]
# - from file: Integrals.Tricks.R;
### draft v.0.3c - v.0.3c-ex2:
# - definite integrals using Moebius Transform:
#   I( log(4*x-3) / (x*(3*x+1)) ) on c(1, 4);
#   I( log(5*x-3) / (x*(x+1)) ) on c(1, 3);
#   I( log(x-3) / (x*(3*x+1)) ) on c(5, 8);
### draft v.0.3a - v.0.3b:
# - Transformations:
#   f((a*x + b) / (c*x + d)) = - f(x);
# - added 2nd example;
# - acos() variant; [v.0.3b]


############################

######################
### Decompositions ###
######################

### Moebius Transform
### Tr(x) = (a*x + b) / (c*x + d)

### f(Tr(x)) = k * f(x)
# - where k = constant, e.g. k = -1;

### Examples:

### Ex 1:
# Tr(x) = (3*x+3) / (5*x - 3)
# Tr(c(1, 3)) = (3, 1);
# f1(x) = f(x) / (5*x - 3);
f1 = function(x) 1/(5*x-3) * sin(((1-2*sqrt(6)/3)*x+1) / ((1+2*sqrt(6)/3)*x + 1))
integrate(f1, lower=1, upper=3)
# I(f1)[1 -> 3] = 0


f2 = function(x) x/(5*x-3) * sin(((1-2*sqrt(6)/3)*x+1) / ((1+2*sqrt(6)/3)*x + 1))
f3 = function(x) 1/5 * sin(((1-2*sqrt(6)/3)*x+1) / ((1+2*sqrt(6)/3)*x + 1))
integrate(f2, lower=1, upper=3)
integrate(f3, lower=1, upper=3)

### relation between sin() & cos()
alpha = (1-2*sqrt(6)/3) / (1+2*sqrt(6)/3);
f1 = function(x) 1/(5*x-3) * sin((1 - alpha) / ((1+2*sqrt(6)/3)*x + 1))
f2 = function(x) 1/(5*x-3) * cos((1 - alpha) / ((1+2*sqrt(6)/3)*x + 1))
integrate(f1, lower=1, upper=3)
integrate(f2, lower=1, upper=3)$value * (- tan(alpha))


### Ex 1b:
f1 = function(x) 1/(5*x-3) * acos(((1-2*sqrt(6)/3)*x+1) / ((1+2*sqrt(6)/3)*x + 1))
f1.exact = function(x) pi/10 * log(5*x-3)
integrate(f1, lower=1, upper=3)
f1.exact(3) - f1.exact(1)

### Ex 1c:
f1 = function(x) log(5*x-3) / (x*(x+1));
f1.exact = function(x) 2 * log(5*x-3) *
	acos(((1-sqrt(8/3))*x+1) / ((1+sqrt(8/3))*x + 1)) +
	- pi * log(5*x-3);
integrate(f1, lower=1, upper=3)
f1.exact(3) - f1.exact(1)
# small numeric instability?

curve(f1(x), from=1, to=3)


#########

#########
### Ex 2:
# Tr(x) = (3*x+1) / (4*x - 3)
# Tr(c(1, 4)) = (4, 1);
# f1(x) = f(x) / (4*x - 3);
f1 = function(x) 1/(4*x-3) * sin(((3-sqrt(13))*x+1) / ((3+sqrt(13))*x + 1))
integrate(f1, lower=1, upper=4)


### Ex 2b:
f1 = function(x) 1/(4*x-3) * acos(((3-sqrt(13))*x+1) / ((3+sqrt(13))*x + 1))
f1.exact = function(x) pi/8 * log(4*x-3)
integrate(f1, lower=1, upper=4)
f1.exact(4) - f1.exact(1)


###
f1 = function(x) log(4*x-3) / (x*(3*x+1));
f1.exact = function(x) 2 * log(4*x-3) *
	acos(((3-sqrt(13))*x+1) / ((3+sqrt(13))*x + 1)) +
	- pi * log(4*x-3);
integrate(f1, lower=1, upper=4)
f1.exact(4) - f1.exact(1)
# small numeric instability?


#########

#########
### Ex 3:
# Tr(x) = (3*x+1) / (x - 3)
# Tr(c(5, 8)) = (8, 5);
# f1(x) = f(x) / (x - 3);
f1 = function(x) 1/(x-3) * sin(((3-sqrt(10))*x+1) / ((3+sqrt(10))*x + 1))
integrate(f1, lower=5, upper=8)


### Ex 3b:
f1 = function(x) 1/(x-3) * acos(((3-sqrt(10))*x+1) / ((3+sqrt(10))*x + 1))
f1.exact = function(x) pi/2 * log(x-3)
integrate(f1, lower=5, upper=8)
f1.exact(8) - f1.exact(5)


###
f1 = function(x) log(x-3) / (x*(3*x+1));
f1.exact = function(x) 2 * log(x-3) *
	acos(((3-sqrt(10))*x+1) / ((3+sqrt(10))*x + 1)) +
	- pi * log(x-3);
integrate(f1, lower=5, upper=8)
f1.exact(8) - f1.exact(5)
# small numeric instability?

