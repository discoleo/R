

########################
########################

########################
### Simple Fractions ###
########################

### I( log(x) / (x + 1)^p ) on [0, Inf]

# qncubed3: Complex Analysis: Integral of log(x)/(x+1)^2
# https://www.youtube.com/watch?v=tPveHNdBWR8

# Note:
# - easy pole with higher multiplicity;

# TODO: full generalization;

integrate(function(x) log(x)/(x+1)^2, 0, Inf)
# == 0

###
integrate(function(x) log(x)/(x+1)^3, 0, Inf)
-1/2

###
integrate(function(x) log(x)/(x+1)^4, 0, Inf)
-1/2

### n = 5
integrate(function(x) log(x)/(x+1)^5, 0, Inf)
- 11 / gamma(5)

### n = 6
integrate(function(x) log(x)/(x+1)^6, 0, Inf)
- 50 / gamma(6)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(- 100/x^5 + 48*log(x)/x^5)/gamma(6)


### n = 7
integrate(function(x) log(x)/(x+1)^7, 0, Inf)
- 274 / gamma(7)

# x = exp(1i*pi);
# - 4i*pi*I + 4*pi^2/(n-1) = 2i*pi*(2*274 - 240*log(x)/x^6)/gamma(7)
# - 2i*I + 2*pi/(n-1) = 1i*(2*274 - 240*log(x)/x^6)/gamma(7)


########################
########################

### I( log(x) / (x + 1)^p ) on [0, 1]

# TODO: truncated Polylog function;

###
integrate(function(x) log(x)/(x+1)^2, 0, 1)
- log(2)

###
integrate(function(x) log(x)/(x+1)^3, 0, 1)
- log(2)/2 + (1/2 - 1)/2

###
integrate(function(x) log(x)/(x+1)^4, 0, 1)
- log(2)/3 + (1/2 + 1/4/2 - 1 - 1/2)/3

###
integrate(function(x) log(x)/(x+1)^5, 0, 1)
- log(2)/4 + (1/2 + 1/4/2 + 1/8/3 - 1 - 1/2 - 1/3)/4

###
integrate(function(x) log(x)/(x+1)^6, 0, 1)
- log(2)/5 + (1/2 + 1/4/2 + 1/8/3 + 1/16/4 - 1 - 1/2 - 1/3 - 1/4)/5


### I( log(x) * sqrt(x + 1) )
integrate(\(x) log(x) * sqrt(x + 1), 0, 1)
- 2/3 * (2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1))
4/3 * (asinh(1) - log(2) - (1 + 2/3)*sqrt(2) + 1 + 1/3)

# Derivation:
integrate(\(x) (sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2
#
integrate(\(x) ((x+1)*sqrt(x + 1) - 1)/x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1)


### I( log(x) / sqrt(x + 1) )
integrate(\(x) log(x) / sqrt(x + 1), 0, 1)
4 * (asinh(1) - log(2) - sqrt(2) + 1)


### I( log(x) * (x + 1)^(1/3) )
integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) + 9/8*log(2^(2/3) + 2^(1/3) + 1) - 9/8*log(3) +
	+ 9/8 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	- 9/8 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));


### I( log(x) / (x + 1)^(2/3) )
integrate(\(x) log(x) / (x + 1)^(2/3), 0, 1)
4 * integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)$value +
	+ 3*2*2^(1/3)*(1 - 1/4) - 3*(1 - 1/4);
# see above for I( log(x) * (x+1)^(1/3) );

#
integrate(\(x) - 3/4 * ((x+1)^(4/3) - 1)/x, 0, 1)
integrate(\(x) - 9/4 * (x^6 - x^2) / (x^3 - 1), 1, 2^(1/3))
-9/16*(2^(4/3) - 1) + integrate(\(x) - 9/4 * (x^3 - x^2) / (x^3 - 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (x + 1) / (x^2 + x + 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (1/2*(2*x+1) + 1/2) / (x^2 + x + 1), 1, 2^(1/3))$value;
	# + 9/8 / sin(pi/3) * atan((x + cos(pi/3)) / sin(pi/3))
	# + integrate(\(x) 9/8 / (x^2 + x + 1), 1, 2^(1/3))$value


# Derivation:
integrate(\(x) ((x+1)^(1/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^3 - x^2)/ (x^3 - 1), 1, 2^(1/3))
#
integrate(\(x) ((x+1)^(4/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^6 - x^2)/ (x^3 - 1), 1, 2^(1/3))


