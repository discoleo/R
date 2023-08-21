

#####################

### on [0,1]

### I( x * log(x+1) / (x+1)^s )
s = 1/sqrt(3)
integrate(\(x) x * log(x+1) / (x+1)^s, 0, 1)
2^(2-s)* log(2)/(2-s) - 2^(1-s)*log(2)/(1-s) +
	- (2^(2-s)-1)/(2-s)^2 + (2^(1-s)-1)/(1-s)^2;


### I( x^2 * log(x+1) / (x+1)^s )
s = 1/sqrt(3)
integrate(\(x) x^2 * log(x+1) / (x+1)^s, 0, 1)
(4/(3-s) - 4/(2-s) + 1/(1-s)) * 2^(1-s)*log(2) +
	- (2^(3-s)-1)/(3-s)^2 + (2^(3-s)-2)/(2-s)^2 - (2^(1-s)-1)/(1-s)^2;


#############

### Power 1/2

### I( log(x) * sqrt(x + 1) )
integrate(\(x) log(x) * sqrt(x + 1), 0, 1)
- 2/3 * (2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1))
4/3 * (asinh(1) - log(2) - (1 + 2/3)*sqrt(2) + 1 + 1/3)

# Derivation:
integrate(\(x) (sqrt(x + 1) - 1) / x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2
#
integrate(\(x) ((x+1)*sqrt(x + 1) - 1) / x, 0, 1)
2*log(2) - 2*asinh(1) + 2*sqrt(2) - 2 + 2/3 * (2^(3/2) - 1)


### I( log(x) / sqrt(x + 1) )
integrate(\(x) log(x) / sqrt(x + 1), 0, 1)
4 * (asinh(1) - log(2) - sqrt(2) + 1)


### Power 1/3

### I( log(x) * (x + 1)^(1/3) )
integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) + 9/8*log(2^(2/3) + 2^(1/3) + 1) - 9/8*log(3) +
	+ 9/8 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	- 9/8 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));

### I( log(x) * (x + 1)^(2/3) )
integrate(\(x) log(x) * (x + 1)^(2/3), 0, 1)
integrate(\(x) - 3/5 * ((x + 1)^(5/3) - 1) / x, 0, 1)
-9*(2^(5/3)/25 + 2^(2/3)/10) + 9*(1/25 + 1/10) +
	+ 9/10 * log(2^(2/3) + 2^(1/3) + 1) - 9/10*log(3) +
	- 9/10 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	+ 9/10 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));


### I( log(x) / (x + 1)^(1/3) )
integrate(\(x) log(x) / (x + 1)^(1/3), 0, 1)
5/2 * integrate(\(x) log(x) * (x + 1)^(2/3), 0, 1)$value +
	+ 3*2*2^(2/3)*(1/2 - 1/5) - 3*(1/2 - 1/5);
# see above for I( log(x) * (x+1)^(2/3) );


### I( log(x) / (x + 1)^(2/3) )
integrate(\(x) log(x) / (x + 1)^(2/3), 0, 1)
4 * integrate(\(x) log(x) * (x + 1)^(1/3), 0, 1)$value +
	+ 3*2*2^(1/3)*(1 - 1/4) - 3*(1 - 1/4);
# see above for I( log(x) * (x+1)^(1/3) );

# Derivations:
integrate(\(x) - 3/4 * ((x+1)^(4/3) - 1) / x, 0, 1)
integrate(\(x) - 9/4 * (x^6 - x^2) / (x^3 - 1), 1, 2^(1/3))
-9/16*(2^(4/3) - 1) + integrate(\(x) - 9/4 * (x^3 - x^2) / (x^3 - 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (x + 1) / (x^2 + x + 1), 1, 2^(1/3))$value
-9/16*(2^(4/3) - 1 + 4*2^(1/3) - 4) +
	+ integrate(\(x) 9/4 * (1/2*(2*x+1) + 1/2) / (x^2 + x + 1), 1, 2^(1/3))$value;
	# + 9/8 / sin(pi/3) * atan((x + cos(pi/3)) / sin(pi/3))
	# + integrate(\(x) 9/8 / (x^2 + x + 1), 1, 2^(1/3))$value

#
integrate(\(x) - 3/5 * ((x + 1)^(5/3) - 1) / x, 0, 1)
integrate(\(x) - 9/5 * (x^7 - x^2)/ (x^3 - 1), 1, 2^(1/3))
-9/25 * (2^(5/3) - 1) - 9/10*(2^(2/3) - 1) +
	+ integrate(\(x) 9/5 * x / (x^2 + x + 1), 1, 2^(1/3))$value
-9*(2^(5/3)/25 + 2^(2/3)/10) + 9*(1/25 + 1/10) +
	+ 9/10 * (log(2^(2/3) + 2^(1/3) + 1) - log(3)) +
	- 9/10 / sin(pi/3) * atan((2^(1/3) + cos(pi/3)) / sin(pi/3)) +
	+ 9/10 / sin(pi/3) * atan((1 + cos(pi/3)) / sin(pi/3));
	# + integrate(\(x) 9/5 * x / (x^2 + x + 1), 1, 2^(1/3))$value


# Derivation:
integrate(\(x) ((x+1)^(1/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^3 - x^2) / (x^3 - 1), 1, 2^(1/3))
#
integrate(\(x) ((x+1)^(4/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^6 - x^2) / (x^3 - 1), 1, 2^(1/3))

#
integrate(\(x) ((x+1)^(2/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^4 - x^2) / (x^3 - 1), 1, 2^(1/3))
#
integrate(\(x) ((x+1)^(5/3) - 1)/x, 0, 1)
integrate(\(x) 3*(x^7 - x^2) / (x^3 - 1), 1, 2^(1/3))


###################

### I( x^p * log(x) / (x + 1)^s ) on [0, 1]
# s = Integer

# TODO: truncated Polylog function;

###
integrate(function(x) log(x) / (x+1)^2, 0, 1)
- log(2)

###
integrate(function(x) log(x) / (x+1)^3, 0, 1)
- log(2)/2 + (1/2 - 1)/2

###
integrate(function(x) log(x) / (x+1)^4, 0, 1)
- log(2)/3 + (1/2 + 1/4/2 - 1 - 1/2)/3

###
integrate(function(x) log(x) / (x+1)^5, 0, 1)
- log(2)/4 + (1/2 + 1/4/2 + 1/8/3 - 1 - 1/2 - 1/3)/4

###
integrate(function(x) log(x) / (x+1)^6, 0, 1)
- log(2)/5 + (1/2 + 1/4/2 + 1/8/3 + 1/16/4 - 1 - 1/2 - 1/3 - 1/4)/5


