########################
##
## Leonard Mada
## [the one and only]
##
## Exact Integration:
##   Polynomial Fractions:
##   Mixed Radicals: 2 Radicals
##
## draft v.0.1a


### Mixed Radicals: 2 Radicals



### I( x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1) )
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1), 0, 1)
1/4 * (beta(3/4, 2-3/4) - beta(1/4, 3/4)) + pi / tan(pi/6) / 6;
1/4 * (- 3/4 * beta(1/4, 3/4)) + pi / tan(pi/6) / 6;


# Derivation:
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 * (1 - x^4)^(1/2) + 1), 0, 1)
integrate(\(x) 1/4 * (x*(1-x)^3)^(1/4) / ((x*(1-x))^(1/2) + 1), 0, 1)
integrate(\(x) x^6 / (x^4 + 1)^3 / (x^2/(x^4+1) + 1), 0, Inf)
integrate(\(x) x^4 / (x^4 + 1)^3 / (x^2/(x^4+1) + 1), 0, Inf)
integrate(\(x) x^4 / (x^4 + 1)^2 / (x^4 + x^2 + 1), 0, Inf)
integrate(\(x) x^2 / (x^4 + 1)^2 - 1/(x^4+1) + 1/(x^4 + x^2 + 1), 0, Inf)
1/4 * (beta(3/4, 2-3/4) - beta(1/4, 3/4)) + pi / tan(pi/6) / 6;


#
integrate(\(x) 1/(x^4 + x^2 + 1), 0, Inf)
- pi / tan(3*pi/6) / 6 + pi / tan(pi/6) / 6;
pi / tan(pi/6) / 6


##################

### I( x^4 * (1 - x^4)^(5/4) / (x^2 + (1 - x^4)^(1/2)) )
integrate(\(x) x^4 * (1 - x^4)^(5/4) / (x^2 + (1 - x^4)^(1/2)), 0, 1)
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 / (1 - x^4)^(1/2) + 1), 0, 1)
pi/16 - beta(1/4, 3/4) / 32;

### I( x^4 * (1 - x^4)^(5/4) * (x^2 + (1 - x^4)^(1/2)) )
integrate(\(x) x^4 * (1 - x^4)^(5/4) * (x^2 + (1 - x^4)^(1/2)), 0, 1)
3/128 * beta(1/4, 3/4)


# Derivation:
integrate(\(x) x^4 * (1 - x^4)^(3/4) / (x^2 / (1 - x^4)^(1/2) + 1), 0, 1)
integrate(\(x) 1/4 * (x*(1-x)^3)^(1/4) / ((x/(1-x))^(1/2) + 1), 0, 1)
integrate(\(x) x^6 / (x^4 + 1)^3 / (1/x^2 + 1), 0, Inf)
integrate(\(x) x^8 / (x^4 + 1)^3 / (x^2 + 1), 0, Inf)
integrate(\(x) 1/2 * x^8 / (x^4 + 1)^2 *
	(1/(x^2 + 1) - (x^2-1)/(x^4 + 1)), 0, Inf)
integrate(\(x) 1/2 * x^8 / (x^4 + 1)^2 / (x^2 + 1) +
	- 1/2 * (x^2-1) / (x^4 + 1) + (x^2-1) / (x^4 + 1)^2 +
	- 1/2 * (x^2-1) / (x^4 + 1)^3, 0, Inf)
integrate(\(x) 1/4 * x^8 / (x^4 + 1) / (x^2 + 1) +
	- 1/4 * x^8 * (x^2-1) / (x^4 + 1)^2 +
	- 1/2 * (x^2-1) / (x^4 + 1) + (x^2-1) / (x^4 + 1)^2 +
	- 1/2 * (x^2-1) / (x^4 + 1)^3, 0, Inf)
integrate(\(x) 1/8 / (x^2 + 1) +
	- 1/4 * (x^2-1) / (x^4 + 1) +
	+ 3/4 * (x^2-1) / (x^4 + 1)^2 +
	- 1/2 * (x^2-1) / (x^4 + 1)^3, 0, Inf)
pi/16  - beta(3/4, 1-3/4)/16 + beta(1/4, 1-1/4)/16 +
	+ beta(3/4, 2-3/4) * 3/16 - beta(1/4, 2-1/4) * 3/16 +
	- beta(3/4, 3-3/4) / 8 + beta(1/4, 3-1/4) / 8;

# Derivation:
integrate(\(x) x^4 * (1 - x^4)^(5/4) * (x^2 + (1 - x^4)^(1/2)), 0, 1)
# integrate(\(x) x^4 * (1 - x^4)^(5/4) / ((1 - x^4)^(1/2) - x^2), 0, 1)
integrate(\(x) x^8 / (x^4 + 1)^3 / (x^2 - 1), 0, 1 - 1E-4)$value +
integrate(\(x) x^8 / (x^4 + 1)^3 / (x^2 - 1), 1+ 1E-4, Inf)$value;
integrate(\(x) (x^6 + x^4 + 1/2*x^2 + 1/2) / (x^4 + 1)^3 +
	- 1/4 * (x^2 + 1) / (x^4 + 1)^2 +
	- 1/8 * (x^2 + 1) / (x^4 + 1), 0, Inf); # + pi/8 / atan(pi/2);
beta(7/4, 3-7/4)/4 + beta(5/4, 3-5/4)/4 +
	+ 1/8 * beta(3/4, 3-3/4) + 1/8 * beta(1/4, 3-1/4) +
	- 1/16 * beta(3/4, 2-3/4) - 1/16 * beta(1/4, 2-1/4) +
	- 1/32 * beta(3/4, 1/4) - 1/32 * beta(1/4, 3/4);
3/128 * beta(1/4, 3/4)

