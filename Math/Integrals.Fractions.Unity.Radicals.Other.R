#####################
##
## Leonard Mada
## [the one and only]
##
## Polynomial Radicals:
## Other Types


### Note:
# - Basic types are in file:
#   Integrals.Fractions.Unity.Radicals.R;


##########################
##########################

### I( sqrt(sqrt(1 + x^4) + 1) )
# Maths 505: A deceivingly tough integral
# https://www.youtube.com/watch?v=nGx8j7-KHxI
# Subst: x = sqrt(tan(x));

### I( sqrt(sqrt(1 + x^4) + 1) )
integrate(\(x) sqrt(sqrt(x^4 + 1) + 1), 0, 1)
integrate(\(x) (1 + x^4) / (1 - x^4)^2 * 2, 0, sqrt(tan(pi/8)))
integrate(\(x) 1/(1-x^2)^2 + 1/(1+x^2)^2, 0, sqrt(tan(pi/8)))
z = sqrt(tan(pi/8))
z/(1-z^4) + atan(z)/2 + log((1+z)/(1-z)) / 4
z/(1-z^4) + atan(z)/2 + atan(1i*z) / 2i


### I( sqrt(sqrt(1 + x^4) - 1) )
integrate(\(x) sqrt(sqrt(x^4 + 1) - 1), 0, 1)
integrate(\(x) sqrt((1-cos(x))/sin(x)) / cos(x)^2 / 2, 0, pi/4)
integrate(\(x) sqrt(sin(x/2)/cos(x/2)) / cos(x)^2 / 2, 0, pi/4)
integrate(\(x) 2 * x^2 * (1 + x^4)/(1-x^4)^2, 0, sqrt(tan(pi/8)))
integrate(\(x) x^2/(1-x^2)^2 + 1/(1+x^2) - 1/(1+x^2)^2, 0, sqrt(tan(pi/8)))
z = sqrt(tan(pi/8))
z^3/(1-z^4) + atan(z)/2 - atan(1i*z) / 2i


###
integrate(\(x) sqrt(sqrt(x^4 + 1) + x^2), 0, 1)
z = sqrt(tan(pi/8))
(z/(1-z^2) + atan(z)) / sqrt(2)

###
integrate(\(x) sqrt(sqrt(x^4 + 1) - x^2), 0, 1)
z = sqrt(tan(pi/8))
(z/(1+z^2) + atan(1i*z) / 1i) / sqrt(2)


### Gen:

###
b = sqrt(3);
integrate(\(x) sqrt(sqrt(x^4 + b^4) + b^2), 0, 1)
integrate(\(x) sqrt(sqrt(x^4 + 1) + 1) * b^2, 0, 1/b)
z = sqrt(tan(atan(1/b^2)/2))
b^2 * (z/(1-z^4) + atan(z)/2 + atan(1i*z) / 2i)

###
integrate(\(x) sqrt(sqrt(x^4 + b^4) - b^2), 0, 1)
z = sqrt(tan(atan(1/b^2)/2))
b^2 * (z^3/(1-z^4) + atan(z)/2 - atan(1i*z) / 2i);


##############

### Pow 2/3

### I( (sqrt(x^6 + 1) + 1)^(2/3) )
integrate(\(x) (sqrt(x^6 + 1) + 1)^(2/3), 0, 1)
integrate(\(x) 1/3 * (cos(x/2)/sin(x/2))^(2/3) / cos(x)^2, 0, pi/4)
integrate(\(x) 2/3 * x^(-2/3) * (1 + 4*x^2/(1-x^2)^2) / (1+x^2), 0, tan(pi/8))
integrate(\(x) 2/3 * x^(-2/3) * (1+x^2) / (1-x^2)^2, 0, tan(pi/8))
integrate(\(x) 2 * (1+x^6) / (1-x^6)^2, 0, (tan(pi/8))^(1/3))
z = (tan(pi/8))^(1/3);
integrate(\(x) 2/3 /(1-x^3) + 2/3 /(1+x^3), 0, z)$value +
	+ 2/3 * z / (1-z^6);
1/9 * (log(z^2+z+1) + 2*sqrt(3)*atan((2*z+1)/sqrt(3)) - 4i*atan(1i*z) +
	- log(z^2-z+1) + 2*sqrt(3)*atan((2*z-1)/sqrt(3)) ) +
	+ 2/3 * z / (1-z^6);
1/9 * (log(z^2+z+1) - log(z^2-z+1) +
	+ 2*sqrt(3) * atan(sqrt(3)*z / (1 - z^2)) - 4i*atan(1i*z) ) +
	+ 2/3 * z / (1-z^6);


### I( (sqrt(x^6 + 1) - 1)^(2/3) )
integrate(\(x) (sqrt(x^6 + 1) - 1)^(2/3), 0, 1)
integrate(\(x) 1/3 * (sin(x/2)/cos(x/2))^(2/3) / cos(x)^2, 0, pi/4)
integrate(\(x) 2/3 * x^(2/3) * (1 + 4*x^2/(1-x^2)^2) / (1+x^2), 0, tan(pi/8))
integrate(\(x) 2/3 * x^(2/3) * (1+x^2) / (1-x^2)^2, 0, tan(pi/8))
integrate(\(x) 2 * x^4 * (1+x^6) / (1-x^6)^2, 0, (tan(pi/8))^(1/3))
z = (tan(pi/8))^(1/3);
integrate(\(x) 2/3 * x/(1+x^3) - 2/3 * x/(1-x^3), 0, up)$value +
	+ diff(x^2/(1-x^3) - x^2/(1+x^3)) / 3;
# TODO


### Gen:
b = sqrt(5)
integrate(\(x) (sqrt(x^6 + b^6) + b^3)^(2/3), 0, 1)
z = (tan(atan(1/b^3)/2))^(1/3);
b^3 / 9 * (log(z^2+z+1) - log(z^2-z+1) + 6*z / (1-z^6) +
	+ 2*sqrt(3) * atan(sqrt(3)*z / (1 - z^2)) - 4i*atan(1i*z) );


### Derivation:

# Fraction Decomposition
up = (tan(pi/8))^(1/3);
integrate(\(x) 2 * (1+x^6) / (1-x^6)^2, 0, up)
integrate(\(x) -2/(1-x^6) + 4 / (1-x^6)^2, 0, up)
integrate(\(x) -2/(1-x^6) + (1/(1-x^3) + 1/(1+x^3))^2, 0, up)
integrate(\(x) 1/(1-x^3)^2 + 1/(1+x^3)^2, 0, up)
# =>
x = c(0,up);
integrate(\(x) 2/3 /(1-x^3) + 2/3 /(1+x^3), 0, up)$value +
	+ diff(x/(1-x^3) + x/(1+x^3)) / 3;

# Fraction 2:
up = (tan(pi/8))^(1/3);
integrate(\(x) 2 * x^4 * (1+x^6) / (1-x^6)^2, 0, up)
integrate(\(x) x^4/(1-x^3)^2 + x^4/(1+x^3)^2, 0, up)
integrate(\(x) x/(1+x^3) - x/(1-x^3) + x/(1-x^3)^2 - x/(1+x^3)^2, 0, up)
# =>
x = c(0,up);
integrate(\(x) 2/3 * x/(1+x^3) - 2/3 * x/(1-x^3), 0, up)$value +
	+ diff(x^2/(1-x^3) - x^2/(1+x^3)) / 3;


# Integration by Parts:
x = c(0,up);
integrate(\(x) 1/(1-x^3), 0, up)
diff(x/(1-x^3)) - integrate(\(x) 3*x^3/(1-x^3)^2, 0, up)$value
diff(x/(1-x^3)) + integrate(\(x) 3/(1-x^3) - 3/(1-x^3)^2, 0, up)$value

#
integrate(\(x) 1/(1+x^3), 0, up)
diff(x/(1+x^3)) + integrate(\(x) 3*x^3/(1+x^3)^2, 0, up)$value
diff(x/(1+x^3)) + integrate(\(x) 3/(1+x^3) - 3/(1+x^3)^2, 0, up)$value

