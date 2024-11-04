

### ATAN( Trig )

### Examples:
# I( atan(2*tan(x)) )
# I( atan(k*sin(x)) )


###################
### Simple Argument

### I( atan(2*tan(x)) )
integrate(\(x) atan(2*tan(x)), 0, pi/2)
5/24 * pi^2 + pracma::polylog(-1/2, 2)/2 + (log(3/2)^2 - log(3)^2)/4;


# Derivation:
# - without the Lerch zeta function;
integrate(\(x) atan(2*tan(x)), 0, pi/2)
pi^2/4 - integrate(\(x) 2 * atan(x) / (4*x^2 + 1), 0, Inf)$value
integrate(\(x) 2 * atan(x) / (x^2 + 4), 0, Inf)
integrate(\(x) log((1+x)/(1-x)) / (1 + 2*x), 0, 1)$value +
integrate(\(x) log((1+x)/(x-1)) / (1 + 2*x), 1, Inf)$value;
5/24 * pi^2 + pracma::polylog(-1/2, 2)/2 + (log(3/2)^2 - log(3)^2)/4;

#
integrate(\(x) log((1+x)/x) / (1 + 2*x), 0, Inf)
pi^2 / 8;

#
integrate(\(x) log((1-x)/x) / (1 + 2*x), 0, 1)$value +
integrate(\(x) log(abs(1-x)/x) / (1 + 2*x), 1, Inf)$value;
- pracma::polylog(-1/2, 2)/2 - pi^2/12 - (log(3/2)^2 - log(3)^2)/4;


######################
######################

### I( atan(k*sin(x)) )
# Maths 505: A surprisingly difficult integral:
#  int 0 to π/2 arctan(2sin(x)) solution using Feynman's trick
# https://www.youtube.com/watch?v=iawTwy_Gvkw

# library(pracma)
# Note:
# - function polylog2 is currently in file:
#   Integrals.Trig.Tan.R;

integrate(\(x) atan(2*sin(x)), 0, pi/2)
polylog2(2/(sqrt(5)+1), 2) - polylog2(-2/(sqrt(5)+1), 2)

### I( atan(k*sin(x)) )
k = sqrt(3);
kk = k / sqrt(k^2 + 1);
integrate(\(x) atan(k*sin(x)), 0, pi/2)
polylog2(tan(asin(kk)/2), 2) - polylog2(- tan(asin(kk)/2), 2)
polylog2(tan(atan(k)/2), 2) - polylog2(- tan(atan(k)/2), 2)


#####################

### I( atan(a*sin(x)) / sin(x) )
# Maths 505: AN INCREDIBLE CALCULUS RESULT: solution using Feynman's technique
# https://www.youtube.com/watch?v=WLq2EThghgc

###
a = sqrt(3)
integrate(\(x) atan(a*sin(x)) / sin(x), 0, pi/2)
pi/2 * asinh(a)


### I( atan(2*cos(x)^2) / cos(x)^2 )
# Maths 505: A RIDICULOUSLY AWESOME INTEGRAL FEAT. φ & π
# https://www.youtube.com/watch?v=6CJ_Sal55uM

integrate(\(x) atan(2*cos(x)^2) / cos(x)^2, 0, pi/2)
integrate(\(x) 4 * sin(x)^2 / (4*cos(x)^4 + 1), 0, pi/2)
integrate(\(x) 4 / (4*x^4 + (x^2+1)^2), 0, Inf)
pi / sqrt((sqrt(5) + 1) / 2)


### I( atan(3*cos(x)^2) / cos(x)^2 )
integrate(\(x) atan(3*cos(x)^2) / cos(x)^2, 0, pi/2)
integrate(\(x) 6 * sin(x)^2 / (9*cos(x)^4 + 1), 0, pi/2)
integrate(\(x) 6 / (9*x^4 + (x^2+1)^2), 0, Inf)
integrate(\(x) 6*x^2 / (9 + (x^2+1)^2), 0, Inf)
integrate(\(x) x^2 * Im(1/(x^2 + 1-3i) - 1/(x^2 + 1+3i)), 0, Inf)
integrate(\(x) Im((1+3i)/(x^2 + 1+3i) - (1-3i)/(x^2 + 1-3i)), 0, Inf)
integrate(\(x) Im((1+3i)^(1/2)/(x^2 + 1) - (1-3i)^(1/2)/(x^2 + 1)), 0, Inf)
pi/2 * Im((1+3i)^(1/2) - (1-3i)^(1/2))


# Alternatives:
integrate(\(x) (1 - cos(x)) / (cos(x)^2 + 2*cos(x) + 2), 0, pi)
# Alt: Weierstrass;
integrate(\(x) (1 - cos(x)) / (cos(x)^2 + 2*cos(x) + 2) +
	(1 + cos(x)) / (cos(x)^2 - 2*cos(x) + 2), 0, pi/2)
integrate(\(x) (6*cos(x)^2 + 4) / (cos(x)^4 + 4), 0, pi/2)
pi / sqrt((sqrt(5) + 1) / 2)


#######################
#######################

### I( atan(tan(x)^2) )
# Maths 505: An outrageous journey of integration: int 0 to pi/4 arctan(cot^2(x))
# https://www.youtube.com/watch?v=VUGlU_dSgPY

# Note:
# - additional variants in file with ATAN-Fractions;

### I( atan(x^2) / (x^2 + 1) )
integrate(\(x) atan(tan(x)^2), 0, pi/4)
integrate(\(x) atan(x^2) / (x^2 + 1), 0, 1)
pi^2/16 - (digamma(3/8 + 1/2) - digamma(3/8)) *
	(digamma(1/8 + 1/2) - digamma(1/8)) / 32;

