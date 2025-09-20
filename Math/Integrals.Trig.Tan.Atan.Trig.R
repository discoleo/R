

### ATAN( Trig )

### Examples:
# I( atan(2*tan(x)) )
# I( atan(k*sin(x)) )


###################
### Simple Argument

### I( atan(2*tan(x)) )
integrate(\(x) atan(2*tan(x)), 0, pi/2)
5/24 * pi^2 + pracma::polylog(-1/2, 2)/2 + (log(3/2)^2 - log(3)^2)/4;

### I( atan(tan(x)/2) )
integrate(\(x) atan(tan(x)/2), 0, pi/2)
pi^2/4 - integrate(\(x) 2 * atan(x) / (x^2 + 4), 0, Inf)$value
1/24 * pi^2 - pracma::polylog(-1/2, 2)/2 - (log(3/2)^2 - log(3)^2)/4


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


#####################

### I( atan(3*tan(x)) )
integrate(\(x) atan(3*tan(x)), 0, pi/2)
pi^2/4 - integrate(\(x) 3 * atan(x) / (9*x^2 + 1), 0, Inf)$value
integrate(\(x) 3 * atan(x) / (x^2 + 9), 0, Inf)
pi^2/4 - log(3)*log(2)/2 +
	- (pracma::polylog(1/4, 2) - pracma::polylog(-1/2, 2)) / 2 +
	- (log(4/3)^2 - log(3/2)^2) / 4;


### I( atan(tan(x)/3) )
integrate(\(x) atan(tan(x)/3), 0, pi/2)
pi^2 / 4 - integrate(\(x) atan(3*tan(x)), 0, pi/2)$value
log(3)*log(2)/2 +
	+ (pracma::polylog(1/4, 2) - pracma::polylog(-1/2, 2)) / 2 +
	+ (log(4/3)^2 - log(3/2)^2) / 4;


# Derivation:

# from atan(x/b) / (x^2 + 9)
b = 5^(1/3)
integrate(\(x) x / (x^2 + b^2) / (x^2 + 9), 0, Inf)
integrate(\(x) x * (1/(x^2 + 9) - 1/(x^2 + b^2)) / (b^2 - 9), 0, Inf)
log(b/3) / (b^2 - 9)

# I( ... db )
integrate(\(x) atan(x) / (x^2 + 9), 0, Inf)
pi^2/12 - log(3)*log(2)/6 +
	- (pracma::polylog(1/4, 2) - pracma::polylog(-1/2, 2)) / 6 +
	- (log(4/3)^2 - log(3/2)^2) / 12

#
integrate(\(b) log(b) / (b^2 - 9), 0, 1)
(pracma::polylog(1/4, 2) - pracma::polylog(-1/2, 2)) / 6 +
	+ (log(4/3)^2 - log(3/2)^2) / 12;


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
# 1. Maths 505: AN INCREDIBLE CALCULUS RESULT: solution using Feynman's technique
# https://www.youtube.com/watch?v=WLq2EThghgc
# Note: Feynman trick;
# 2. Michael Penn: an integral with some of my favorite techniques
# https://www.youtube.com/watch?v=CWPX9oCnd1I
# Note: as double integral, then classic;

###
k = sqrt(3)
integrate(\(x) atan(k*sin(x)) / sin(x), 0, pi/2)
pi/2 * asinh(k)


### I( atan(cos(x)) / sin(x) )
integrate(\(x) (pi/4 - atan(cos(x))) / sin(x), 0, pi/2)
Catalan / 2;


### I( atan(sin(x)) / sin(x)^2 )
integrate(\(x) atan(sin(x)) / sin(x)^2 - 1/x, 0, pi/2, rel.tol=1E-13)
- sqrt(2) * log(sqrt(2)+1) + log(2) + 1 - log(pi/2);

### I( atan(sin(x)^2) / sin(x)^2 )
integrate(\(x) atan(sin(x)^2) / sin(x)^2, 0, pi/2, rel.tol=1E-13)
pi * 2^(3/4) / 4 / cos(pi/8);


### I( sin(2*x) * atan(b*sin(2*x)) )
# Maths 505: A captivating integral problem
# https://www.youtube.com/watch?v=hqC_wLN7UYo
# Subst: x => sin(y)^2 & Feynman;

integrate(\(x) atan(sqrt(x*(1-x))), 0, 1)
pi/2*(sqrt(5) - 2)

###
b = 1/3
integrate(\(x) sin(2*x) * atan(b*sin(2*x)), 0, pi/2)
pi/2*(sqrt(b^2+1) - 1) / b

###
b = 1/3
integrate(\(x) sin(2*x) * atan(b/sin(2*x)), 0, pi/2)
pi/2*(b + 1 - sqrt(b^2+1))



### I( atan(k*cos(x)^2) / cos(x)^2 )
# Maths 505: A RIDICULOUSLY AWESOME INTEGRAL FEAT. φ & π
# https://www.youtube.com/watch?v=6CJ_Sal55uM

### I( atan(2*cos(x)^2) / cos(x)^2 )
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
pi/2 * Im((1+3i)^(1/2) - (1-3i)^(1/2));
pi/2i * sqrt(2 - 2*sqrt(3^2 + 1) + 0i);


### Gen: I( atan(k*cos(x)^2) / cos(x)^2 )
k = 5^(1/3)
integrate(\(x) atan(k*cos(x)^2) / cos(x)^2, 0, pi/2)
pi/2i * sqrt(2 - 2*sqrt(k^2 + 1) + 0i);


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

