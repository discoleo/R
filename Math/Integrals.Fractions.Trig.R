

Euler   = 0.57721566490153286060651209008240243079;

### I on [0, 1]
# I( x^p / (x^n + 1))
int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}

#############

### Based on:
# Integrals.Exercises.R &
# Integrals.Fractions.Unity.Derived.R
# - but those files need major refactoring;
# - see also:
#   Maths 505: A STUNNING integral from the JEE (main) exam
#   https://www.youtube.com/watch?v=mjM2GVJusU4


###
integrate(function(x) 1 / (cos(x)^4 + sin(x)^4), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1) / (x^4 + 1), lower=0, upper=1)
pi/(2*sqrt(2))


###
integrate(function(x) 1 / (cos(x)^6 + sin(x)^6), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^2 / (x^6 + 1), lower=0, upper=1)
pi/2


###
integrate(function(x) 1 / (cos(x)^8 + sin(x)^8), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^3 / (x^8 + 1), lower=0, upper=1)
sum(c(1,3,3,1) * sapply(c(6,4,2,0), int.FrU01, n=8))
pi*(1/tan(pi/16) + 1/tan(7*pi/16) + 3/tan(pi*3/16) + 3/tan(pi*5/16))/16
# can be expanded explicitly

# TODO: alternative?
integrate(function(x) (x^2 + 4) / (x^4 + 4*x^2 + 2), lower=0, upper=-Inf)
# integrate(function(x) 2 / (x^2 + 4 + 2*sqrt(2)), lower=0, upper=Inf)$value +
pi/sqrt(4+2*sqrt(2)) +
integrate(function(x) (4 - sqrt(2)) / (x^4 + 4*x^2 + 2), lower=0, upper=Inf)$value


###
n = 10; n2 = 2*n;
integrate(function(x) 1 / (cos(x)^n + sin(x)^n), lower=0, upper=pi/4)
integrate(function(x) (x^2 + 1)^(n/2 - 1) / (x^n + 1), lower=0, upper=1)
pi*(1/tan(pi/n2) + 1/tan(9*pi/n2) + 4/tan(pi*3/n2) + 4/tan(pi*7/n2) +
	+ 6/tan(pi*5/n2))/n2


#################
#################

###
integrate(\(x) sin(x) / (sin(x)^6 + cos(x)^6), 0, pi/4)
integrate(\(x) sin(x) / (1 - 3*sin(x)^2*cos(x)^2), 0, pi/4)
integrate(\(x) 1 / (1 - 3*x^2*(1-x^2)), sqrt(2)/2, 1)
integrate(\(x) 1 / (3*x^4 - 3*x^2 + 1), sqrt(2)/2, 1)
integrate(\(x) sqrt(3) / (x^4 - 3*x^2 + 3), sqrt(6)/2, sqrt(3))

# alternative:
integrate(\(x) sin(x) / (sin(x)^6 + cos(x)^6), 0, pi/4)
integrate(\(x) x * (x^2 + 1)^(3/2) / (x^6 + 1), 0, 1)
integrate(\(x) 1/2 * (x + 1)^(3/2) / (x^3 + 1), 0, 1)
integrate(\(x) x^4 / ((x^2-1)^3 + 1), 1, sqrt(2))

# =>
integrate(\(x) sqrt(3) / (x^4 - 3*x^2 + 3), sqrt(6)/2, sqrt(3));
1i*log( (1 - sqrt(r[1]/3)) * (sqrt(2)/2 + sqrt(r[1]/3))) / sqrt(r[1]) +
	- 1i*log( (1 - sqrt(r[2]/3)) * (sqrt(2)/2 + sqrt(r[2]/3))) / sqrt(r[2]) +
	- log(6) * sin(pi/12) / 3^(1/4) +
	+ 2*pi * cos(pi/12) / 3^(5/4);


# Derivation:
r = (sqrt(3) + c(-1i,1i)) * sqrt(3)/2;
rd = 1/2 - c(-1i,1i)*sqrt(3)/6;
rr = - c(-1i,1i)*sqrt(3)/6;
rfr = 1 + c(-1, 1) * sqrt(3)*1i;
integrate(\(x) sqrt(3) / (x^4 - 3*x^2 + 3), sqrt(6)/2, sqrt(3));
1i*log( (sqrt(3) - sqrt(r[1])) / (sqrt(3) + sqrt(r[1])) ) / (2*sqrt(r[1])) +
	- 1i*log( (sqrt(3) - sqrt(r[2])) / (sqrt(3) + sqrt(r[2])) ) / (2*sqrt(r[2])) +
	- 1i*log( (sqrt(6)/2 - sqrt(r[1])) / (sqrt(6)/2 + sqrt(r[1])) ) / (2*sqrt(r[1])) +
	+ 1i*log( (sqrt(6)/2 - sqrt(r[2])) / (sqrt(6)/2 + sqrt(r[2])) ) / (2*sqrt(r[2]))
1i*log( (1 - sqrt(r[1]/3))^2 / (1 - r[1]/3) ) / (2*sqrt(r[1])) +
	- 1i*log( (1 - sqrt(r[2]/3))^2 / (1 - r[2]/3) ) / (2*sqrt(r[2])) +
	- 1i*log( (sqrt(2)/2 - sqrt(r[1]/3))^2 / (1/2 - r[1]/3) ) / (2*sqrt(r[1])) +
	+ 1i*log( (sqrt(2)/2 - sqrt(r[2]/3))^2 / (1/2 - r[2]/3) ) / (2*sqrt(r[2]));
1i*log( (1 - sqrt(r[1]/3))^2 / rd[1] ) / (2*sqrt(r[1])) +
	- 1i*log( (1 - sqrt(r[2]/3))^2 / rd[2] ) / (2*sqrt(r[2])) +
	- 1i*log( (sqrt(2)/2 - sqrt(r[1]/3))^2 / rr[1] ) / (2*sqrt(r[1])) +
	+ 1i*log( (sqrt(2)/2 - sqrt(r[2]/3))^2 / rr[2] ) / (2*sqrt(r[2]));
1i*log( (1 - sqrt(r[1]/3))^2 / (sqrt(2)/2 - sqrt(r[1]/3))^2 *
		1/rfr[1] ) / (2*sqrt(r[1])) +
	- 1i*log( (1 - sqrt(r[2]/3))^2 / (sqrt(2)/2 - sqrt(r[2]/3))^2 *
		1/rfr[2] ) / (2*sqrt(r[2]));
1i*log( (1 - sqrt(r[1]/3)) / (sqrt(2)/2 - sqrt(r[1]/3))) / sqrt(r[1]) +
	- 1i*log( (1 - sqrt(r[2]/3)) / (sqrt(2)/2 - sqrt(r[2]/3))) / sqrt(r[2]) +
	- 1i * log(rfr[1]) * (1/(2*sqrt(r[1])) + 1/(2*sqrt(r[2]))) +
	+ 1i * log(4) / (2*sqrt(r[2]));
1i*log( (1 - sqrt(r[1]/3)) * (sqrt(2)/2 + sqrt(r[1]/3))) / sqrt(r[1]) +
	- 1i*log( (1 - sqrt(r[2]/3)) * (sqrt(2)/2 + sqrt(r[2]/3))) / sqrt(r[2]) +
	- 1i * log(6) * (sqrt(r[1]) - sqrt(r[2])) / (2*sqrt(3)) +
	+ pi * (sqrt(r[1]) + sqrt(r[2])) / (3*sqrt(3));


# Fraction Decomposition:
r = (3 + c(-1,1)*sqrt(3)*1i)/2;
r = (sqrt(3) + c(-1i,1i)) * sqrt(3)/2;
x = 5^(1/3)
1 / (x^4 - 3*x^2 + 3) # ==
(1/(x^2 - r[1]) - 1/(x^2 - r[2])) / sqrt(3) * 1i
(1i/(x - sqrt(r[1])) - 1i/(x + sqrt(r[1]))) / (2*sqrt(3*r[1])) +
	- (1i/(x - sqrt(r[2])) - 1i/(x + sqrt(r[2]))) / (2*sqrt(3*r[2]))

