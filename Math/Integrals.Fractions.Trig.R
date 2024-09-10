

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
# the pattern of: (...)^(n/2 - 1)
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
integrate(\(x) x^2 / (x^4 - 3*x^2 + 3), 1, sqrt(2))

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


########################
########################

### I( 1 / (sin(x)^5 + cos(x)^5) )
integrate(\(x) 1 / (sin(x)^5 + cos(x)^5), 0, pi/4)
# sin(x) + cos(x) => x
integrate(\(x) - 4/sqrt(2 - x^2) / (x^5 - 5*x), 1, sqrt(2))
# x => sqrt(2)*x
integrate(\(x) - 2*sqrt(2) / sqrt(1 - x^2) / (4*x^5 - 5*x), 1/sqrt(2), 1)
# x => sin(x):
integrate(\(x) - 2*sqrt(2) / (4*sin(x)^5 - 5*sin(x)), pi/4, pi/2)
# cos(x) => x
integrate(\(x) - 2*sqrt(2) / ((1 - x^2) * (4*(1-x^2)^2 - 5)), 0, 1/sqrt(2))
integrate(\(x) - 2*sqrt(2) / ((1 - x^2) * (4*x^4 - 8*x^2 - 1)), 0, 1/sqrt(2))
# Fraction Decomposition
integrate(\(x) 1/sqrt(10) / (x^2 - 1) *
	(1/(x^2 - 1 - sqrt(5)/2) - 1/(x^2 - 1 + sqrt(5)/2)), 0, 1/sqrt(2))
integrate(\(x) sqrt(2)/5 *
	(1/(x^2 - 1 - sqrt(5)/2) + 1/(x^2 - 1 + sqrt(5)/2) - 2/(x^2-1)), 0, 1/sqrt(2))
#
2*sqrt(sqrt(5) + 2)/5 * (pi/2 - atan(sqrt(sqrt(5) - 2))) +
	+ sqrt(sqrt(5) - 2)/5 *
		log((sqrt(sqrt(5) + 2) - 1)^2 / (sqrt(5) + 1)) +
	- 2*sqrt(2)/5 * log(sqrt(2)-1)

# TODO:
# - compact / generalized formula ???

### on [0, pi/2]
integrate(\(x) 1 / (sin(x)^5 + cos(x)^5), 0, pi/2)
integrate(\(x) 2 / (sin(x)^5 + cos(x)^5), 0, pi/4)
# see above;


###
lim = c(pi/7, pi/3)
integrate(\(x) 1 / (1 - sin(x)), lim[1], lim[2])
diff(tan(lim/2 + pi/4))

###
# explore alternative method: but the same;
integrate(\(x) 1 / (sin(x)^5 + cos(x)^5), 0, pi/4)
integrate(\(x) { y = sin(x) + cos(x); 1 / (y^5 - 5/2*(y^2-1)*y^3 + 5/4*(y^2-1)^2*y) }, 0, pi/4)
integrate(\(x) { y = sin(x + pi/4); - 2*sqrt(2) / (4*y^5 - 5*y) }, 0, pi/4)
integrate(\(x) { y = sin(x); - 2*sqrt(2) / (4*y^5 - 5*y) }, pi/4, pi/2)


### Pow = 7
integrate(\(x) 1 / (sin(x)^7 + cos(x)^7), 0, pi/4)
integrate(\(x) { y = sin(x) + cos(x); y2 = y^2 - 1;
	1 / (y^7 - 7/2*y2*y^5 + 7/2*y2^2*y^3 - 7/8*y2^3*y) }, 0, pi/4)
integrate(\(x) { y = sin(x + pi/4);
	4*sqrt(2) / (8*y^7 - 4*7*y^5 + 2*7*y^3 + 7*y) }, 0, pi/4)
integrate(\(x) { y = sin(x); 4*sqrt(2) / (8*y^7 - 4*7*y^5 + 2*7*y^3 + 7*y) }, pi/4, pi/2)
integrate(\(x) 8 / ((x^2 - 2) * (x^6 + x^4 - 9*x^2 - 1)), 0, 1)

# Fraction Decomposition
r = - 4*cos(1:3 * 2*pi/7) - 1;
x = r;
x^3 + x^2 - 9*x - 1 # == 0

integrate(\(x) 8 / ((x^2 - 2) * (x^6 + x^4 - 9*x^2 - 1)), 0, 1)
integrate(\(x) 8 / (x^2 - 2) *
	(1/(x^2 - r[1]) - 1/(x^2 - r[2])) / ((r[1] - r[2]) * (x^2 - r[3])), 0, 1)
integrate(\(x) 8 / (x^2 - 2) *
	( 1/(x^2 - r[1]) / ((r[1] - r[2])*(r[1] - r[3])) +
	+ 1/(x^2 - r[2]) / ((r[2] - r[1])*(r[2] - r[3])) +
	+ 1/(x^2 - r[3]) / ((r[3] - r[1])*(r[3] - r[2])) ), 0, 1)
integrate(\(x)
	( 8/(x^2 - r[1]) / ((r[1] - r[2])*(r[1] - r[3])*(r[1] - 2)) +
	+ 8/(x^2 - r[2]) / ((r[2] - r[1])*(r[2] - r[3])*(r[2] - 2)) +
	+ 8/(x^2 - r[3]) / ((r[3] - r[1])*(r[3] - r[2])*(r[3] - 2)) +
	- 8/(x^2 - 2)/7 ), 0, 1)

# prod(r[i] - r[-i]) == - (r^2 - 18*r - 3)/r
# sum(- r/(3*r^2 - 6*r - 1) ) == -1/7
# TODO

####################
####################

### on [0, pi/2]

### Base:
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) sin(x)^p * cos(x)^(n*k-p-2) / (sin(x)^n + cos(x)^n)^k, 0, pi/2)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
p = sqrt(2); q = sqrt(3); n = sqrt(5);
integrate(\(x) sin(x)^p * cos(x)^q / (sin(x)^n + cos(x)^n)^((p+q+2)/n), 0, pi/2)
gamma((p+1)/n) * gamma((q+1)/n) / gamma((p+q+2)/n) / n

###
p = sqrt(2); q = sqrt(3); n = sqrt(5); k = sqrt(3) - sqrt(2);
integrate(\(x) sin(x)^p * cos(x)^q / (sin(x)^n + cos(x)^n)^k, 0, pi/2)
integrate(\(x) x^p / ((x^n + 1)^k * (x^2 + 1)^((p+q-n*k+2)/2)), 0, Inf)

# TODO: ???

