########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### Polynomial Radicals
###
### draft v.0.1c


### Types:

# I( x^p / (x^n + 1)^k ) on [0, Inf]
# I( x^p / (x^n + 1)^k ) on [0, 1]


### History

# - moved Integrals with Radicals
#   to this file from file:
#   Integrals.Fractions.Unity.Definite.R;


######################

### Helper Functions

Euler = 0.57721566490153286060651209008240243079;

### I on [0, 1]
# I( x^p / (x^n + 1))
int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}
int.FrDU01 = function(n, p=0) {
	if( p!= 0) {
		# r = int.FrDU01(1, p) - int.FrDU01(n, p);
		r = digamma((p+1)) - digamma((p+1)/n) - log(n);
		return(r);
	}
	digamma(1/n) + Euler + log(n);
}

### I on [0, Inf]
int.FrUInf = function(n, p=0, pow=1, coeff=1) {
	k = 1/pow;
	tmp = sapply(p, function(p) {
		gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
	});
	tmp = sum(coeff * tmp);
	return(tmp);
}

##########################
##########################

################
### Radicals ###
################

### [0, Inf]
p = 1 - sqrt(2)
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n


###############

### [0, 1]

### I( 1 / (x^3 + 1)^(1/3) )
integrate(\(x) 1 / (x^3 + 1)^(1/3), 0, 1)
- (digamma(1/3) + Euler)/3 + 1/2*log((2^(2/3) + 2^(1/3) + 1)/3) +
	- 1/sqrt(3) * atan((2^(1/3) + 1/2) * 2/sqrt(3)) +
	+ 1/sqrt(3) * atan((1 + 1/2) * 2/sqrt(3));

### I( x / (x^3 + 1)^(1/3) )
# - simplified formula: see below;
integrate(\(x) x / (x^3 + 1)^(1/3), 0, 1)
1/2^(1/3) - gamma(2/3)^2/gamma(4/3) / 6

### I( 1 / (x^3 + 1)^(2/3) )
# - simplified formula: see below;
integrate(\(x) 1 / (x^3 + 1)^(2/3), 0, 1)
gamma(1/3)^2/gamma(2/3) / 6

### I( x / (x^3 + 1)^(2/3) )
integrate(\(x) x / (x^3 + 1)^(2/3), 0, 1)
- (digamma(2/3) + Euler)/3 +
	+ integrate(\(x) (x^2 - 1) / (x^3 - 1), 1, 2^(1/3))$value


### Pow = 5

### I( 1 / (x^5 + 1)^(1/5) )
integrate(\(x) 1 / (x^5 + 1)^(1/5), 0, 1)
- (digamma(1/5) + Euler)/5 +
	+ integrate(\(x) x^3 * (x-1) / (x^5 - 1), 1, 2^(1/5))$value;

### I( x^2 / (x^5 + 1)^(1/5) )
integrate(\(x) x^2 / (x^5 + 1)^(1/5), 0, 1)
1/2^(6/5) - gamma(3/5)^2/gamma(1/5) / 4

### I( 1 / (x^5 + 1)^(2/5) )
integrate(\(x) 1 / (x^5 + 1)^(2/5), 0, 1)
gamma(1/5)^2/gamma(2/5) / 10

### I( x / (x^5 + 1)^(2/5) )
integrate(\(x) x / (x^5 + 1)^(2/5), 0, 1)
- (digamma(2/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - x^2) / (x^5 - 1), 1, 2^(1/5))$value;

### I( x^2 / (x^5 + 1)^(3/5) )
integrate(\(x) x^2 / (x^5 + 1)^(3/5), 0, 1)
- (digamma(3/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - x) / (x^5 - 1), 1, 2^(1/5))$value;

### I( x^3 / (x^5 + 1)^(3/5) )
integrate(\(x) x^3 / (x^5 + 1)^(3/5), 0, 1)
1/2^(3/5) - gamma(4/5)^2/gamma(3/5) / 2

### I( x / (x^5 + 1)^(4/5) )
integrate(\(x) x / (x^5 + 1)^(4/5), 0, 1)
gamma(2/5)^2/gamma(4/5) / 10

### I( x^3 / (x^5 + 1)^(4/5) )
integrate(\(x) x^3 / (x^5 + 1)^(4/5), 0, 1)
- (digamma(4/5) + Euler)/5 +
	+ integrate(\(x) (x^4 - 1) / (x^5 - 1), 1, 2^(1/5))$value;


### Pow = 7

### I( 1 / (x^7 + 1)^(1/7) )
integrate(\(x) 1 / (x^7 + 1)^(1/7), 0, 1)
- (digamma(1/7) + Euler)/7 +
	+ integrate(\(x) x^5 * (x - 1) / (x^7 - 1), 1, 2^(1/7))$value


### Note:
# - the fraction decomposition of x^p / (x^n - 1)
#   (for n, p = integers) is described in file:
#   Integrals.Fractions.Unity.R;


### Derivation:
# Limit: (p+1) -> 0
n = 3; p = -1 + 1E-6; k = 3;
integrate(\(x) 1 / (x^3 + 1)^(1/3)/x - 1/x * exp(-x), 0, Inf)
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n - gamma(p+1)
- (digamma(1/3) - 2*Euler)/3

#
integrate(\(x) 1/x * exp(-x) - 1/(x*(x+1)), 0, Inf)

#
integrate(\(x) x / (x^2 + x + 1), 1, 2^(1/3))
1/2*log((2^(2/3) + 2^(1/3) + 1)/3) - 1/sqrt(3)*atan((2^(1/3) + 1/2)*2/sqrt(3)) +
	+ 1/sqrt(3)*atan((1 + 1/2)*2/sqrt(3))

# simplifications:
- 1/sqrt(3) * atan((2^(1/3) + 1/2) * 2/sqrt(3)) +
	+ 1/sqrt(3) * atan((1 + 1/2) * 2/sqrt(3));
- 1/sqrt(3) * atan((2^(1/3)-1)*2/sqrt(3) / (1 + 4/3*(2^(1/3) + 1/2)*(1 + 1/2)))
- 1/sqrt(3) * atan((2^(1/3) - 1) / (2^(1/3) + 1) / sqrt(3));


##################

### Special Cases:

### I( 1 / (x^(2*n) + 1)^(1/n) )

###
n = sqrt(5)
integrate(\(x) 1 / (x^(2*n) + 1)^(1/n), 0, 1)
gamma(1/(2*n))^2 / gamma(1/n) / (4*n)


### Series:

### I( x / (x^3 + 1)^(1/3) )
integrate(\(x) 1 / (x^(3/2) + 1)^(4/3), 0, 1)
gamma(2/3)^2/gamma(4/3) / 3
#
integrate(\(x) x / (x^3 + 1)^(4/3), 0, 1)
gamma(2/3)^2/gamma(4/3) / 6
#
integrate(\(x) x / (x^3 + 1)^(1/3), 0, 1)
1/2^(1/3) - gamma(2/3)^2/gamma(4/3) / 6

### I( 1 / (x^3 + 1)^(2/3) )
integrate(\(x) 1 / (x^3 + 1)^(2/3), 0, 1)
gamma(1/3)^2/gamma(2/3) / 6


### Trivial:
# see section further below;
integrate(\(x) 1 / (x^3 + 1)^(4/3), 0, 1)
1/2^(1/3)

###
integrate(\(x) 1 / (x^(3/2) + 1)^(5/3), 0, 1)
integrate(\(x) 2 * x / (x^3 + 1)^(5/3), 0, 1)
1/2^(2/3)
#
integrate(\(x) x / (x^3 + 1)^(5/3), 0, 1)
1/2^(5/3)


### Series: Pow = 5

###
integrate(\(x) x^2 / (x^5 + 1)^(6/5), 0, 1)
gamma(3/5)^2/gamma(6/5) / 10
#
integrate(\(x) x^2 / (x^5 + 1)^(1/5), 0, 1)
1/2^(6/5) - gamma(3/5)^2/gamma(6/5) / 20


###
integrate(\(x) x^3 / (x^5 + 1)^(8/5), 0, 1)
gamma(4/5)^2/gamma(8/5) / 10
#
integrate(\(x) x^3 / (x^5 + 1)^(3/5), 0, 1)
1/2^(3/5) - 3*gamma(4/5)^2/gamma(8/5) / 10


### Other Series:

### Type: (x^n + 1)^(2/n)
integrate(\(x) 1 / (x^6 + 1)^(1/3), 0, 1)
gamma(1/6)^2 / gamma(1/3) / 12

###
integrate(\(x) 1 / (x^(6/5) + 1)^(5/3), 0, 1)
integrate(\(x) 5*x^4 / (x^6 + 1)^(5/3), 0, 1)
gamma(5/6)^2 / gamma(5/3) / (6/5)/2


#######################
#######################

### Integer Powers:

# TODO:

###
integrate(function(x) 1/(x^3 + 1)^5, 0, 1, rel.tol=1E-8)
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3 + 2/3 * int.FrU01(3, 0))));
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3))) +
	+ 11/12 * 8/9 * 5/6 * 2/3 * int.FrU01(3, 0);
# TODO:
# (3!/(3*2^4) + 11 * 2!/(3^2*2^3) + 11*8 * 1!/(3^3*2^2) + 11*8*5 * 0!(3^4*2)) / gamma(5)
1/2^4 / 12 + 1/2^3 * 11/12 / 9 + 1/2^2 * 11/12*8/9 / 6 + 1/2 * 11/12*8/9*5/6 / 3 +
	+ int.FrU01(3, 0) * gamma(5 - 1/3)/gamma(1 - 1/3)/gamma(5);

###
# I[3,6] = 1/2^5 / (3*5) + 14/15 * I[3,5];
integrate(function(x) 1/(x^3 + 1)^6, 0, 1, rel.tol=1E-8)
1/2^5 / 15 + 1/2^4 * 14/15 / 12 + 1/2^3 * 14/15*11/12 / 9 +
	+ 1/2^2 * 14/15*11/12*8/9 / 6 + 1/2 * 14/15*11/12*8/9*5/6 / 3 +
	+ int.FrU01(3, 0) * gamma(6 - 1/3)/gamma(1 - 1/3)/gamma(6);


#############

### Radicals:

###
n = 5
I0 = integrate(\(x) sqrt(x^n + 1), 0, 1)
Ii = integrate(\(x) 1/sqrt(x^n + 1), 0, 1)
(n+2)/2*I0$value - n/2*Ii$value - sqrt(2) # == 0


###
n = 5
I1 = integrate(\(x) (x^n + 1)^(1/3), 0, 1)
I2 = integrate(\(x) (x^n + 1)^(2/3), 0, 1)
In1 = integrate(\(x) 1/(x^n + 1)^(1/3), 0, 1)
In2 = integrate(\(x) 1/(x^n + 1)^(2/3), 0, 1)
(n+3)/3*I1$value - n/3*In2$value - 2^(1/3) # == 0
(2*n+3)/3*I2$value - 2*n/3*In1$value - 2^(2/3) # == 0


### [old]
p = sqrt(5) - 2
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
IInf = gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n
print(IInf)

# on [0, 1]
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=1)$value
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=1, upper=Inf)$value
# TODO: find formula for -0.3248496;


#####################
#####################

### Diff Type

### I( x^p * (1 - x^n)^r )
### on [0, 1]

###
p = sqrt(3); r = sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^r, lower=0, upper=1)
gamma((p+1)/n) * gamma(r+1) / gamma((p+1)/n + r + 1) / n

### r > -1; p > - 1;
p = sqrt(7); r = - 1/sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p * (1 - x^n)^r, lower=0, upper=1)
gamma((p+1)/n) * gamma(r+1) / gamma((p+1)/n + r + 1) / n


### Special Cases:

### I( (1 - x^n)^r )
r = sqrt(7); n = sqrt(5)
integrate(function(x) (1 - x^n)^r, lower=0, upper=1)
gamma(1/n) * gamma(r+1) / gamma(1/n + r + 1) / n


### Fractions:
# works for r < 2; (e.g. up to r = 1.95)
p = sqrt(7); r = sqrt(3); n = 1/sqrt(5)
integrate(function(x) x^p / (1 - x)^r - n^r * x^p / (1 - x^n)^r, lower=0, upper=1)
gamma(p+1) * gamma(1 - r) / gamma(p - r + 2) +
	- n^r * gamma((p+1)/n) * gamma(1 - r) / gamma((p+1)/n - r + 1) / n;


ff = \(x) sapply(x, \(r) {
	integrate(function(x) x^p / (1 - x)^r - n^r * x^p / (1 - x^n)^r, 0, 1)$value;
})
fg = \(r) gamma(p+1) * gamma(1 - r) / gamma(p - r + 2) +
	- n^r * gamma((p+1)/n) * gamma(1 - r) / gamma((p+1)/n - r + 1) / n;
curve(ff(x), 1.5, 1.95, col="green")
curve(fg(x), add=TRUE, col="red", lty=2)


### Lim: r -> 1
p = sqrt(7); n = 1/sqrt(5)
integrate(function(x) x^p / (1 - x) - n * x^p / (1 - x^n), lower=0, upper=1)
digamma((p+1)/n) - digamma((p+1)) + log(n);


#####################

#############
### Other ###
#############

### Special Cases:

n = 8
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 9
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 10
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = sqrt(11)
integrate(function(x) 1/(1 - x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
integrate(function(x) 2/(2 - x^n)^(1+1/n), lower=0, upper=1)
# == 1


### Gen: I( 1 / (a - x^n)^(1/n + 1) )
a = sqrt(7); n = sqrt(5);
integrate(function(x) 1 / (a - x^n)^(1/n + 1), lower=0, upper=1)
1/(a-1)^(1/n) / a


### Gen: I( 1 / (x^n + a)^(1/n+1) )
n = sqrt(7); a = sqrt(5);
integrate(function(x) 1 / (x^n + a)^(1/n + 1), lower=0, upper=1)
1/(a+1)^(1/n) / a


### Gen: I( x^p / (x^n + a)^((p+1)/n + 1) )
p = sqrt(3); n = sqrt(7); a = sqrt(5);
integrate(function(x) x^p / (x^n + a)^((p+1)/n + 1), lower=0, upper=1)
1/(a+1)^((p+1)/n) / (a*(p+1))

