########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### Polynomial Fractions: Unity
### Definite Integrals
###
### draft v.0.2a



### Definite Integrals

# I( 1 / (x^n + 1) )
# I( x^p / (x^n + 1) )
# I( x^p / (x^n + 1)^r ) on [0, Inf]
# I( 1/(1 - x) - n/(1 - x^n) ) on [0, 1]
# I( x^p * (1 - x^n)^r ) on [0, 1]

### Intervals:
# [0, Inf], [0, 1];

# - contains various Shortcut Formulas;
# - moved specific code to this file from file:
#   Integrals.Fractions.Unity.R;
# - [refactor] moved initial derivations to file:
#   Integrals.Fractions.Unity.Definite.Deriv.R;


### Sections:
# A.) Helper Formulas:
# - Trigonometric identities;
# B.) Definite Integrals


######################

### Helper Functions

constEuler = 0.57721566490153286060651209008240243079;
Euler = constEuler;


### I on [0, 1]
# I( x^p / (x^n + 1))

# - code based on the Digamma function:
#   enables continuous values of n & p;
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

int.FrUInf = function(n, p=0, pow=1, coeff=1) {
	k = 1/pow;
	tmp = sapply(p, function(p) {
		gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
	});
	tmp = sum(coeff * tmp);
	return(tmp);
}

# [old code]
# - works only with n, p = Integers;
# - formulas are slightly different for n = Odd vs Even;
intUnityI01 = function(n) {
	if(n %% 2 == 0) {
		return(intUnityI01Even(n));
	}
	n2 = n %/% 2;
	id = seq(n2);
	sum(cos(2*pi*id/n)*log(cos(pi*id/n))) * 2 / n +
		# simplification:
		# + sum(id*sin(2*pi*id/n)) * 2*pi / n^2;
		+ pi/(2*n) / sin(pi/n);
}
intUnityI01Even = function(n) {
	if(n %% 2 == 1) return(intUnityI01(n));
	isNotM4 = n %% 4 != 0;
	id = seq(1, n, by=2);
	if(isNotM4) id = id[ - (n+2)/4]; # cs == 0;
	cs = cos(id*pi/n); # sn = sin(id*pi/n);
	# Base-Formula:
	# fr = sum(2*cs*x + 2)/(x^2 + 2*cs*x + 1)) / n;
	# r = sum(cs*log(x^2 + 2*cs*x + 1)) + sum(2*atan((x - cs)/sn));
	r = sum(cs*log(cs + 1)) +
		# + sum(2*sn*atan((1 - cs)/sn)) - sum(2*sn*atan((0 - cs)/sn));
		# simplification: sum(2*sn*atan((1 - cs)/sn));
		# simplification: - 2*pi*sum(id*sn)/n;
		+ pi/2 / sin(pi/n);
	return(r/n);
}
### x^p:
intUnityI01WX = function(m, p=1) {
	if(n %% 2 == 0) {
		return(intUnityI01EvenWX(n, p=p));
	}
	n2 = n %/% 2;
	id = seq(n2);
	cs = cos(2*pi*id*(p+1)/n); csH = cos(pi*id/n);
	sign = if(p %% 2 == 0) 1 else -1;
	pi/(2*n) / sin(pi*(p+1)/n) + sign * sum(cs*log(csH)) * 2 / n;
}
intUnityI01EvenWX = function(n, p = 1) {
	if(n %% 2 == 1) return(intUnityI01WX(n, p=p));
	id = seq(1, n, by=2);
	cs = cos(id*pi*(p+1)/n); csH = cos(pi*id/n);
	sign = if(p %% 2 == 0) 1 else -1;
	r = pi/2 / sin(pi*(p+1)/n) + sign * sum(cs*log(csH + 1));
	return(r/n);
}

### Fractional n & p
# - can be converted to integers;

# I( 1 / (x^(n/2) + 1) )
intUnityI01Half = function(n) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	cs2 = 2*cs^2 - 1; sn2 = 2*sn*cs;
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( cs2*log(cs + 1) ) +
		- 2/n * sum( sn2 * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( sn2 * (1/2 - 2*seq((n-1)/2)/n) );
	return(2*int);
}
# I( 1 / (x^(n/p) + 1) )
# - works for p < n; added proper sign;
intUnityI01POld = function(n, p=3) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	csp = cos(p*seq((n-1)/2)*2*pi/n);
	snp = sin(p*seq((n-1)/2)*2*pi/n);
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( csp*log(cs + 1) ) +
		- 2/n * sum( snp * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( snp * (1/2 - 2*seq((n-1)/2)/n) );
	if(p %% 2 == 1) int = - int;
	return(p*int);
}
# simplification:
intUnityI01P = function(n, p=3) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	csp = cos(p*seq((n-1)/2)*2*pi/n);
	snp = sin(p*seq((n-1)/2)*2*pi/n);
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( csp*log(cs + 1) ) +
		- 2/n * sum( snp * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( snp/2 );
	if(p %% 2 == 1) int = - int;
	int = int + pi/sin(pi*p/n)/n;
	return(p*int);
}

### Plot
# sum( cos(pi*seq(...)*p/n) * log(cos(pi*seq(...)/n)) )
# [solved] (digamma(((p+1)/n+1)/2) - digamma((p+1)/n/2) - pi/sin(pi*(p+1)/n) ) / 4
plot.logcos = function(n, dx = 1/8, len=129, points=TRUE,
		type = "l", col.px = "red", title=TRUE, ...) {
	p = seq(-1 + dx, n - 1 - dx, length.out=len);
	r = sapply(p, function(p) {
		integrate(function(x) x^p / (x^n + 1), 0, 1)$value;
	})
	r = r - pi / sin(pi*(p+1)/n) / (2*n);
	plot(p, r, type=type, ...);
	if(title) title(paste0("n = ", n));
	# Explicit Points
	if(points) {
		isOdd = (n %% 2 == 1)
		if(isOdd) {
			# even seq:
			id = 2 * seq(1, (n-1)/2);
		} else {
			id = seq(1, n, by=2);
			# csH = 1 + cos(pi*id/n);
		}
		logcsH = log(cos(pi*id/(2*n)));
		px = seq(0, n-2, by=1);
		py = sapply(px, function(p) {
			sign = if(p %% 2 == 0) 1 else -1;
			sign * sum(cos(pi*id*(p+1)/n) * logcsH);
		})
		points(px, 2*py / n, col=col.px);
		return(py);
	}
	invisible()
}


####################
####################

#################
### Section A ###
#################

### Useful Formulas
# - Trigonometric identities;

### n = ODD
n = 11
p = 3
#
id = seq((n-1)/2)
sum(id * sin(2*pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/4/sin(pi*p/n)

### n = ODD
n = 11
p = 4
#
id = seq((n-1)/2)
sum(id * sin(2*pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/4/sin(pi*p/n)

### n = 4*k + 2
k = 3
p = 1
#
n = 4*k + 2;
id = seq(1, n, by=2);
sum(id * sin(pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/2/sin(pi*p/n)

### n = 4*k
# - same formula;
k = 3
p = 1
#
n = 4*k;
id = seq(1, n, by=2);
sum(id * sin(pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/2/sin(pi*p/n)


##############

### sum(sn)
n = 10 # EVEN!
id = seq(1, n, by=2);
sn = sin(pi*id/n);
sum(sn) - 1/sin(pi/n) # = 0


######################
######################

###################
###  Section B  ###
###  Integrals  ###
###################

### Infinite Integrals
# - shortcut formulas are available;

### I( x^k/(x^n + 1) )

###
n = 5
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n

###
n = 7
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^5/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^4/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n
#
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(3*pi/n) / n


####################

### Radicals

### Generalized Formula
### [0, Inf]

# Derivation:
# - see file: Integrals.ComplexAnalysis.R;
#   [but is based on beta-function]

p = sqrt(2); n = 5*sqrt(11); k = 2 + sqrt(3);
integrate(function(x) x^p/(x^n + 1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n

###
p = 1 - sqrt(2); n = sqrt(11); k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n

### Composite:
integrate(\(x) (x^(1/3) + 2*x^(1/4) - 1)/sqrt(x^5 + 1), 0, Inf)
int.FrUInf(5, c(1/3, 1/4, 0), pow = 1/2, coeff = c(1,2,-1))


### Non-Standard Powers

###
n = sqrt(5)
#
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n


###
n = 3; p = sqrt(2);
#
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n

###
n = 5; p = sqrt(3);
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n


###
n = sqrt(11); p = sqrt(3);
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n


##################

### Composite

###
integrate(\(x) 1 / ((x^2+1) * (x^5+1)), 0, Inf)
integrate(\(x) 1/2 * ((x^4 - x^3 - x^2 + x + 1)/(x^5+1) - (x-1)/(x^2+1)), 0, Inf)
pi/2 * sum(c(-1,-1,1,1)/5/sin(pi * c(4,3,2,1)/ 5), 1/2/sin(pi/2))

###
p = sqrt(2)
# Monomials x^4 & x: cannot be ignored!
integrate(\(x) x^p / ((x^2+1) * (x^5+1)), 0, Inf)
pi/2 * sum(c(1,-1,-1,1,1)/5/sin(pi * (p + c(5,4,3,2,1))/ 5),
	c(-1,1)/2/sin(pi*(p+c(2,1))/2))


############################
############################

### Special Finite Integrals
# - on [0, 1];
# - shortcut formulas are available:
#   based on the digamma function;


###
n = sqrt(31)
p = sqrt(5)
integrate(function(x) x^p / (x^n + 1), 0, 1)
(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n)


###
n = 9
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 14
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 15
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 16
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 18
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)


### Derivation: Specific Cases
# Note:
# - derivation moved to file:
#   Integrals.Fractions.Unity.Definite.Deriv.R;


##################

### "Continuity" & Approximation
# - "quasi-analytic" continuity possible;
# - [solved] using digamma function;

###
n = seq(3, 22)
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")

### [old] Only Odd Powers
n = seq(3, 21, by=2)
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")

###
n = c(seq(3,20), seq(21, 101, by=3))
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")

###
plot.logcos(9)
plot.logcos(10)

n = 10
plot.logcos(n)
p = seq(0, n - 2)
(digamma(((p+1)/n+1)/2) - digamma((p+1)/n/2) - pi/sin(pi*(p+1)/n) ) / 4


### Special Fractional Powers

# [old]
# - full generalization using digamma function;

n = 7
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;

###
n = 8
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;


##########################
##########################

##########################
### Half-Beta Function ###
# on [1/2, 1]

p = 1/5; n = 11/3;
integrate(\(x) x^p / (x^n + 1), 0, 1)
integrate(\(x) 1/n * x^((p+1)/n - 1) / (x + 1), 0, 1)
integrate(\(x) 1/n * (1/x - 1)^((p+1)/n - 1) / x, 1/2, 1)
integrate(\(x) 1/n * x^(-(p+1)/n) * (1 - x)^((p+1)/n - 1), 1/2, 1)
int.FrU01(n=n, p=p)
(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);

###
p = 1/7
integrate(\(x) x^(-p) * (1 - x)^(p-1), 1/2, 1)
(digamma((p+1)/2) - digamma(p/2)) / 2;


###
p = 1/5; n = 11/3;
integrate(\(x) x^p / (x^n + 1)^2, 0, 1)
integrate(\(x) 1/n * x^(1 - (p+1)/n) * (1 - x)^((p+1)/n - 1), 1/2, 1)
#
integrate(\(x) x^(1 - p) * (1 - x)^(p - 1), 1/2, 1)
(digamma((p + 0)/2) - digamma((p - 1)/2)) * (p - 1) / 2 - 1/2
#
integrate(\(x) x^(- p) * (1 - x)^p, 1/2, 1)
(digamma((p + 1)/2) - digamma(p/2)) * p / 2 - 1/2

#
integrate(\(x) x^p / (x^n + 1)^2, 0, 1)
(digamma(((p+1)/n + 0)/2) - digamma(((p+1)/n - 1)/2)) * ((p+1)/n - 1) / (2*n) - 1/(2*n)

####################
####################

### 1 / (x^n - 1)

# I() on [0, 1]
n = 5
#
id = seq(floor((n - 1)/2))
cs = cos(2*pi*id/n); # ONLY 1 * cos()!
sn = sin(2*pi*id/n);
sn2 = sin(pi*id/n);
integrate(\(x) (x^3 + 2*x^2 + 3*x + 4) / (x^4 + x^3 + x^2 + x + 1), 0, 1)
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)

###
n = 7; # n = integer: odd & even!
#
id = seq(floor((n - 1)/2))
cs = cos(2*pi*id/n);
sn = sin(2*pi*id/n);
sn2 = sin(pi*id/n);
integrate(\(x) (x^5 + 2*x^4 + 3*x^3 + 4*x^2 + 5*x + 6) /
	(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1), 0, 1)
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)


# - cs * log(x^2 - 2*cs*x + 1) + 2*sn*atan((x - cs)/sn)
- sum(cs * log(2 - 2*cs)) + 2*sum(sn*atan((1 - cs)/sn) - sn*atan(- cs/sn))
log(2)/2 - sum(cs * log(1 - cs)) + 2*sum(sn*atan((1 - cs)/sn) - sn*atan(- cs/sn))
log(2) - 2*sum(cs * log(sn2)) + pi/2 * cos(pi/n) / sin(pi/n)


### Fraction Decomposition:
x = sqrt(11)
#
n = 5
1/(x - 1) - n/(x^n - 1)
(x^3 + 2*x^2 + 3*x + 4) / (x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1))
#
n = 6
1/(x - 1) - n/(x^n - 1)
(x^4 + 2*x^3 + 3*x^2 + 4*x + 5) / (x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1)) + 1/(x+1)
#
n = 7
1/(x - 1) - n/(x^n - 1)
(x^5 + 2*x^4 + 3*x^3 + 4*x^2 + 5*x + 6) / (x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1))
#
n = 8
1/(x - 1) - n/(x^n - 1)
(x^6 + 2*x^5 + 3*x^4 + 4*x^3 + 5*x^2 + 6*x + 7) /
	(x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1)
id = seq(floor((n - 1)/2))
cs = 2*cos(2*pi*id/n)
sum((2 - cs*x)/(x^2 - cs*x + 1)) + 1/(x+1)


####################
####################

### Radicals

### [0, Inf]
p = 1 - sqrt(2)
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n


###############

### [0, 1]
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


####################

# TODO:
# Integer Powers:

integrate(function(x) 1/(x^3 + 1)^5, 0, 1, rel.tol=1E-8)
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3 + 2/3 * int.FrU01(3, 0))));
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3))) +
	+ 11/12 * 8/9 * 5/6 * 2/3 * int.FrU01(3, 0);
# TODO:
# (3!/(3*2^4) + 11 * 2!/(3^2*2^3) + 11*8 * 1!/(3^3*2^2) + 11*8*5 * 0!(3^4*2)) / gamma(5)
1/2^4/12 + 11/12*(1/2^3/9 + 8/9*(1/2^2/6 + 5/6*(1/2/3))) +
	+ int.FrU01(3, 0) * gamma(5 - 1/3)/gamma(1 - 1/3)/gamma(5);


#######################
#######################

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


##########################
##########################

Euler = 0.57721566490153286060651209008240243079;

### Digamma Function
# see: Lines That Connect: Extending the Harmonic Numbers to the Reals
# https://www.youtube.com/watch?v=9p_U_o1pMKo


###
n = 7
integrate(\(x) (1 - (1-x)^n)/x, 0, 1)
integrate(\(x) (1 - x^n) / (1 - x), 0, 1)
sum(1/seq(n))
digamma(n + 1) + Euler


###
n = sqrt(5)
integrate(\(x) (1 - (1-x)^n)/x, 0, 1)
integrate(\(x) (1 - x^n) / (1 - x), 0, 1)
digamma(n + 1) + Euler


####################
####################

### Any n:
n = sqrt(7)
integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)
digamma(1/n) + Euler + log(n)


### I( x^p/(x-1) - n*x^p/(x^n - 1) )
p = sqrt(3); n = sqrt(11);
# can be: p > n;
integrate(\(x) x^p/(x-1) - n*x^p/(x^n - 1), 0, 1)
int.FrDU01(n, p)


# [alternative] [n = 5]
pi/tan(pi/n)/2 + 1/4*log(n) + sqrt(5)/2*atanh(1/sqrt(5))

### ODD n:
n = 5
id = seq(1, floor(n/2))
cs = cos(2*pi*id/n)
sn = sin(2*pi*id/n);sh = sin(pi*id/n);
integrate(\(x) 1/(1 - x) - n / (1 - x^n), 0, 1)
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) - sn * atan(-cs/sn))

sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) - sn * atan(-cs/sn))
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn) + sn * (pi/2 - 2*id*pi/n))
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn)) +
	+ 4*pi * sum(id*sn)/n - pi * sum(sn)
sum(cs * log(2 - 2*cs)) - 2 * sum(sn * atan((1 - cs)/sn)) +
	+ pi/sin(pi/n) - pi/tan(pi/n/2)/2
2*sum(cs * log(sh)) - log(2) +
	+ pi/sin(pi/n)/2 - pi/tan(pi/n/2)/2
digamma(1/n) + Euler + log(n) +
	+ pi/tan(pi/n)/2 - pi/tan(pi/n/2)/2 + pi/sin(pi/n)/2; # sum(TRIG) == 0;

# Helper:
pi * sum(sn)
pi/tan(pi/n/2)/2
#
sum(id*sn)
n/sin(pi/n)/4


#####################
#####################

### Radicals

### I( x^p * (1 - x^n)^r )
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

