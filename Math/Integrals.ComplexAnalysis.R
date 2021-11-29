####################
###
### Leonard Mada
###
### Integrals: Complex Analysis
### draft 0.1d


### Integrals:
# 1) I( (1 - cos(x^n)) / x^(n+1) )
# 2) I( x*sin(x)/(x^n + 1) )
# 3) I( x^k / (x^2+1) ), where -1 < k < 0;
#    Note: from [0, Inf];

########################
########################

### I( (1 - cos(x^n)) / x^(n+1) )

# Ref:
# Complex Analysis: Integral of (1-cos(x))/x^2 using Contour Integration
# https://www.youtube.com/watch?v=PmQiZ5ipAAA

### Extension: x^n
# Method 1: Complex Analysis;
# Method 2: u-Subst on Base-Integral;

integrate.cosFr = function(n, iter=2000, tol=1E-8, print=TRUE) {
	r = integrate(function(x) (1-cos(x^n))/x^(n+1), lower=-Inf, upper=Inf, subdivisions=iter, abs.tol=tol);
	if(print) print(r$abs.error);
	return(r$value);
}
integrate.cosFrAbs = function(n, iter=2000, tol=1E-8, print=TRUE) {
	r = integrate(function(x) (1-cos(abs(x)^n))/abs(x)^(n+1), lower=-Inf, upper=Inf, subdivisions=iter, abs.tol=tol);
	if(print) print(r$abs.error);
	return(r$value);
}
###
# r = pi / n;  # for n = odd;

###
n = 3
integrate.cosFr(n) * n

###
n = 4
integrate.cosFrAbs(n) * n

###
n = 5
integrate.cosFr(n) * n

###
n = 1/2
integrate.cosFrAbs(n) * n

#####################
#####################

### I [-Inf, + Inf] ( x*cos(x) / (x^(2*n) + 1) )

unityMinus = function(n) {
	id = seq(1, by=2, length.out=n);
	pin = pi / n;
	complex(re=cos(pin*id), im=sin(pin*id))
}
### Exact Formula:
divUnityMinus = function(n) {
	m = unityMinus(n)
	isPlus = (Im(m) > 0); # only Upper Half;
	mp = m[isPlus]; mn = m[ ! isPlus];
	m.div = sapply(seq_along(mp), function(id) {
		mr = c(mn, mp[-id]);
		prod(mp[id] - mr);
	})
	
	r = 2*pi*sum(mp*exp(1i*mp) / m.div);
	return(r);
}

################

# - only even n!

### n = 4
n = 4;

integrate(function(x) x*sin(x)/(x^n+1), lower=-Inf, upper=Inf)
divUnityMinus(n)


### n = 6
n = 6;

integrate(function(x) x*sin(x)/(x^n+1), lower=-Inf, upper=Inf)
divUnityMinus(n)


### n = 8
n = 8;

integrate(function(x) x*sin(x)/(x^n+1), lower=-Inf, upper=Inf)
divUnityMinus(n)


### n = 10
n = 10;

integrate(function(x) x*sin(x)/(x^n+1), lower=-Inf, upper=Inf)
divUnityMinus(n)


#####################
#####################

### I[0, Inf]( x^k / (x+1) )
# where -1 < k < 0;
# I = - pi / sin(k*pi)

# Ref: Complex Analysis: Integral of (x^n)/(x+1) using Contour Integration
# https://www.youtube.com/watch?v=zgLNBdtQT5Q


###
k = - sqrt(2) / 2;
integrate(function(x) x^k /(x+1), lower=0, upper=Inf)
- pi / sin(k*pi)


### Extension:
# I[0, Inf]( x^k / (x^2+1) )

### Ex 1:
k = - sqrt(2) / 2;
integrate(function(x) x^k /(x^2+1), lower=0, upper=Inf)
pi*cos((k-1)/2*pi) / sin(k*pi)


### Ex 2:
k = sqrt(2) - 2;
integrate(function(x) x^k /(x^2+1), lower=0, upper=Inf)
pi*cos((k-1)/2*pi) / sin(k*pi)


### Ex 3:
k = sqrt(3) - 2;
integrate(function(x) x^k /(x^2+1), lower=0, upper=Inf)
pi*cos((k-1)/2*pi) / sin(k*pi)

