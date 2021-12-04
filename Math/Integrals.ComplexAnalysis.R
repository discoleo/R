####################
###
### Leonard Mada
###
### Integrals: Complex Analysis
### draft 0.1k


### Integrals:
# 1) I( (1 - cos(x^n)) / x^(n+1) )
# 2) I( x*sin(x)/(x^n + 1) )
# 3) I( x^k / (x^2+1) ), where -1 < k < 0;
#    Note: from [0, Inf];
# 3.b.) Generalization
#  - various (integer) powers of n:
#    I( x^k / (x^n+1) ), where -1 < k < 0;


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


### Test: Base-Case
k = - sqrt(2) / 2;
integrate(function(x) x^k /(x+1), lower=0, upper=Inf)
- pi / sin(k*pi)


###############
### Extensions:

### I( x^k / (x^2+1) )
# with -1 < k < 1;

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


### Ex 4:
k = sqrt(2)/2;
integrate(function(x) x^k /(x^2+1), lower=0, upper=Inf)
pi*cos((k-1)/2*pi) / sin(k*pi)


### Variant:
k = sqrt(2)/2;
integrate(function(x) cosh(k*log(x)) /(x*cosh(log(x))), lower=0, upper=1)
integrate(function(x) cosh(k*x) / cosh(x), lower=0, upper=1000) # BUG: upper=Inf;
pi*cos((k-1)/2*pi) / sin(k*pi)


###################

### Generalization:
### x^k /(x^n+1)

### Ext: n = 3
# with -1 < k < 2

### Ex 1:
k = - sqrt(3)/3;
integrate(function(x) x^k /(x^3+1), lower=0, upper=Inf)
pi/3 * sum(exp(1i*(k+1)*pi/3 * c(1,3,5))) / sin(k*pi) / exp(k*pi*1i)
pi/3 * sum(-1, exp(1i*pi/3*(c(1,5)*(k+1) - 3*k))) / sin(k*pi)
- pi/3 * sum(1, 2*cos(2*(k+1)*pi/3)) / sin(k*pi)

### Ex 2:
k = sqrt(2)
integrate(function(x) x^k /(x^3+1), lower=0, upper=Inf)
- pi/3 * sum(1, 2*cos(2*(k+1)*pi/3)) / sin(k*pi)


### Ext: n = 5
k = - sqrt(3)/3;
integrate(function(x) x^k /(x^5+1), lower=0, upper=Inf)
pi/5 * sum(exp(1i*(k+1)*pi/5 * c(1,3,5,7,9))) / sin(k*pi) / exp(k*pi*1i)
pi/5 * sum(-1, exp(1i*pi/5*(c(1,3,7,9)*(k+1) - 5*k))) / sin(k*pi)
# TODO: simplify;


#####################
#####################

### I[0, Inf]( log(x) / (x^n + 1) )
# I = - (pi/n)^2 * cos(pi/n) / sin(pi/n)^2

# Ref: Complex Analysis: Integral of ln(x)/(x^n+1) using Contour Integration
# https://www.youtube.com/watch?v=Sj8IJOBK33w

### Note:
# I( 1 / (x^n + 1) ), on (0, Inf)
# = pi/n / sin(pi/n);
# - for the exact formula on an arbitrary domain, see:
#   Integrals.Fractions.Unity.R; (n = only integers)


### Test: Base-Case
n = sqrt(2);
integrate(function(x) log(x) / (x^n + 1), lower=0, upper=Inf)
- (pi/n)^2 * cos(pi/n) / sin(pi/n)^2


### Power = 2
# Integration by parts =>
n = sqrt(2);
integrate(function(x) log(x) / (x^n + 1)^2, lower=0, upper=Inf)
- pi/n^2 / sin(pi/n) - pi^2*(n-1)/n^3 * cos(pi/n) / sin(pi/n)^2


### log(x) / (x^n - 1)
# I[0, 1] == I[1, Inf];
n = 2
integrate(function(x) log(x) / (x^n - 1), 0, 1)
pi^2 / 8

###
n = 3
integrate(function(x) log(x) / (x^n - 1), 0, 1)$value +
	integrate(function(x) log(x) / (x^n - 1), 1, Inf)$value
pi^2 / n^3 * 4

###
n = 4
integrate(function(x) log(x) / (x^n - 1), 0, 1)$value +
	integrate(function(x) log(x) / (x^n - 1), 1, Inf)$value
pi^2 / n^3 * 8

###
n = 6
integrate(function(x) log(x) / (x^n - 1), 0, 1)$value +
	integrate(function(x) log(x) / (x^n - 1), 1, Inf)$value
pi^2 / n^3 * 24


###
n = sqrt(2);
n = 3
m  = complex(re=cos(2*pi/n), im=sin(2*pi/n))
m1 = complex(re=cos(pi/n), im=sin(pi/n))
integrate(function(x) log(x^n + 1) / (x^n + 1), lower=0, upper=Inf)
# TODO: wrong
-2i*(pi/n) * exp(pi*1i/n) / (1-m)
m1/((m1+1)*(m1^2-1))


##########################
##########################

### I( x^p / sqrt(x^n + 1) )
### I( x^p / (x^n + 1)^(1/k) )


### Derivation:

n = 4; k = 3; p = 1/5;
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf) # y = x^n =>
integrate(function(x) 1/n * x^((p+1)/n - 1) / (x+1)^(1/k), lower=0, upper=Inf)
# y = x / (x+1) => y = 1 - 1/(x+1)=>
integrate(function(x) 1/n * x^((p+1)/n - 1) * (1-x)^(1/k - (p+1)/n - 1), lower=0, upper=1)
beta((p+1)/n, 1/k - (p+1)/n) / n
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n

# Note:
# Simple Contour: does NOT work!
m  = complex(re=cos(2*pi/n), im=sin(2*pi/n))
mn = complex(re=cos(pi/n), im=sin(pi/n))


###
n = 3
integrate(function(x) 1 /sqrt(x^n+1), lower=0, upper=Inf)
gamma(1/n)*gamma(1/2 - 1/n)/n / sqrt(pi)


###
n = 5
integrate(function(x) 1 /sqrt(x^n+1), lower=0, upper=Inf)
gamma(1/n)*gamma(1/2 - 1/n) / n / sqrt(pi)


###
n = 7
integrate(function(x) 1 /sqrt(x^n+1), lower=0, upper=Inf)
gamma(1/n)*gamma(1/2 - 1/n) / n / sqrt(pi)


###
n = 9
integrate(function(x) 1 /sqrt(x^n+1), lower=0, upper=Inf)
gamma(1/n)*gamma(1/2 - 1/n) / n / sqrt(pi)


### sqrt(N)
# N > 4;
n = sqrt(5)
integrate(function(x) 1 /sqrt(x^n+1), lower=0, upper=Inf)
gamma(1/n)*gamma(1/2 - 1/n) / n / sqrt(pi)


###############
###############

### Higher Order Radicals

###
n = 5
integrate(function(x) 1 / (x^n+1)^(1/3), lower=0, upper=Inf)
gamma(1/n)*gamma(1/3 - 1/n) / gamma(1/3) / n

###
n = 6
integrate(function(x) 1 / (x^n+1)^(1/3), lower=0, upper=Inf)
gamma(1/n)*gamma(1/3 - 1/n) / gamma(1/3) / n

###
n = 7
integrate(function(x) 1 / (x^n+1)^(1/3), lower=0, upper=Inf)
gamma(1/n)*gamma(1/3 - 1/n) / gamma(1/3) / n

###
n = sqrt(11)
integrate(function(x) 1 / (x^n+1)^(1/3), lower=0, upper=Inf)
gamma(1/n)*gamma(1/3 - 1/n) / gamma(1/3) / n


#####
n = 7
k = 6
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n


#####
n = sqrt(7)
k = sqrt(5)
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n


#####
n = sqrt(7)
k = 1/sqrt(5)
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n


### Selected Cases:

###
n = 3/2; k = 1/2; # 1/k = 2;
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n
pi / sin(pi/3) * 2/9

###
n = 3/2; k = 2/3; # 1/k = 3/2;
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n
gamma(2/3)*gamma(5/6) / gamma(1/2) * 4/3


###
n = 3/4;
p = n-1; k = 2/3;
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n;
2/n;


#############
### Variants:
n = sqrt(7)
k = sqrt(5)
integrate(function(x) 1 / (x^n+1)^(1/k), lower=0, upper=Inf)
integrate(function(x) (x^(n/k-2) + 1) / (x^n+1)^(1/k), lower=0, upper=1)
integrate(function(x) 2^(1-1/k)*cosh((n/(2*k) - 1)*log(x)) / x / cosh(n/2*log(x))^(1/k), lower=0, upper=1)
integrate(function(x) 2^(1-1/k)*cosh((n/(2*k) - 1)*x) / cosh(n*x/2)^(1/k), lower=0, upper=100) # BUG: upper = Inf
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) / n

### =>
integrate(function(x) cosh((n/(2*k) - 1)*x) / cosh(n*x/2)^(1/k), lower=0, upper=100) # BUG: upper = Inf
gamma(1/n)*gamma(1/k - 1/n) / gamma(1/k) * 2^(1/k - 1) / n


### full variant:
### cosh(p*x) / cosh(n*x)^(1/k)
n = sqrt(19)
k = sqrt(3); p = sqrt(2)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
integrate(function(x) (x^(n/k - p - 2) + x^p) / (x^n+1)^(1/k), lower=0, upper=1)
integrate(function(x) 2^(1-1/k)*cosh((n/(2*k) - p - 1)*log(x)) / x / cosh(n/2*log(x))^(1/k), lower=0, upper=1)
integrate(function(x) 2^(1-1/k)*cosh((n/(2*k) - p - 1)*x) / cosh(n*x/2)^(1/k), lower=0, upper=100) # BUG: upper = Inf
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n

### =>
integrate(function(x) cosh((n/(2*k) - p - 1)*x) / cosh(n*x/2)^(1/k), lower=0, upper=100) # BUG: upper = Inf
gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) * 2^(1/k - 1) / n

### =>
n = sqrt(19); k = sqrt(3); p = sqrt(2);
integrate(function(x) cosh(p*x) / cosh(n*x)^(1/k), lower=0, upper=100) # BUG: upper = Inf
gamma(1/(2*k) - p/(2*n)) * gamma(1/(2*k) + p/(2*n)) / gamma(1/k) * 2^(1/k - 2) / n

