
########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration of Polynomial Fractions
### - Roots of unity: Higher Powers
###   Integral( 1 / (x^n - 1)^p )dx
### - Roots of minus unity:
###   Integral( 1 / (x^n + 1)^p )dx
### - Polynomial fractions:
###   Integral( P(x) / (x^n - 1)^p )dx
###
### draft v.0.1f



### History

# v 0.1d:
# - added the base-cases for k < 0:
#   Integral 1 / (x^k * (x^n - 1)) dx;
# v 0.1c:
# - all:
#   Integral x^k / (x^n - 1)^p dx; [0 <= k, but works with k < 0 as well]
#   Integral x^(n+k) / (x^n - 1)^p dx;
# v.0.1b - v.0.1b-cor:
# - added:
#   Integral x / (x^n - 1)^p dx;
#   Integral x^2 / (x^n - 1)^p dx; [v.0.1b-bis]
#   Integral x^(n+1) / (x^n - 1)^p dx;
#   Integral x^(n+2) / (x^n - 1)^p dx; [v.0.1b-bis]
# - corrections in comments; [v.0.1b-cor]
# v.0.1a-fix:
# - bug fix: wrong sign in formulas;
# - solved:
#   Integral 1 / (x^n - 1)^p dx;
#   Integral x^n / (x^n - 1)^p dx;
# v.0.1a:
# - initial release;


###################

### Roots of unity: Baseline
#   Integral( P(x) / (x^n - 1) )dx
#   Integral( P(x) / (x^n + 1) )dx
# - see file: Integrals.Fractions.Unity.R


### Terminology

# F(k, n, p) = x^k / (x^n - 1)^p; # only the Fraction!
# I(k, n, p) = Integral( x^k / (x^n - 1)^p ) dx;


#######################


################
### Examples ###
################

### I(k, n, 1)
# - is computed in Integrals.Fractions.Unity.R;


##################
### I(n, n, p + 1)
### I(0, n, p + 1)

# (x^n - 1 + 1) / (x^n - 1)^(p+1)
# = 1/(x^n - 1)^p + 1/(x^n - 1)^(p+1)
# =>
# I(n, n, p+1) = I(0, n, p) + I(0, n, p+1)
# I(n, n, p+1) - I(0, n, p+1) = I(0, n, p)

# d/dx F(1, n, p)
# = (x^n - 1 - n*p*x^n) / (x^n - 1 )^(p+1)
# = -(n*p - 1)*F(n, n, p+1) - F(0, n, p+1)
# =>
# (n*p - 1)*I(n, n, p+1) + I(0, n, p+1) = - F(1, n, p);
# Note: Integral = Fraction!

### I(n, n, p+1)
# I(n, n, p+1) = 1/(n*p) * (I(0, n, p) - F(1, n, p))

### I(0, n, p+1)
# I(0, n, p+1) = -1/(n*p) * ((n*p - 1)*I(0, n, p) + F(1, n, p))

F.f = function(x, k, n, p) {
	if(k == 0) {
		r = 1 / (x^n - 1)^p
	} else {
		r = x^k / (x^n - 1)^p
	}
	r
}
F.range = function(lim, k, n, p) {
	F.f(lim[2], k,n,p) - F.f(lim[1], k,n,p)
}
I.f = function(lim, k, n, p) {
	integrate(F.f, lower=lim[1], upper=lim[2], k=k, n=n, p=p)
}

### Test
n = 5
p = 1
lim = c(2, 4)
### I(0, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=0, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 1, n, p) + (n*p - 1)*I.f(lim, 0, n, p)$value)


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(0, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=0, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 1, n, p) + (n*p - 1)*I.f(lim, 0, n, p)$value)


########################
########################


####################
### I(n+1, n, p + 1)
### I(1, n, p + 1)

# x*(x^n - 1 + 1) / (x^n - 1)^(p+1)
# = x/(x^n - 1)^p + x/(x^n - 1)^(p+1)
# =>
# I(n+1, n, p+1) = I(1, n, p) + I(1, n, p+1)
# I(n+1, n, p+1) - I(1, n, p+1) = I(1, n, p)

# d/dx F(2, n, p)
# = (2*x^(n+1) - 2*x - n*p*x^(n+1)) / (x^n - 1 )^(p+1)
# = -(n*p - 2)*F(n+1, n, p+1) - 2*F(1, n, p+1)
# =>
# (n*p - 2)*I(n+1, n, p+1) + 2*I(1, n, p+1) = - F(2, n, p);
# Note: Integral = Fraction!

### I(n+1, n, p+1)
# I(n+1, n, p+1) = 1/(n*p) * (2*I(1, n, p) - F(2, n, p))

### I(1, n, p+1)
# I(1, n, p+1) = -1/(n*p) * ((n*p - 2)*I(1, n, p) + F(2, n, p))


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(n+1, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k= n+1, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 2, n, p) - 2*I.f(lim, 1, n, p)$value)


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(0, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=1, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 2, n, p) + (n*p - 2)*I.f(lim, 1, n, p)$value), n, p + 1)
### I(0, n, p + 1)

# x*(x^n - 1 + 1) / (x^n - 1)^(p+1)
# = x/(x^n - 1)^p + x/(x^n - 1)^(p+1)
# =>
# I(n+1, n, p+1) = I(1, n, p) + I(1, n, p+1)
# I(n+1, n, p+1) - I(1, n, p+1) = I(1, n, p)

# d/dx F(2, n, p)
# = (2*x^(n+1) - 2*x - n*p*x^(n+1)) / (x^n - 1 )^(p+1)
# = -(n*p - 2)*F(n+1, n, p+1) - 2*F(1, n, p+1)
# =>
# (n*p - 2)*I(n+1, n, p+1) + 2*I(1, n, p+1) = - F(2, n, p);
# Note: Integral = Fraction!

### I(n+1, n, p+1)
# I(n+1, n, p+1) = 1/(n*p) * (2*I(1, n, p) - F(2, n, p))

### I(1, n, p+1)
# I(1, n, p+1) = -1/(n*p) * ((n*p - 2)*I(1, n, p) + F(2, n, p))


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(n+1, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k= n+1, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 2, n, p) - 2*I.f(lim, 1, n, p)$value)


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(1, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=1, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 2, n, p) + (n*p - 2)*I.f(lim, 1, n, p)$value)


########################


####################
### I(n+2, n, p + 1)
### I(2, n, p + 1)

# x^2*(x^n - 1 + 1) / (x^n - 1)^(p+1)
# = x^2/(x^n - 1)^p + x^2/(x^n - 1)^(p+1)
# =>
# I(n+2, n, p+1) = I(2, n, p) + I(2, n, p+1)
# I(n+2, n, p+1) - I(2, n, p+1) = I(2, n, p)

# d/dx F(3, n, p)
# = (3*x^(n+2) - 3*x^2 - n*p*x^(n+2)) / (x^n - 1 )^(p+1)
# = -(n*p - 3)*F(n+2, n, p+1) - 3*F(2, n, p+1)
# =>
# (n*p - 3)*I(n+2, n, p+1) + 3*I(2, n, p+1) = - F(3, n, p);
# Note: Integral = Fraction!

### I(n+2, n, p+1)
# I(n+2, n, p+1) = 1/(n*p) * (3*I(2, n, p) - F(3, n, p))

### I(2, n, p+1)
# I(2, n, p+1) = -1/(n*p) * ((n*p - 3)*I(2, n, p) + F(3, n, p))


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(n+2, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k= n+2, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 3, n, p) - 3*I.f(lim, 2, n, p)$value)


### Test
n = 5
p = 3
lim = c(1.1, 4)
### I(2, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=2, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, 3, n, p) + (n*p - 3)*I.f(lim, 2, n, p)$value)


########################


####################
### I(n+k, n, p + 1)
### I(k, n, p + 1)

# x^k*(x^n - 1 + 1) / (x^n - 1)^(p+1)
# = x^k/(x^n - 1)^p + x^k/(x^n - 1)^(p+1)
# =>
# I(n+k, n, p+1) = I(k, n, p) + I(k, n, p+1)
# I(n+k, n, p+1) - I(k, n, p+1) = I(k, n, p)

# d/dx F(k+1, n, p)
# = ((k+1)*x^(n+k) - (k+1)*x^k - n*p*x^(n+k)) / (x^n - 1 )^(p+1)
# = -(n*p - k - 1)*F(n+k, n, p+1) - (k+1)*F(k, n, p+1)
# =>
# (n*p - k - 1)*I(n+k, n, p+1) + (k+1)*I(k, n, p+1) = - F(k+1, n, p);
# Note: sum(Integrals) = Fraction!

### I(n+k, n, p+1)
# I(n+k, n, p+1) = 1/(n*p) * ((k+1)*I(k, n, p) - F(k+1, n, p))

### I(k, n, p+1)
# I(k, n, p+1) = -1/(n*p) * ((n*p - k - 1)*I(k, n, p) + F(k+1, n, p))

# Note:
# - the 2 formulas seem equivalent;
#   [are valid for any k]

### Test
n = 5
p = 3
k = 3
lim = c(1.1, 4)
### I(n+k, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k= n+k, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, k+1, n, p) - (k+1)*I.f(lim, k, n, p)$value)


### Test
n = 5
p = 3
k = 3
lim = c(1.1, 4)
### I(k, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=k, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, k+1, n, p) + (n*p - k - 1)*I.f(lim, k, n, p)$value)


### Test
n = 5
p = 3
k = sqrt(5) - sqrt(3) # works as well;
lim = c(1.1, 4)
### I(k, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=k, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, k+1, n, p) + (n*p - k - 1)*I.f(lim, k, n, p)$value)



#############################
#############################

### 1 / (x^k * (x^n - 1))
# Base-case: p = 1;

###  1 / (x * (x^n - 1))
# x^(n-1)/(x^n - 1) - 1/x = 1 / (x * (x^n - 1))
# I(-1, n, 1) = I(n-1, n, 1) - log(x)
# Note: I(n-1, n, 1) is also of type log(...);

###  1 / (x^n * (x^n - 1))
# 1/(x^n - 1) - 1/x^n = 1 / (x^n * (x^n - 1))
# I(-n, n, 1) = I(0, n, 1) + 1/(n-1) * 1/x^(n-1)

###  1 / (x^k * (x^n - 1))
# x^(n-k)/(x^n - 1) - 1/x^k = 1 / (x^k * (x^n - 1))
# I(-k, n, 1) = I(n-k, n, 1) + 1/(k-1) * 1/x^(k-1)


### Test: k = -1, p == 1
n = 5
lim = c(1.1, 4)
integrate(F.f, lower=lim[1], upper=lim[2], k=-1, n=n, p = 1)
(I.f(lim, n-1, n, 1)$value - log(lim[2]/lim[1]))


### Test: k = -n, p == 1
n = 5
lim = c(1.1, 4)
integrate(F.f, lower=lim[1], upper=lim[2], k=-n, n=n, p = 1)
(I.f(lim, 0, n, 1)$value + (1/lim[2]^(n-1) - 1/lim[1]^(n-1))/(n-1))


### Test: p == 1
n = 5
k = 3 # -k;
lim = c(1.1, 4)
integrate(F.f, lower=lim[1], upper=lim[2], k=-k, n=n, p = 1)
(I.f(lim, n-k, n, 1)$value + (1/lim[2]^(k-1) - 1/lim[1]^(k-1))/(k-1))


### Higher powers of p
# - the formulas in Section [A] can be applied with k < 0;

### Test
n = 5
p = 3
k = -sqrt(5) # works as well [even if the base-integral cannot be yet computed exactly];
lim = c(1.1, 4)
### I(k, n, 2)
integrate(F.f, lower=lim[1], upper=lim[2], k=k, n=n, p = p + 1)
-1/(n*p) * (F.range(lim, k+1, n, p) + (n*p - k - 1)*I.f(lim, k, n, p)$value)


######################
######################

### Mixed Fractions

### I( 1 / ((x^n + 1) * (x^(2*n) + 1)) )
n = sqrt(7)
integrate(\(x) 1 / ((x^n + 1) * (x^(2*n) + 1)), 0, 1)
(digamma((1/n + 1)/2) - digamma(1/n/2)) / (4*n) +
	- (digamma((1/n + 3)/4) - digamma((1/n + 1)/4)) / (8*n) +
	+ (digamma((1/n + 2)/4) - digamma((1/n + 0)/4)) / (8*n);


### I( x^p / ((x^n + 1) * (x^(2*n) + 1)) )
n = sqrt(7); p = sqrt(3)
integrate(\(x) x^p / ((x^n + 1) * (x^(2*n) + 1)), 0, 1)
(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (4*n) +
	- (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (8*n) +
	+ (digamma(((p+1)/n + 2)/4) - digamma(((p+1)/n + 0)/4)) / (8*n);


### Special Case: [0, Inf]
n = sqrt(7)
integrate(\(x) 1 / ((x^2 + 1) * (x^n + 1)), 0, Inf)
pi/4; # always


### I( 1 / sqrt(x^4 + 1) ) on [0, 1]
n = sqrt(5); # any n, including n= 0;
integrate(\(x) 1 / (sqrt(x*(x^2 + 1)) * (x^n + 1)), 0, Inf)
integrate(\(x) 1 / sqrt(x*(x^2 + 1)) * 1/2, 0, Inf)
integrate(\(x) 1 / sqrt(x*(x^2 + 1)), 0, 1)
integrate(\(x) 2 / sqrt(x^4 + 1), 0, 1)
gamma(1/4)^2 / gamma(1/2) / 4;


### Series: 3*n

### I( 1 / ((x^n + 1) * (x^(3*n) + 1)) )
n = sqrt(7)
integrate(\(x) 1 / ((x^n + 1) * (x^(3*n) + 1)), 0, 1)
1/(6*n) + (2*n-1)*(digamma((1/n + 1)/2) - digamma(1/n/2)) / (6*n^2) +
	- (digamma((1/n + 5)/6) - digamma((1/n + 2)/6)) / (18*n) +
	+ (digamma((1/n + 3)/6) - digamma((1/n + 0)/6)) / (18*n);


### Derivation:
integrate(\(x) 1 / ((x^n + 1) * (x^(3*n) + 1)), 0, 1)
integrate(\(x) 1/3/(x^n + 1)^2 - 1/3*(x^n + 1)/(x^(3*n) + 1) +
	+ 1/(x^(2*n) - x^n + 1) / (x^n + 1), 0, 1)
integrate(\(x) 1/3/(x^n + 1)^2 + 1/3/(x^n + 1) +
	- 1/3*(x^(2*n) - 1)/(x^(3*n) + 1), 0, 1)

###
integrate(\(x) 1/(x^n + 1)^2, 0, 1)
1/(2*n) + (n-1)*(digamma((1/n + 1)/2) - digamma(1/n/2)) / (2*n^2)

