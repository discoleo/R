
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
### draft v.0.1c



### History

# v 0.1c:
# - all:
#   Integral x^k / (x^n - 1)^p dx; [0 <= k;]
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


