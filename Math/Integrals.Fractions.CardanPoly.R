
### Leonard Mada
###
### Integrals: Polynomial Fractions
### Cardan-Type Polynomials
###
### draft v.0.2b-t3


############

### History

# draft v.0.2b - v.0.2b-t3:
# - improved P5 variant;
# - added P7 variant (needs thorough testing);
# - added P9 variant [v.0.2b-p9];
# - completed the tests for P3 integration [v.0.2b-t3];
# draft v.0.2a:
# - systematic approach to these polynomials;
# - fraction decomposition for P3;
# draft v.0.1:
# - this is the nicer decomposition, using conjugate roots;
# - basic decomposition of P5;
# - a previous solution used all the individual order 1 polynomials:
#  -- the old approach yields a compact solution as well;
#  -- but the current approach seems better;


###############

### Introduction

### Polynomials
# - Cardan-type polynomials:
#   see https://github.com/discoleo/R/blob/master/Math/Polynomials.CardanGeneralisation.R;

# Examples:
# P3: Integral 1 / (x^3 - 3*c*x - 2*d) dx;
# P5: Integral 1 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d) dx;
# P7: Integral 1 / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d) dx;


### Terminology
# F(k) = x^k / Q(x);
# I(k) = Integral x^k / Q(x) dx;
# where Q(x) = Cardan-type polynomial;


####################

### helper Functions
qdiv.f = function(x, c, d, k=0, n=3) {
	# the Fraction: F(k);
	if(k == 0) {
		xk = 1
	} else {
		xk = x^k
	}
	if(n == 3) {
		div = (x^3 - 3*c*x - 2*d)
	} else {
		div = 0
	}
	xk / div;
}
lnfr.f = function(lim, c, d, k=0, n=3) {
	# Note: log(Q(x)) = log(1/Fr_Q(x))
	# => lim[1] / lim[2];
	log(qdiv.f(lim[1], c=c, d=d, k=k, n=n) / qdiv.f(lim[2], c=c, d=d, k=k, n=n))
}
ln.fr = function(lim, r) {
	# log(x - r)
	log( (lim[2] - r) / (lim[1] - r))
}
I.f = function(lim, c, d, k=0, n=3) {
	# Integral of the fraction;
	# TODO: exact integral (based on fraction decomposition);
	integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=k, n=n)$value
}

#################


##########
### P3 ###
##########

### 1/(x^3 - 3*c*x - 2*d)

### Fraction Decomposition
# F(0) = a / (x - r) - (a*x + b) / (x^2 + r*x + r^2 - 3*c)
# where r = r[0] = p + q;

### ()*x:
# 2*r*a - b = 0 =>
# b = 2*r*a;
### Free Term:
# (r^2 - 3*c)*a + r*b = 1
# (r^2 - 3*c)*a + 2*r^2*a = 1
# 3*(r^2 - c)*a = 1

# a = 1/3 * 1/(r^2 - c) = 1/3 * (p - q)/(p^3 - q^3)
# b = 2/3 * r*(p - q)/(p^3 - q^3)

decompose.p3 = function(c, d) {
	det = sqrt(d^2 - c^3 + 0i)
	p = (d + det)^(1/3); q = (d - det)^(1/3);
	r = p + q;
	# Coeffs Fraction
	a = 1/3 * (p - q)/(p^3 - q^3)
	b = 2*r*a;
	return(list(r=r, a=a, b=b))
}

### Test: Fraction
x = 5 # arbitrary x;
c = 1
d = 2
#
q.dc = decompose.p3(c, d)
#
1/(x^3 - 3*c*x - 2*d)
qdiv.f(x, c, d, k=0, n=3)
q.dc$a / (x - q.dc$r) - (q.dc$a*x + q.dc$b) / (x^2 + x * q.dc$r + q.dc$r^2 - 3*c)


#############
### Integrals

### I(0)
# - can be computed trivially using the partial fraction decomposition;

### I(2)
# d(Q(x)) / Q(x) =
# 3*(x^2 - c) / Q(x)
# =>
# I(2) = ln(Q(x)) / 3 + c * I(0)

### I(1)
# F: 1/(x - r) = (x^2 + r*x + r^2 - 3*c) / Q(x);
# I(1) = (log(x - r) - I(2) - (r^2 - 3*c)*I(0)) / r

#########
### Test
c = 1
d = 2
lim = c(3, 5)
#
q.dc = decompose.p3(c, d);
r = q.dc$r;
### I(0)
integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=0, n=3)
I.f(lim, c,d, k=0, n=3)
### I(2)
integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=2, n=3)
1/3 * lnfr.f(lim, c=c, d=d, k=0, n=3) + c * I.f(lim, c=c, d=d, k=0, n=3)
### I(1)
integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=1, n=3)
(ln.fr(lim, r) - I.f(lim, c,d, k=2, n=3) - (r^2 - 3*c)*I.f(lim, c,d, k=0, n=3)) / r



##################
##################

unity.sum.f = function(n) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	n.half = (n-1)/2
	# m.m = matrix(m^(1:n.half), ncol=1)
	# m.m = cbind(m.m, 1/m.m)
	# m.sum = apply(m.m, 1, sum)
	m.all = m^(0:(n-1))
	m.half = m^(1:n.half)
	m.sum = m.half + 1/m.half
	return(list(m=m.all, m.sum=m.sum))
}
decompose.fr = function(coeff, n) {
	m = unity.sum.f(n)
	m.sum = m$m.sum
	# Roots
	c = coeff[1]; d = coeff[2];
	det = d^2 - c^n;
	if(det < 0) {
		det = sqrt(det + 0i);
	} else det = sqrt(det)
	#
	p = (d + det)^(1/n)
	q = (d - det)^(1/n)
	r = p*m$m + q/m$m
	# Coefficents of Partial Fractions
	b0 = 1/n * (p-q)/(p^n - q^n) # TODO: check if always valid!
	b = -2 * b0 * r[1] # ALL b are the same;
	a = b0 * m.sum # TODO: check if correct;
	return(list(r=r, b0=b0, a=a, b=b))
}

##########
### P5 ###
##########

### 1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d)

# Parameters: free to change
n = 5
c = 1
d = 3
# Roots & Decomposition
fr = decompose.fr(c(c, d), n=n)
fr
### Test
x = 3 # any value - for testing the fraction;
#
n.half = (n+1) / 2
1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d) # ==
fr$b0/(x - fr$r[1]) + sum( (fr$a*x + fr$b) / ((x - fr$r[2:n.half]) * (x - fr$r[n:(n.half+1)])) )



#########
### TODO:
1/(x^5 - 5*c^2*x^4 + 5*c*(2*d)^2*x^2 - (2*d)^4) # ==


5*(x^4 - 3*c*x^2 + c^2) / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d)


#################
### Integrals ###

### 1.) Base-Integral

# 1 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d)
# TODO: trivial



### 2.) Polynomials

# x^3 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d)
c = 1
d = 3
#
lower = 1
upper = 2
#
integrate(function(x) x^3 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d), lower=lower, upper=upper)
integrate(function(x) (2*d)^3 / (x^5 - 5*c^2*x^4 + 5*c*(2*d)^2*x^2 - (2*d)^4), lower=2*d/lower, upper=2*d/upper)

# TODO: exact integral;


##########
### P7 ###
##########

### 1/(x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)

# Parameters: free to change
n = 7
c = 1
d = 3
# Roots & Decomposition
fr = decompose.fr(c(c, d), n=n)
fr
### Test
x = 3 # any value - for testing the fraction;
#
n.half = (n+1) / 2
1/(x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d) # ==
fr$b0/(x - fr$r[1]) + sum( (fr$a*x + fr$b) / ((x - fr$r[2:n.half]) * (x - fr$r[n:(n.half+1)])) )


################
################

##########
### P9 ###
##########

### 1/(x^9 - 9*c*x^7 + 27*c^2*x^5 - 30*c^3*x^3 + 9*c^4*x - 2*d)

# Parameters: free to change
n = 9
c = 1
d = 3
# Roots & Decomposition
fr = decompose.fr(c(c, d), n=n)
fr
### Test
x = 3 # any value - for testing the fraction;
#
n.half = (n+1) / 2
1/(x^9 - 9*c*x^7 + 27*c^2*x^5 - 30*c^3*x^3 + 9*c^4*x - 2*d) # ==
fr$b0/(x - fr$r[1]) + sum( (fr$a*x + fr$b) / ((x - fr$r[2:n.half]) * (x - fr$r[n:(n.half+1)])) )





##################
##################

##################
### Derivation ###

### n = 5
# some of the steps used in the derivation:
# for n = 5

R = p + q
R11 = (m*p + m*q + m^4*p + m^4*q) # = R*(m + m^4)
R12 = (m^2*p + m^2*q + m^3*p + m^3*q) # = R*(m^2 + m^3)
R21 = (m*p^2 + m*q^2 + m^4*p^2 + m^4*q^2) # = (R^2 - 2*c)*(m + m^4)
R22 = (m^2*p^2 + m^2*q^2 + m^3*p^2 + m^3*q^2) # = (R^2 - 2*c)*(m^2 + m^3)

# =>

a1 + a2 + b0 # = 0
a1*R11 + a2*R12 + b0*R + b1 + b2 # = 0
- a1*(R21 + c - c*m^2 - c*m^3) - a2*(R22 + c - c*m - c*m^4) + b0*(R^2 - 5*c) + b1*R11 + b2*R12 # = 0
a1*(c*R12 - R^3 + 3*c*R) + a2*(c*R11 - R^3 + 3*c*R) + b0*(R^3 - 5*c*R) - b1*(R21 + c - c*m^2 - c*m^3) - b2*(R22 + c - c*m - c*m^4) # = 0
b0 * (R^4 - 5*c*R^2 + 5*c^2) + b1 * (c*R12 - R^3 + 3*c*R) + b2 * (c*R11 - R^3 + 3*c*R) # = 1

# =>
a1 + a2 + b0 # = 0
a1*R*(m + m^4) + a2*R*(m^2 + m^3) + b0*R + b1 + b2 # = 0
- a1*(R^2*(m + m^4) + 2*c - c*m - c*m^4) - a2*(R^2*(m^2 + m^3) + 2*c - c*m^2 - c*m^3) + b0*(R^2 - 5*c) + b1*R11 + b2*R12 # = 0
=> a1*(R^2*(m^2 + m^3) - 3*c + c*m + c*m^4) + a2*(R^2*(m + m^4) - 3*c + c*m^2 + c*m^3) + b1*R11 + b2*R12 # = 0
a1*(c*R12 - R^3 + 3*c*R) + a2*(c*R11 - R^3 + 3*c*R) + b0*(R^3 - 5*c*R) - b1*(R21 + c - c*m^2 - c*m^3) - b2*(R22 + c - c*m - c*m^4) # = 0
b0 * (R^4 - 5*c*R^2 + 5*c^2) + b1 * (c*R12 - R^3 + 3*c*R) + b2 * (c*R11 - R^3 + 3*c*R) # = 1


### Solution to Linear System:

c.m = matrix(
c(1,1,1,0,0,
R11, R12, R, 1, 1,
- (R21 + c - c*m^2 - c*m^3), - (R22 + c - c*m - c*m^4), (R^2 - 5*c), R11, R12,
(c*R12 - R^3 + 3*c*R), (c*R11 - R^3 + 3*c*R), (R^3 - 5*c*R), - (R21 + c - c*m^2 - c*m^3), - (R22 + c - c*m - c*m^4),
0, 0, (R^4 - 5*c*R^2 + 5*c^2), (c*R12 - R^3 + 3*c*R), (c*R11 - R^3 + 3*c*R)
), ncol=5, byrow=T)

c.m
c.sol = solve(c.m, c(rep(0,4), 1))
c.sol

###
a = c.sol[1:2]
b0 = c.sol[3]
b = c.sol[4:5]

# the parametric solution can be derived as well;
# [see section: Examples]
