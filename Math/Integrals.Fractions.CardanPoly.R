
### Leonard Mada
###
### Integrals: Polynomial Fractions
### Cardan-Type Polynomials
###
### draft v.0.2b


############

### History

# draft v.0.2b:
# - improved P5 variant;
# - new P7 variant (needs thorough testing);
# draft v.0.2a:
# - systematic approach to these polynomials;
# - fraction decomposition for P3;
# draft v.0.1:
# - this is the nicer decomposition, using conjugate roots;
# - basic decomposition of P5;
# - a previous solution used all the individual order 1 polynomials:
#  -- the old approach yields a compact solution as well;
#  -- but the current approach seems better;


##########

### Terminology

# F(k) = x^k / Q(x);
# I(k) = Integral x^k / Q(x) dx;
# where Q(x) = Cardan-type polynomial;

qdiv.f = function(x, c, d, k=0, n=3) {
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
I.f = function(lim, c, d, k=0, n=3) {
	integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=k, n=n)$value
}

#################


##########
### P3 ###
# 1/(x^3 - 3*c*x - 2*d)

### Fraction Decomposition
# F(0) = a / (x - r) - (a*x + b0) / (x^2 + r*x + r^2 - 3*c)
# where r = r[0] = p + q;

### ()*x:
# 2*r*a - b0 = 0 =>
# b0 = 2*r*a;
### Free Term:
# (r^2 - 3*c)*a + r*b0 = 1
# (r^2 - 3*c)*a + 2*r^2*a = 1
# 3*(r^2 - c)*a = 1

# a = 1/3 * 1/(r^2 - c) = 1/3 * (p - q)/(p^3 - q^3)
# b0 = 2/3 * r*(p - q)/(p^3 - q^3)


### Test: Fraction
x = 5 # arbitrary x;
c = 1
d = 2
#
det = sqrt(d^2 - c^3 + 0i)
p = (d + det)^(1/3); q = (d - det)^(1/3); r = p + q;
a = 1/3 * (p - q)/(p^3 - q^3)
b0 = 2*r*a;
#
1/(x^3 - 3*c*x - 2*d)
qdiv.f(x, c, d, k=0, n=3)
a / (x - r) - (a*x + b0) / (x^2 + r*x + r^2 - 3*c)

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

### Test
c = 1
d = 2
lim = c(3, 5)
# I(2)
integrate(qdiv.f, lower=lim[1], upper=lim[2], c=c, d=d, k=2, n=3)
1/3 * lnfr.f(lim, c=c, d=d, k=0, n=3) + c * I.f(lim, c=c, d=d, k=0, n=3)
# I(1)
# TODO



##################
##################

##########
### P5 ###
# 1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d)

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


### Roots of unity
n = 5 # b0 is currently limited to n = 5!
m = unity.sum.f(n)
m.sum = m$m.sum

# Parameters: free to change
c = 1
d = 3
# Roots
det = sqrt(d^2 - c^n)
p = (d + det)^(1/n)
q = (d - det)^(1/n)
r = p*m$m + q/m$m
r
# Coefficents of Partial Fractions
# b0 = 1 / (1/(r[1]*r[1]) + 2/(r[2]*r[5]) + 2/(r[3]*r[4]) ) / r[1] / (2*d)
b0 = 1/5 * (p-q)/(p^5-q^5) # TODO: check if always valid!
b = -2 * b0 * r[1] # ALL b are the same;
a = b0 * m.sum # TODO: check if correct;
c(a, b0, b) # displays only 1 coeff. b
### Test
x = 3 # any value - for testing the fraction;
#
1/(x^5 - 5*c*x^3 + 5*c^2*x - 2*d) # ==
b0/(x-p-q) + sum( (a*x + b) / ((x - r[2:3]) * (x - r[5:4])) )



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

### Roots of unity
n = 7
m = unity.sum.f(n)
m.sum = m$m.sum

# Parameters: free to change
c = 1
d = 3
# Roots
det = d^2 - c^n; det = if(det >= 0) sqrt(det) else sqrt(det + 0i);
p = (d + det)^(1/n)
q = (d - det)^(1/n)
r = p*m$m + q/m$m
r
# Coefficents of Partial Fractions
b0 = 1/n * (p-q)/(p^n-q^n) # TODO: check if always valid!
b = -2 * b0 * r[1] # TODO: check if ALL b are the same;
a = b0 * m.sum # TODO: check if correct
c(a, b0, b) # displays only 1 coeff. b
### Test
x = 3 # any value - for testing the fraction;
#
1/(x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d) # ==
b0/(x-p-q) + sum( (a*x + b) / ((x - r[2:4]) * (x - r[7:5])) )




##################
##################

##################
### Derivation ###

# some of the steps used in the derivation:


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
