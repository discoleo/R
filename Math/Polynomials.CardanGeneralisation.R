
### Generalisation of Cardan's Formula
### to Higher Order Polynomials
###
### (C) Leonard Mada
### draft v.0.2
###

### Citation:
### Leonard Mada. Generalisation of Cardan's Formula to Higher Order Polynomials.
### https://github.com/discoleo/R/blob/master/Math/Polynomials.CardanGeneralisation.R
### [I am open to collaboration, if anyone wants to publish properly]


### History:
# draft v.0.2:
# - renamed root() to solve.crd();
# - added section on Special cases: Bugs;
# Theory for Cardan-type Polynomials:
# - consistent theory developed during 2018-2019;
# - disseminated, but never oficially published;

####################

####################
### Introduction ###

# Cardan type polynomials:
# - can be easily solved (see below);
# - have very nice properties:
#  -- see Polynomial wavelets:
#   https://github.com/discoleo/R/blob/master/Math/Polynomial.Wavelet.R
#  -- polynomial fractions: very nice and compact integrals:
#   TODO: material on integrals of fractions;
# - these polynomials are a very special case of a broader class: Class 1 polynomials;

# Note:
# Terminology used is self developed and may differ from common mathematical terminology;
# [e.g. Class 1, Class 2, Class 3 polynomials]

### TODO:
# Prepare materials for:
# - exact integrals of the polynomial fractions;
# - Class 1 polynomials;
#   [the full generalisation]
# - other Classes (Class 2, Class 3);


### Theory

# let c, d = two parameters;
# p = (d + sqrt(d^2 - c^n))^(1/n)
# q = (d - sqrt(d^2 - c^n))^(1/n)

### Roots
### Base root:
# r = p + q
### All roots:
# r[j] = p*m^j + q*m^(-j)
# [all these n roots are presented in detail below]

### Properties of p & q:
# p^n + q^n = 2*d
# p*q = c

# r is the root of a polynomial of order n with *rational* coefficients based on (c, d);
# - if (c, d) are integers, the coefficients will also be integers, e.g.:

### Coefficients:
# [for n = odd]
# b[2*j-1] = C[j] * c^p,
# with p = (n+1)/2 - j, j = 1..(n-1)/2 and C[j] some constant;
# b[even] = 0 for n = odd;
# b0 = 2*d for n = odd;

### All Roots:
# r[j] = p*m^j + q*m^(-j), for j = 0..(n-1);
# where m^n = 1, m = root of unity of order n;

# Note:
# - this is how Cardan's formula should be tought!
# - it is fully generalisable;
# - all roots are properly covered & easily computable;

### Related/Interesting Materials:
# 1. Two Brilliant Ways to Solve This Non Factorable Cubic Equation
#    https://www.youtube.com/watch?v=Jz-Z0jfs2V4


##################

### Examples:

### "Quadratic" / Cardan-type Polynomials
# Odd-Order
f3  = function(x) x^3 - 3*c*x - 2*d
f5  = function(x) x^5 - 5*c*x^3 + 5*c^2*x - 2*d
f7  = function(x) x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d
f9  = function(x) x^9 - 9*c*x^7 + 27*c^2*x^5 - 30*c^3*x^3 + 9*c^4*x - 2*d
f11 = function(x) x^11 - 11*c*x^9 + 44*c^2*x^7 - 77*c^3*x^5 + 55*c^4*x^3 - 11*c^5*x - 2*d
f13 = function(x) x^13 - 13*c*x^11 + 65*c^2*x^9 - 156*c^3*x^7 + 182*c^4*x^5 - 91*c^5*x^3 + 13*c^6*x - 2*d
f15 = function(x) x^15 - 15*c*x^13 + 90*c^2*x^11 - 275*c^3*x^9 + 450*c^4*x^7 - 378*c^5*x^5 + 140*c^6*x^3 - 15*c^7*x - 2*d
# Even-Order
f2 = function(x) x^2 - 2*c - 2*d
f4 = function(x) x^4 - 4*c*x^2 + 2*c^2 - 2*d
f6 = function(x) x^6 - 6*c*x^4 + 9*c^2*x^2 - 2*c^3 - 2*d
f8 = function(x) x^8 - 8*c*x^6 + 20*c^2*x^4 - 16*c^3*x^2 + 2*c^4 - 2*d
f10 = function(x) x^10 - 10*c*x^8 + 35*c^2*x^6 - 50*c^3*x^4 + 25*c^4*x^2 -2*c^6 - 2*d
f12 = function(x) x^12 - 12*c*x^10 + 54*c^2*x^8 - 112*c^3*x^6 + 105*c^4*x^4 - 36*c^5*x^2 + 2*c^6 - 2*d
f14 = function(x) x^14 - 14*c*x^12 + 77*c^2*x^10 - 210*c^3*x^8 + 294*c^4*x^6 - 196*c^5*x^4 + 49*c^6*x^2 -2*c^7 - 2*d

### Roots:

# Note:
# - these formulas may fail with certain complex type coefficients
#   (see section on Special Cases);
# - there is a lot of boiler-plate code to handle:
#   ( - real_value)^(1/n);
solve.crd = function(c, d, n=3) {
	det = d^2 - c^n
	det = if(Im(det) == 0 && Re(det) >= 0) sqrt(Re(det))
		else if(Im(det) == 0) complex(re=0, im=sqrt(-Re(det)))
		else sqrt(det)
	if(Im(det) != 0 || Im(d) != 0) {
		pq = (c(d + det, d - det))^(1/n)
	} else {
		det = Re(det)
		pq = c(d + det, d - det)
		if(pq[1] >=0 && pq[2] >= 0) {
			pq = pq^(1/n)
		} else if(n %% 2 == 1) {
			pq.sign = sign(pq)
			pq = pq.sign * (pq * pq.sign)^(1/n)
		} else {
			# TODO:
			pq = pq^(1/n)
		}
	}
	#
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	m = m^(0:(n-1))
	r = pq[1] * m + pq[2] / m
	return(r)
}

round0 = function(x, tol=1E-12) {
	isZero = abs(Im(x)) < tol
	x[isZero] = Re(x[isZero])
	#
	isZero = abs(Re(x)) < tol
	isComplex = (Im(x)) != 0
	x[isZero & isComplex] = complex(re=0, im=Im(x[isZero & isComplex]))
	x[isZero & ! isComplex] = 0
	#
	return(x)
}

#############

### A.) Tests

n = 5
# Parameters:
c = 1
d = 2
### Roots
x = solve.crd(c, d, n)
x
# residual error
err = f5(x)
round0(err) # tol < 1E-12


### B.) Special cases:

### === Polynomial of Order 3 ===

# B.1.) Ordinary Case:
n = 3
#
c = 3
d = 2
### Roots
x = solve.crd(c, d, n)
x
# residual error
err = f3(x)
round0(err) # tol < 1E-12

# Note:
# - Cardan's solution covers the entire polynomial space for polynomials of order 3;
# - it only covers a special subclass for polynomials of higher order;
# - TODO:
#   presentation of Class 1 polynomials (which is a generalisation);

# B.2.) Inverted Root:
# x^3 - 3*c*x^2 - 2*d
n = 3 # fixed
c = 1
d = 2
# x^3 - 3*c*x^2 - 2*d
# x => 2*d/x
# 4*d^2 - 6*c*d*x - x^3
# x^3 + 6*c*d*x - 4*d^2
# d is now 2*d^2; c = -2*c*d;
x = solve.crd( -2*c*d, 2*d^2, n)
x = 2*d/x
x
# err
err = x^3 - 3*c*x^2 - 2*d
round0(err)


# B.3.) Shifted Root
# x^3 + b2*x^2 + b1*x + b0
# Step 1: shift by -b2/3;
# Step 2: solve using presented formula;
# Step 3: shift roots back;
# TODO: function/automation;
# [is already available in some of the other materials]


###########################

### C.) Special cases: Bugs

# the simple formula fails when c = m^j,
# where m^n = 1, m = root of unity of order n;

# Note:
# p*q = (c^n)^(1/n) = (m^(j*n))^(1/n)
# = 1^(1/n) = 1,
# BUT: 1 != m^j;

### Example: n = 3
# x^3 - 3*m*x - 2 = 0, where m^3 = 1;
# c = m, d = -1;
# => p*q = 1, but 1 != m;

m3 = complex(re=cos(2*pi/3), im=sin(2*pi/3))

# Solution:
# Step 1:
# x = y*m^2 =>
# y^3 - 3*y - 2 = 0
# => solve for y;
# Step 2:
# x = y * m^2;
y = solve.crd(1, 1, n=3)
x = y * m3^2
### Test
err = x^3 - 3*m3*x - 2
round0(err)
x



####################
####################

##################
### Properties ###

# these polynomials have nice properties;

### P.1.) Polynomial Wavelets
# see https://github.com/discoleo/R/blob/master/Math/Polynomial.Wavelet.R

### P.2.) Integrals of Polynomial Fractions
# the polynomial fractions have very nice compact exact integrals;

### TODO:
# - implement & describe the modern approach;
#  [the initial code was based on decomposition using fractions of order 1]


### n = 3:
### Integral: 1 / (x^3 - 3*c*x - 2*d) dx
c = 1
d = 3
### Exact integral:
# Step 1: compute p & q
x = solve.crd(c, d, 3)
x
x = x[1] # only real root
# Step 2:
# Exact integral:
exact = function(x, lower, upper) {
	A = -1/3 / (x^2 - c) # / (p^2 + q^2 + c);
	B = 2*A*x;
	C = B - A*x/2;
	D = (x^2 - 4*c) * 3/4
	D.sqr = sqrt(D)
	exact.f = function(y) A/2 * log((y^3 - 3*c*y -2*d) / (y-x)^3) + C/D.sqr * atan((y+x/2)/D.sqr)
	exact.f(upper) - exact.f(lower)
}
#
lower = 3
upper = 5
exact(x, lower, upper)
# Test
integrate(function(x) 1 / (x^3 - 3*c*x - 2*d), lower=lower, upper=upper)
# 0.05688931 with both methods;


### n = 5
### Integral: 1 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d) dx

# TODO:
# reimplement using modern approach;
# [initial approach used a full fraction decomposition; unpublished;]


#################
#################

#################
### Generator ###

# Even-powered Polynomials
poly.cardan.gen = function(n) {
	if(n %% 2 == 0) {
		coeff = list("c0"=1, "c2"=c(1, -2))
	}
	for(i in seq(4, n, by=2)) {
		coeff = expand(i, coeff)
	}
	return(coeff)
}
expand = function(n, coeff) {
	last = n/2
	ch = choose(n, seq(n/2, n-1, by=1))
	tmp = 0
	for(id in last:1) {
		tmp = tmp + c(rep(0,last-id), coeff[[id]]) * ch[id]
	}
	ch = c(1, -tmp)
	coeff[[paste("c", n, sep="")]] = ch
	return(coeff)
}

###
coeff = poly.cardan.gen(16)
coeff



########################
########################

### some variants

x^5 - 5*c^2*x^4 + 5*c*(2*d)^2*x^2 - (2*d)^4
x^5 - 5*c^2*x^4 + 20*c*d^2*x^2 - 16*d^4

# x^5 - 5*B^2/A*x^4 + 20*B*x^2 - 16*A
A = 3
B = 1
#
d = A^(1/4)
c = B/d^2
#
x = solve.crd(c, d, n=5)
x = 2*d/x
err = x^5 - 5*B^2/A*x^4 + 20*B*x^2 - 16*A
round0(err)
x


### complicated Polynomial with same solution:
(5/A * (B*x^2 - 2*A)^2 - 4*A)^(1/5) # ==
(2*A/B - 1/B * sqrt((x^5 + 4*A)*A/5))^(1/2)

err = (5/A * (B*x^2 - 2*A)^2 - 4*A)^2 -
(2*A/B - 1/B * sqrt((x^5 + 4*A)*A/5))^5
round0(err)



