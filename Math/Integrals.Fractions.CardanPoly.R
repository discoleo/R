
### Leonard Mada
###
### Integrals: Polynomial Fractions
### Cardan-Type Polynomials
###
### draft v.0.2d-rTr3


############

### History

# draft v.0.2d - v.0.2d-rTr3:
# - work on P5 polynomial terms;
# - initial work on the P7 polynomial (v.0.2d-bis);
# - more formulas using n-th derivatives & other Transforms:
#  -- d[n](Q(x)) / Q(x) (P7 & generalizable) (v.0.2d-der);
#  -- Tr(Q(x)) / Q(x) (P7 & generalizable) (v.0.2d-rTr, v.0.2d-rTr2);
#  -- general formulas for some of the root-transforms (v.0.2d-rTr3);
# draft v.0.2c:
# - added fraction decomposition for polynomials of even power;
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


####################

####################
### Introduction ###

### Polynomials
# - Cardan-type polynomials:
#   see https://github.com/discoleo/R/blob/master/Math/Polynomials.CardanGeneralisation.R;

### Examples:
# P3: Integral 1 / (x^3 - 3*c*x - 2*d) dx;
# P5: Integral 1 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d) dx;
# P7: Integral 1 / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d) dx;
### Even order:
# P6: Integral 1 / (x^6 - 6*c*x^4 + 9*c^2*x^2 - 2*c^3 - 2*d) dx;
# P8: Integral 1 / (x^8 - 8*c*x^6 + 20*c^2*x^4 - 16*c^3*x^2 + 2*c^4 - 2*d) dx;


### Terminology
# F(k) = x^k / Q(x);
# I(k) = Integral x^k / Q(x) dx;
# where Q(x) = Cardan-type polynomial;


#################
### Relations ###

### Derivatives

### sum( 1/(x - r) )
# = d(Q(x)) / Q(x);
# where d() = 1st derivative;

### sum( 1 / ((x - r[i1])*(x - r[i2])) )
# = 1/2! * d[2](Q(x)) / Q(x);
# where i1 < i2 and d[2] = 2nd derivative;

### sum( 1 / ((x - r[i1])*(x - r[i2])*(x - r[i3])) )
# = 1/3! * d[3](Q(x)) / Q(x);
# where i1 < i2 < i3;

# [...]


### Transforms

### sum( r[i] / (x - r[i]) )
# = (E1*x^(n-1) - 2*E2*x^(n-2) + 3*E3*x^(n-3) - 4*E4*x^(n-4) + ... + (-1)^n * (n-1)*E[n-1]*x + (-1)^(n+1) * n*E[n]) / Q(x)

### sum( (r[i1] + r[i2]) / ((x - r[i1])*(x - r[i2])) )
# = (1*(n-1)*E1*x^(n-2) - 2*(n-2)*E2*x^(n-3) + 3*(n-3)*E3*x^(n-4) + ... + (-1)^(n-1) * (n-2)*2*E[n-2]*x + (-1)^n * (n-1)*1*E[n-1]) / Q(x)
# = sum( (-1)^(i+1) * (n-i)*i*E[i]*x^(n-i-1) ) / Q(x);
# where i = 1:(n-1);


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
expand.idgrid = function(n, k) {
	# expand grid to compute various root combinations
	id.l = rep(list(1:n), k)
	id.gr = expand.grid(id.l)
	for(i in 1:(k-1)) {
	for(j in (i+1):k) {
		id.gr = id.gr[id.gr[,i] < id.gr[,j],]
	}
	}
	return(id.gr)
}
mult.pfr = function(r, k=2, type=1) {
	# type: 0 => 1/..., 1 => sum(r[i])/...,
	#       k => prod(r[i])/...;
	n = length(r); id.all = 1:n
	gr = if(k == 1) matrix(id.all, ncol=1) else expand.idgrid(n, k=k);
	p = rep(0, n)
	rows = 1:nrow(gr);
	for(id in rows) {
		is.id = id.all %in% gr[id,]
		r.inv = r[ ! is.id]
		# print(r[is.id])
		p.m = rev(Poly(r.inv))
		if(type == 0) {
		} else if(type == 1) {
			p.m = p.m * sum(r[is.id]) # Sum
		} else if(type == k) {
			p.m = p.m * prod(r[is.id]) # Product
		}
		p.m = c(p.m, rep(0, n - length(p.m)))
		p = p + p.m
	}
	round0.p(p)
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

### Coefficients
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

##########
### Pn ###
##########

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
	# TODO: rootn() when (d +/- det) < 0;
	p = (d + det)^(1/n)
	q = (d - det)^(1/n)
	r = p*m$m + q/m$m
	# Coefficients of Partial Fractions
	b0 = 1/n * (p-q)/(p^n - q^n) # TODO: check if always valid!
	b = -2 * b0 * r[1] # ALL b are the same;
	a = b0 * m.sum # TODO: check if always correct;
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
### other

# TODO:
1/(x^5 - 5*c^2*x^4 + 5*c*(2*d)^2*x^2 - (2*d)^4) # ==


5*(x^4 - 3*c*x^2 + c^2) / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d)


#################
### Integrals ###

### 1.) I(0)
# Base-Integral: using fraction decomposition;
# Integral 1 / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d) dx;
# TODO: trivial


### 2.) Polynomials

### 1/(x - r)
# = (x^2 - r*x*(m+m^4) + r^2 - (m^4+m+3)*c)*(x^2 - r*x*(m^2+m^3) + r^2 - (m^3+m^2+3)*c) / Q(x)
# = (x^4 + r*x^3 + (r^2 - 5*c)*x^2 + (r^3 - 5*c*r)*x + r^4 - 5*c*r^2 + 5*c^2) / Q(x)
# I(...) = log(x - r);

### (x^4 - c*x^2 + c^2) / (x^5 - 5*c*x^3 + 5*c^2*x - 2*d)
# I(...) = 1/5 * log(Q(x));

### Remaining terms of decomposition;
### 1 / ((x - r[2])*(x - r[5]))
# = (x - r)*(x^2 - r*x*(m^2+m^3) + r^2 - (m^3+m^2+3)*c) / Q(x)
# = (x^3 + r*x^2*(m^4+m) - (r^2*(m^4+m) + (m^3+m^2+3)*c)*x - r^3 + (m^3+m^2+3)*c*r) / Q(x)

### 1 / ((x - r[3])*(x - r[4]))
# = (x - r)*(x^2 - r*x*(m+m^4) + r^2 - (m^4+m+3)*c) / Q(x)
# = (x^3 + r*x^2*(m^3+m^2) - (r^2*(m^3+m^2) + (m^4+m+3)*c)*x - r^3 + (m^4+m+3)*c*r) / Q(x)

# Sum =>
# 1 / ((x - r[2])*(x - r[5])) + 1 / ((x - r[3])*(x - r[4]))
# = (2*x^3 - r*x^2 + (r^2 - 5*c)*x - 2*r^3 + 5*c*r) / Q(x)
# Diff =>
# 1 / ((x - r[2])*(x - r[5])) - 1 / ((x - r[3])*(x - r[4]))
# = (m^4 - m^3 - m^2 + m) * (r*x^2 - (r^2 - c)*x - c*r) / Q(x)

### TODO:
# - solve individual Terms;


### [old]

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


###################
###################

### Even Powers ###

##########
### P6 ###
##########

### x^6 - 6*c*x^4 + 9*c^2*x^2 - 2*c^3 - 2*d

# Parameters: free to change
n = 6
c = 1
d = 3
# Roots & Decomposition
fr = decompose.fr(c(c, d), n=n)
fr
n.half = n / 2
### Test
x = 3 # any value - for testing the fraction;
# Note: b0 = -b!
1/(x^6 - 6*c*x^4 + 9*c^2*x^2 - 2*c^3 - 2*d) # ==
-fr$b/(x^2 - fr$r[1]^2) + sum( (fr$a*x + fr$b) / ((x - fr$r[2:n.half]) * (x - fr$r[n:(n.half+2)])) )



##########
### P8 ###
##########

### x^8 - 8*c*x^6 + 20*c^2*x^4 - 16*c^3*x^2 + 2*c^4 - 2*d

# Parameters: free to change
n = 8
c = 1
d = 3
# Roots & Decomposition
fr = decompose.fr(c(c, d), n=n)
fr
n.half = n / 2
### Test
x = 3 # any value - for testing the fraction;
# Note: b0 = -b!
1/(x^8 - 8*c*x^6 + 20*c^2*x^4 - 16*c^3*x^2 + 2*c^4 - 2*d) # ==
-fr$b/(x^2 - fr$r[1]^2) + sum( (fr$a*x + fr$b) / ((x - fr$r[2:n.half]) * (x - fr$r[n:(n.half+2)])) )



##################
##################


###################
### Integrals   ###
### Polynomials ###

##########
### P7 ###
##########

### x^k/(x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)

### Equations

### sum( 1/(x - r) )
# = d(Q(x)) / Q(x)
# counterintuitive, but logical;
### sum( r/(x - r) )
# = (E1*x^6 - 2*E2*x^5 + 3*E3*x^4 - 4*E4*x^3 + 5*E5*x^2 - 6*E6*x + 7*E7) / Q(x)
# = (14*c*x^5 - 56*c^2*x^3 + 42*c^3*x + 14*d) / Q(x)
# = 14*(c*x^5 - 4*c^2*x^3 + 3*c^3*x + d) / Q(x)


### sum( 1/ ((x - r[i])*(x - r[j])) ), where i < j
# = (21*x^5 + 15*E1*x^4 + 10*E2*x^3 + 6*E3*x^2 + 3*E4*x + E5) / Q(x)
# = (21*x^5 - 70*c*x^3 + 42*c^2*x) / Q(x)
# = 1/2 * d2(Q(x)) / Q(x);
### sum ( x / ...)
# = 1/2 * x * d2(Q(x)) / Q(x);
### sum( (r1 + r2)/((x - r1)*(x - r2)))
# = (6*E1, -10*E2, 12*E3, -12*E4, 10*E5, -6*E6) / Q(x)
# = (70*c*x^4 - 14*12*c^2*x^2 + 42*c^3) / Q(x)
# = 7*(10*c*x^4 - 24*c^2*x^2 + 6*c^3) / Q(x)


### sum( 1/ ((x - r[i])*(x - r[j])*(x - r[k])) ), where i < j < k
# = 1/3 * (105*x^4 - 210*c*x^2 + 42*c^2) / Q(x)
# = (35*x^4 - 70*c*x^2 + 14*c^2) / Q(x)
# = 1/6 * d3(Q(x)) / Q(x);
### sum ( x / ...)
# = 1/6 * x * d3(Q(x)) / Q(x);
# = 

### Test
n = 7
d = 3
c = 1
x = 3 # Test value
#
r = decompose.fr(c(c, d), n=n)
# T1
# trivial d(Q(x)) / Q(x);
sum( 1 / (x - r$r) )
7*(x^6 - 5*c*x^4 + 6*c^2*x^2 - c^3) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)
#
sum( r$r / (x - r$r) )
14*(c*x^5 - 4*c^2*x^3 + 3*c^3*x + d) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)

# T2
id.gr = expand.idgrid(n, 2)
sum( 1 / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])) )
7*(3*x^5 - 10*c*x^3 + 6*c^2*x) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)
#
sum( x / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])) )
7*(3*x^6 - 10*c*x^4 + 6*c^2*x^2) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)
#
sum( (r$r[id.gr[,1]] + r$r[id.gr[,2]]) / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])) )
7*(10*c*x^4 - 24*c^2*x^2 + 6*c^3) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)

#
id.gr = expand.idgrid(n, 3)
# T3
sum( 1 / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])*(x - r$r[id.gr[,3]])) )
(35*x^4 - 70*c*x^2 + 14*c^2) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)
#
sum( x / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])*(x - r$r[id.gr[,3]])) )
(35*x^5 - 70*c*x^3 + 14*c^2*x) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)
#
sum( x^2 / ((x - r$r[id.gr[,1]])*(x - r$r[id.gr[,2]])*(x - r$r[id.gr[,3]])) )
(35*x^6 - 70*c*x^4 + 14*c^2*x^2) / (x^7 - 7*c*x^5 + 14*c^2*x^3 - 7*c^3*x - 2*d)


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

#######################

### Special derivations

### Elementary Polynomials

### P7: r[i] * prod(x - r[-i])

### P7: x^4-term
sum( r1*(E2 - r1*(E1 - r1)) )
sum( r1*E2 - r1^2*E1 + r1^3 )
E1*E2 - E1*(E1^2 - 2*E2) + E1^3 - 3*E1*E2 + 3*E3
3*E3

### P7: x^3-term
sum( r1*(E3 - r1*(E2 - r1*(E1 - r1))) )
sum( r1*E3 - r1^2*(E2 - r1*(E1 - r1)) )
sum( r1*E3 - r1^2*E2 + r1^3*E1 - r1^4 )
E1*E3 - (E1^2 - 2*E2)*E2 + (E1^3 - 3*E1*E2 + 3*E3)*E1 - (E1^4 - 4*E1^2*E2 + 4*E3*E1 + 2*E2^2 - 4*E4)
4*E4


#######
### P7: (r[i1] + r[i2]) * prod(x - r[-i])
(6*E1, -10*E2, 12*E3, -12*E4, 10*E5, -6*E6)

### P7: x^5-term
6*E1

### P7: x^4-term
- sum( (r[1] + r[2])*(E1 - (r[1] + r[2])) )
- sum( (r[1] + r[2])*E1 - (r[1] + r[2])^2 )
- (6*E1^2 - 6*r^2 - 2*E2)
-10*E2

### P7: x^3-term
12*E3

### P7: x^2-term
-12*E4

### P7: x^1-term
10*E5

### P7: x^0-term
-6*E6




