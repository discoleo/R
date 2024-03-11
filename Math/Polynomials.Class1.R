
### Major Theories in Polynomials

### Polynomials of Class 1
### Roots with "Simple" Radicals
###
### Leonard Mada
###
### draft v.0.1e-expl2

### based on work during:
### 2018 - 2020


###############
### History ###
###############

### v.0.1e: [2021-12-10]
# - re-started to explore techniques to solve the quintic of Class 1;
### v.0.1d:
# - moved mpfr-functions to new file:
#   Polynomials.Helper.mpfr.R;
### v.0.1c:
# - added example with triple nested radicals;
#   [but overflows very easily!]
### v.0.1b:
# - added examples for entanglements
#   with multiple radicals;
#   [both nested & non-nested]
### v.0.1a:
# - initial draft posted on Github;
# - based largely on work during 2018-2020;


##############
### Theory ###
##############

### TODO:
# - full theory!

### Simple Radicals

### Base root:
# r = sum(s[j]*k^j),
# where k = K^(1/n) and s[j] and K are parameters;
# explicit summation:
# r = s[n-1]*k^(n-1) + s[n-2]*k^(n-2) + ... + s2*k^2 + s1*k + s0

### All roots:
# *ALL* roots can be derived from the Base root!
# r[j] = sum(s[id]*k[j]^id)
# where k[j] = k * m^j, m^n = 1;
# and index id goes from 0 to n-1;


####################

### Helper Functions

source("Polynomials.Helper.R")


###############


# Coefficients:
coef.p5 = function(b, K) {
	s1 = b[1]; s2 = b[2]; s3 = b[3]; s4 = b[4];
	c(- K*s1^5 + K^2*(5*(s1*s2^3*s3 - s1^2*s2*s3^2 - s1^2*s2^2*s4 + s1^3*s3*s4) - s2^5) +
			+ K^3*(5*(s1*s2*s4^3 - s1*s3^2*s4^2 + s2*s3^3*s4 - s2^2*s3*s4^2) - s3^5) - K^4*s4^5,
		- 5*K*(K^2*s3*s4^3 + K*(s1*s2*s3*s4 + s1*s3^3 - s1^2*s4^2 - s2^2*s3^2 + s2^3*s4) + s1^3*s2),
		- 5*K*(s1*s2^2 + s1^2*s3 + K*(s2*s4^2 + s3^2*s4)),
		- 5*K*(s1*s4 + s2*s3), 0, 1);
}
coef.p5K = function(b) {
	s1 = b[1]; s2 = b[2]; s3 = b[3]; s4 = b[4];
	c(
		K  = - 5*(s1*s4 + s2*s3),
		K2 = - 5*(s2*s4^2 + s3^2*s4), K = - 5*(s1*s2^2 + s1^2*s3),
		K3 = - 5*s3*s4^3, K2 = -5*(s1*s2*s3*s4 + s1*s3^3 - s1^2*s4^2 - s2^2*s3^2 + s2^3*s4),
			K = - 5*s1^3*s2,
		K4 = - s4^5,
			K3 = (5*(s1*s2*s4^3 - s1*s3^2*s4^2 + s2*s3^3*s4 - s2^2*s3*s4^2) - s3^5),
			K2 = (5*(s1*s2^3*s3 - s1^2*s2*s3^2 - s1^2*s2^2*s4 + s1^3*s3*s4) - s2^5), K = - s1^5
	);
}

### Other
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
# round to 0
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(any(isZero)) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(any(isZero)) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}

###
expand.grid.mod = function(p, n, K) {
	n = n - 1;
	lst = lapply(seq(n), function(id) seq(0, p-1));
	names(lst) = paste0("s", seq(n));
	lst$K = K;
	lst = expand.grid(lst);
	return(lst);
}

#####################

#################
### Exercises ###
#################


### E.1.) x^5 + 750*x + 3750
K = 24
#
m = unity(5, all=T)
k = K^(1/5)
x = sapply(m, function(m) { k = k*m; k^4/4 - k^3/2 - k^2/2 - k; } )
err = x^5 + 750*x + 3750;
round0(err)
# Note:
# x^5 + 6*5^3*x + 6*5^4; # NO 5^5!


### E.2.) x^5 - 25/2 * x^2 - 125/4
# - let P(x) = x^5 - 25/2 * x^2 - 125/4;
# - let r0 = s*k^4 - s*k^3 + k^2 + k, where k = K^(1/5)
# E.1.a.) Show that for K=2 and s=1/2, r0 is a root of P(x) = 0;
# E.1.b.) Show that r1 = s*(k*m)^4 - s*(k*m)^3 + (k*m)^2 + k*m,
#         is also a root of P(x) = 0, where m^5 = 1, m = root of unity;
# E.1.c.) Find the remaining roots of this polynomial;

K = 2
s = 1/2
#
k = K^(1/5)
x = s*k^4 - s*k^3 + k^2 + k
err = x^5 - 25/2 * x^2 - 125/4
err # small precision error
# Students are invited to provide the mathematical proof!


##################

### some Examples:
m = complex(re=cos(2*pi/5), im=sin(2*pi/5))
m = m^(0:4)

###
# Parameter: can be modified
K = 5
# Roots
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 5*K * x^3 - 5*(K^2 + K) * x^2 - 5*K^3 * x - K^4 - K^3 - K - 5*(K^3 - K^2)

# for K=2
K = 2 # fixed
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 10*x^3 - 30*x^2 - 40*x - 46

# for K=3
K = 3 # fixed
k = m * K^(1/5)
x = k^4 + k^3 + k
x^5 - 15*x^3 - 60*x^2 - 135*x - 201

###
K = -3 # can be modified
# test
k = if(K >= 0) K^(1/5) else  - (-K)^(1/5)
x = k^4 - 3*k^2 - k
x
x^5 + 5*K*x^3 + 15*K*(K + 3)*x^2 + 5*K*(28*K - 3)*x - K^4 + 15*K^3 + 198*K^2 + K


###
K = 2
s = 1/2
#
k = m*K^(1/5)
x = s*k^4 - s*k^3 + k^2 + k
#
x^5 - 5*K*(K*s^2 + K*s^3 - s + 1)*x^2 + 5*K*(K^2*s^4 + K*s^3 + 3*K*s^2 - K*s - 1)*x +
- K - K^2 - 10*K^2*s - 10*K^2*s^2 + 10*K^3*s^3 - 10*K^3*s^4 + K^3*s^5 - K^4*s^5
# for K = 2, s = 1/2
# x^5 - 25/2 * x^2 - 125/4


#################
#################

### Solution to Quintic of Class 1:

m = cos(2*pi/5) + 1i * sin(2*pi/5); m = m^seq(0, 4);

### Modular Arithmetic
K = 3; # arbitrary
s = c(-1,2,-3,1); # fixed
#
k = K^(1/5)
x = sapply(m, function(m) sum(s*(m*k)^seq(4)));
round0.p(poly.calc(x))
err = x^5 + 35*K*x^3 - 5*K*(11*K - 7)*x^2 + 5*K*(3*K^2 - 4*K + 2)*x - K^4 + 68*K^3 - 7*K^2 + K;
round0(err)

### Mod 2:
x = 1; # only Test;
(x^5 + 35*K*x^3 - 5*K*(11*K - 7)*x^2 + 5*K*(3*K^2 - 4*K + 2)*x - K^4 + 68*K^3 - 7*K^2 + K) %% 2
# r = k^4 + k^3 + k (mod 2)
(x^5 - 5*K*x^3 - 5*K*(K + 1)*x^2 - 5*K^3*x - K^4 - 6*K^3 + 5*K^2 - K) %% 2;
# which reduces to:
(x^5 + K*x^3 + K^3*x + K^4) %% 2; # (mod 2)

### Mod 3:
x = 1; # x = 2; # only Test;
(x^5 + 35*K*x^3 - 5*K*(11*K - 7)*x^2 + 5*K*(3*K^2 - 4*K + 2)*x - K^4 + 68*K^3 - 7*K^2 + K) %% 3;
# r = k^4 - k^2 - k (mod 3)
(x^5 + 5*K*x^3 + 5*K*(K + 1)*x^2 + 5*K*(2*K - 1)*x - K^4 + 5*K^3 - 4*K^2 + K) %% 3;
# which reduces to:
(x^5 - K*x^3 - K*(K + 1)*x^2 + K*(K + 1)*x + K^2) %% 3; # (mod 3)

### Mod 4:
p = 4;
s.lst = expand.grid.mod(p, n = 5, K = c(2,3,4,5));
pp = apply(s.lst, 1, function(x) { b = coef.p5(unlist(x[-5]), K = x[5]) %% p; b = rev(b); c(x, b); } )
pp = t(pp);
colnames(pp)[6:11] = paste0("x", seq(5,0,-1));

table(pp[, 8])
table(pp[, 9])
pp[ pp[,9] == 3, ]


#################

###
# for arbitrary K:
# (but polynomials need to be updated)
K = 3
s = c(1,-2,0,1)
#
k = m * K^(1/5)
x = sapply(k, function(k) sum(s*k^seq(4)));
round0.p(poly.calc(x))
-246 + 435*x + 30*x^2 - 15*x^3 + x^5
# P[^2, reduced] = P[r^2 - 2*...]
round0.p(poly.calc(x^2 - 2*(s[1]*s[4] + s[2]*s[3])*K))
866610 + 135405*x + 1440*x^2 + 735*x^3 + x^5
# P[^3, reduced] = P[r^3 - 3*...]
round0.p(poly.calc(x^3 - 3*((s3^2*s4 + s2*s4^2)*K^2 + (s1*s2^2 + s1^2*s3)*K)))
-1474293000 + 78324440*x - 134730*x^2 + 15660*x^3 + x^5

# We can construct an overdetermined polynomial system:
s1 = s[1]; s2 = s[2]; s3 = s[3]; s4 = s[4];
# Primary coefficients:
# b3/5 == b4[P[^2]] / 10
(s[1]*s[4] + s[2]*s[3])*K
# - b2/5 == b4[P[^3]] / 15
(s3^2*s4 + s2*s4^2)*K^2 + (s1*s2^2 + s1^2*s3)*K
# ...
# Higher coefficients of P[^n]:
s1^3*s2*K + (s1*s3^3 + s2^3*s4 + 3*prod(s))*K^2 +
	+ s3*s4^3*K^3 # = (b3[*] - b2^2/5)/5
(-735 - 15^2/5) / 5 / 2

# TODO: ???
# - but very hard to crack!

s1 = s[1]; s2 = s[2]; s3 = s[3]; s4 = s[4];
# - b3
((s1*s2^2+s1^2*s3)*K + (s3^2*s4+s2*s4^2)*K^2) * 5
# - b3[*]
(2*s1^3*s2*K + ((s1*s4 + s2*s3)^2 + 2*s1*s3^3 + 2*s2^3*s4 + 6*prod(s))*K^2 +
	+ 2*s3*s4^3*K^3) * 5

###
r = toPoly.pm("s4*k^4 + s3*k^3 + s2*k^2 + s1*k")
# P[r]
p = diff.pm(r*r, toPoly.pm("5*(s[1]*s[4] + s[2]*s[3])*k^5"))
p = mult.pm(p, r)
p3 = p[p$k %% 5 == 0, ]
p3$coeff = p3$coeff / 3;
p3$k = p3$k / 5; names(p3)[1] = "K";
print.pm(p3)
# P[r^3]
p = pow.pm(r, 3)
p3 = p[p$k %% 5 == 0, ]
p3$k = p3$k / 5; names(p3)[1] = "K";
eval.pm(p3, c(K, s)) * -5


####################

### TODO:
# - clean;
# - expand;

library(polynom)

# alternative function to compute coefficients of polynomial
elemPoly = function(x, start=0, adjustSign=TRUE) {
	coeff = sapply(c(start, seq_along(x)), function(n){
		sum(apply(combn(x, n), 2, prod))
	})
	if(adjustSign) {
		if(start %% 2 ==0) {
			adj = c(1,-1)
		} else {
			adj = c(-1,1)
		}
		len = length(coeff)
		adj = rep(adj, len/2)
		if(len %% 2 == 1) {
			adj = c(adj, adj[1])
		}
		coeff[abs(coeff) < 1E-10 ] = 0
		coeff = coeff * adj
	}
	return(coeff)
}

poly.radical.create = function(root_coeff, K, n=length(root_coeff), doRound=TRUE) {
	k = ifelse(K < 0,
			ifelse(n %% 2 == 0, complex(re=0, im=(-K)^(1/n)), -(-K)^(1/n) ),
			K^(1/n))
	len = length(root_coeff) - 1
	n_1 = n - 1
	# root rotation
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	### Roots
	x = sapply(0:n_1, function(i) sum(root_coeff * (k*m^i)^(0:len)) )
	
	### Polynomial
	# does give a warning with complex roots
	# poly = poly.calc(x)
	poly = elemPoly(x) # * (-1)^((n:0)+ifelse(n %% 2 == 0, 0, 1) )
	# round: dangerous, but this is a minimal polynomial with integer coeffs,
	# if K and coeff are integers
	if(doRound) {
		poly = round(poly)
	}
	poly.coeff = as.numeric(round(poly))
	poly.str = paste(poly.coeff, "*x^", n:0, sep="")
	poly.str = poly.str[poly.coeff != 0]
	
	poly.list = list(
		r = x,
		poly.coeff = poly,
		poly = paste(poly.str, collapse="+"),
		n = n
	)
	return(poly.list)
}

### Examples

# E.1.) root = k^4 - 2*k^2 - k, where k = K^(1/5)
K = 2
p = poly.radical.create(c(0, -1, -2, 0, 1), K, n=5)
p
# 114 + 160*x + 80*x^2 + 10*x^3 + x^5
err = sapply(1:p$n, function(id) sum(p$poly * p$r[id]^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)


# E.2.) root = k^5 - k^3 - 3*k^2 + k, where k = K^(1/6)
K = 2
p = poly.radical.create(c(0, 1, -3, -1, 0, 1), K, n=6)
p
# 3726 + 2268*x + 594*x^2 + 108*x^3 - 18*x^4 + x^6
err = sapply(p$r, function(root) sum(p$poly * root^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)


# E.3.) root = k^6 - k^3 - 3*k^2 + k, where k = K^(1/7)
### NOTE: overflows easily with loss of precision!
K = 2
p = poly.radical.create(c(0, 1, -2, -1, 0, 0, 1), K, n=7)
p
# 818 + 672*x - 126*x^2 + 308*x^3 + 98*x^4 - 14*x^5 + x^7
err = sapply(p$r, function(root) sum(p$poly * root^(p$n:0)) )
err
# with package polynom
p.p = polynomial(rev(p$poly))
p.p
predict(p.p, p$r)


#####################
#####################

#####################
### Entanglements ###
#####################

### Multiple Radicals

#######################
### A.) Simple Radicals

m3 = unity(3, all=T)
m5 = unity(5, all=T)
m = expand.grid(m3, m5)
colnames(m) = c("m3", "m5")
# TODO: use grid

### Entanglement with s0
K1 = 2
K2 = 3
#
k1 = K1^(1/3)*m3
k2 = K2^(1/5)*m5
#
s0 = k1^2 - k1
x = sapply(s0, function(s0) k2^3 - k2 + s0)
#
-1833890 + 2488620*x - 655200*x^2 - 283735*x^3 + 322905*x^4 - 85131*x^5 + 22195*x^6 + 3015*x^7 + 
+ 1075*x^9 - 582*x^10 + 495*x^11 - 55*x^12 + 30*x^13 + x^15


### Entanglement with s1
x = sapply(s0, function(s0) k2^3 - k2^2 + s0*k2)
-5706018 - 260010*x + 2972700*x^2 + 3948615*x^3 + 3385800*x^4 + 1893861*x^5 + 650700*x^6 +  
+ 279990*x^7 + 67230*x^8 + 18225*x^9 + 4986*x^10 + 900*x^11 + 180*x^12 + 45*x^13 + x^15


### Entanglement with multiple s[j]
x = sapply(k1, function(k) k2^3 + (k^2 - k) * k2^2 + (k-1)*k2)
-4067604 - 10811880*x - 11835720*x^2 - 8053695*x^3 - 4910625*x^4 - 2417553*x^5 - 619515*x^6 +  
44145*x^7 + 109890*x^8 + 54675*x^9 + 12708*x^10 + 405*x^11 - 315*x^12 + x^15


#######################
### A.) Nested Radicals

n1 = 5
n2 = 3
m5 = unity(n1, all=T)
m3 = unity(n2, all=T)

### 3x5: f2( ( f1(inner_root^1/5) )^1/3 )
K = 2
s1 = c(0, 1, -2, 2) # coeffs of nested root
s2 = c(0, -1, 1) # coeffs of outer root
#
id1 = 0:(length(s1)-1)
id2 = 0:(length(s2)-1)
k = K^(1/n1) * m5 # all radical units
r = sapply(k, function(k) sum(s1 * k^id1) )
r = as.vector(r^(1/n2)) # all nested roots
### all outer roots
r = sapply(r, function(x) x * m3) # base units for outer roots
r = as.vector(r)
x = sapply(r, function(x) sum(s2 * x^id2)) # actual roots
#
-37410 + 160950*x + 182700*x^2 + 421930*x^3 + 247860*x^4 + 184140*x^5 + 74020*x^6 +
+ 30990*x^7 + 8100*x^8 + 3800*x^9 + 780*x^10 + 360*x^11 + 80*x^12 + x^15


### 5x3: f2( ( f1(inner_root^1/3) )^1/5 )
K = 2
s1 = c(1, -2, 1) # coeffs of nested root
s2 = c(0, -1, -1, 2) # coeffs of outer root
#
id1 = 0:(length(s1)-1)
id2 = 0:(length(s2)-1)
k = K^(1/n2) * m3
r = sapply(k, function(k) sum(s1 * k^id1) )
r = as.vector(r^(1/n1))
# all roots
r = sapply(r, function(x) x * m5)
r = as.vector(r)
x = sapply(r, function(x) sum(s2 * x^id2))
#
578700 - 279900*x - 807900*x^2 + 2096575*x^3 - 563775*x^4 + 3355995*x^5 - 412805*x^6 +  
+ 735885*x^7 + 101730*x^8 + 25075*x^9 + 1212*x^10 + 225*x^11 - 15*x^12 + 30*x^13 + x^15



##############################
### B.) Triple Nested Radicals

library(polynom)

n1 = 3
n2 = 3
n3 = 3
m_1 = unity(n1, all=T)
m_2 = unity(n2, all=T)
m_3 = unity(n3, all=T)

### 3x3x3: f3( f2( f1(inner_root^1/3)^1/3 )^(1/3) )
K = 2
# !!! Overflows easily !!!
s1 = c(0, -1, 1) # coeffs of inner root: s0 + s1*k + s2*k^2;
s2 = c(0, -1, 1) # coeffs of middle root
s3 = c(0, -1, 1) # coeffs of outer root
#
id1 = 0:(length(s1)-1)
id2 = 0:(length(s2)-1)
id3 = 0:(length(s3)-1)
# for ALL n1*n2*n3 roots:
m.gr = expand.grid(m_1, m_2, m_3)

# Level 1 (innermost)
k = K^(1/n1) * m.gr[,1]
x = sapply(k, function(k) sum(s1*k^id1))
# Level 2 (middle)
k = x^(1/n2) * m.gr[,2]
x = sapply(k, function(k) sum(s2*k^id2))
# Level 3 (outer)
k = x^(1/n3) * m.gr[,3]
x = sapply(k, function(k) sum(s3*k^id3))

### Polynomial
round0.p(poly.calc(x))

# TODO: check if correct [using a robust approach]
# - coefficients overflow easily!
err = -3410 - 31050*x - 129600*x^2 - 311076*x^3 - 453222*x^4 - 402894*x^5 - 270270*x^6 - 250614*x^7 +
- 90882*x^8 + 509070*x^9 + 992844*x^10 + 732240*x^11 + 141462*x^12 - 56754*x^13 + 66582*x^14 +
+ 108162*x^15 + 32238*x^16 - 9234*x^17 - 1464*x^18 + 4752*x^19 + 1620*x^20 - 252*x^21 - 108*x^22 + x^27
round0(err) # some rounding errors in evaluation!



#####################
#####################


### Experimental
# NOT yet functional

m5 = unity(5, all=T)

K = 2
n = 3
#
k = (n^5 + K)^(1/5)*m5
b0 = k - n

p = mult.p(c(b0[2]*b0[5], -k[2]-k[5], 1), c(b0[3]*b0[4], -k[3]-k[4], 1))
p = round0.p(mult.p(p, c(-b0[1], 1)))
p



