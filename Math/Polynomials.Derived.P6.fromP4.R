
########################
###
### Leonard Mada
### [the one and only]
###
### P6 Polynomials:
### Derived from Root Permutations
###
### draft v.0.1f


### Generate P6
# - by Permuting & Entangling roots,
#   e.g. the roots of P4 or P3;
# - P3: new roots = f(r[i], r[j]) based on the 6 permutations;
# - P4: new roots = (r[i] - r[j])^2,
#   and using only the 6 distinct values;
# - deriving a P9 from composition of 2 P3:
#   P3(P3) => P9 = P3[base] * P6
#   and factorization of P9 using initial/base P3;
# - this will be discussed more extensively in:
#   Polynomials.Derived.R;

### TODO:
# - compute elementary polynomials;


#####################

###############
### History ###

# draft v.0.1d-v.0.1f:
# - added some derivations of type P3(P3):
#   including more examples;
# draft v.0.1c:
# - added the P3 permutations;
#   [initially in Polynomials.Derived.P6.R]
# draft v.0.1b:
# - added various special cases:
#   these are however decomposable into P[2]*P[4] or P[3]*P[3];
# draft v.0.1a:
# - moved from Polynomials.Derived.P6.R [v.0.3e [pre-z]]
#   to separate file;


################
### Sections ###

### Section A:
# - P4 Permutations;
### Section B:
# - P3 Permutations;
### Section C:
# - Derived P6 from P3;


####################

library(polynom)
library(pracma)

### helper functions
unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
### Generator
permute.p = function(coeff) {
	r = roots(coeff)
	r.g = expand.grid(r, r)
	r.g = r.g[ r.g[,1] != r.g[,2] , ]
	x = round0((r.g[,1] - r.g[,2])^2)
	x = unique(x) # TODO: may fail;
	p = round0.p(poly.calc(x))
	return(list(x=x, p=p))
}

###################

### P4 => P6

###
coeff = c(1,0,0,1, 1)
p = permute.p(coeff)
x = p$x
229 + 216*x - 112*x^2 + 26*x^3 + 8*x^4 + x^6

# TODO:
# x^6 + (-3*E1^2 + 8*E2)*x^5 + (...)*x^4 + ...
#
# b4 = (3*r^4 - 4*r[i]*r[j]^3 + 8*r[i]^2*r[j]^2 - 2*r[i]*r[j]*r[k]^2 + 12*E4)
#
# (3+4)*(E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4)
# - 4*E1*(E1^3 - 3*E1*E2 + 3*E3) +
# + 8*r[i]^2*r[j]^2 + # <== TODO
# - 2*r112 + 12*E4;

### Test/Derivation
E = coeff[-1]
r2 = E[1]^2 - 2*E[2];
r3 = -(E[1]^3 - 3*E[1]*E[2] + 3*E[3]);
r4 = E[1]^4 - 4*E[1]^2*E[2] + 4*E[1]*E[3] + 2*E[2]^2 - 4*E[4];
r13 = E[1]*(E[1]^3 - 3*E[1]*E[2] + 3*E[3]) - (E[1]^4 - 4*E[1]^2*E[2] + 4*E[1]*E[3] + 2*E[2]^2 - 4*E[4]);
r112 = E[2] * r2 - r13
# r[i]*r[j]*r[k]^2 = E2 * r^2 - r[i]*r[j]^3
#   = E2*(E1^2 - 2*E2) - E1*(E1^3 - 3*E1*E2 + 3*E3) + (E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4);
# r[i]*r[j]^3 = E1*(E1^3 - 3*E1*E2 + 3*E3) - (E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4);
# r^4 = E1^4 - 4*E1^2*E2 + 4*E1*E3 + 2*E2^2 - 4*E4;
# r^3 = -(E1^3 - 3*E1*E2 + 3*E3);
# r^2 = E1^2 - 2*E2;
#
# pow = 3
# sum(r4[1]*r4[-1]^pow, r4[2]*r4[-2]^pow, r4[3]*r4[-3]^pow, r4[4]*r4[-4]^pow)
#
# b3 = (-r^6 + 2*r[i]*r[j]^5 - 7*r[i]^2*r[j]^4 + 8*r[i]^3*r[j]^3 + 2*r[i]*r[j]^2*r[k]^3 + 16*r[ijk]*r[m]^3 +
#  -10*r[i]^2*r[j]^2*r[k]^2)*x^3 + ...

# Derivation:
# r^4 = E1*(E1^3 - 3*E1*E2 + 3*E3) - E2*(E1^2 - 2*E2) + E1*E3 - 4*E4
# (x-(r1-r2)^2)*(x-(r1-r3)^2)*(x-(r1-r4)^2)*(x-(r2-r3)^2)*(x-(r2-r4)^2)*(x-(r3-r4)^2)
# (x-(a-b)^2)*(x-(a-c)^2)*(x-(a-d)^2)*(x-(b-c)^2)*(x-(b-d)^2)*(x-(c-d)^2)
r.g = expand.grid(1:4,1:4,1:4) # doubled: r.g[,2] vs r.g[,3]
r.g = r.g[ r.g[,1] != r.g[,2] , ]
r.g = r.g[ r.g[,1] != r.g[,3] , ]
r.g = r.g[ r.g[,2] != r.g[,3] , ]
sum(r[r.g[,1]]^2*r[r.g[,2]]*r[r.g[,3]]) / 2

###
coeff = c(1,0,0,1,-1)
p = permute.p(coeff)
x = p$x
-283 - 216*x - 112*x^2 + 26*x^3 - 8*x^4 + x^6

# and a inverse-symmetric P12 for fun:
x = sapply(x, function(r) roots(c(1, 1-r, -1))) # parameters: +/- r & b0 = +/- 1;
1 - 6*x + x^2 + 16*x^3 - 80*x^4 + 370*x^5 - 436*x^6 - 370*x^7 - 80*x^8 - 16*x^9 +
+ x^10 + 6*x^11 + x^12


###
coeff = c(1,0,1,1, 2)
p = permute.p(coeff)
x = p$x
1825 - 250*x - 335*x^2 + 86*x^3 + 38*x^4 + 8*x^5 + x^6

###
p = sapply(-6:6, function(b) print(permute.p(c(1, 0,0,b,1))$p))


p = sapply(-6:6, function(b) print(permute.p(c(1, 0,0,-b,b))$p))


#########
### P12 !
coeffs = c(1,2,0,0, 1)
alpha = 2
r4 = roots(coeffs)
r.g = expand.grid(r4, r4)
x = round0((alpha*r.g[,1] - r.g[,2]) )
x = unique(x[ x != 0 ])
round0.p(poly.calc(x)) / polynomial(rev(coeffs))
# a P12 when alpha != 1


#####################
### Special Cases ###

# all specia cases are decomposable:
# => P[2]*P[4] or
# => P[3]*P[3];


### m5 => P[2]*P[4]
m5 = unity(5, all=T)
m5 = m5[-1]

x.r = m5^3 - 2*m5^2 - m5
r.gr = expand.grid(x.r, x.r)
x = round0((r.gr[,1] - r.gr[,2])^2)
x = x[ x != 0]
15125 + 152500*x + 396700*x^2 + 62250*x^3 + 3740*x^4 + 100*x^5 + x^6


### ()^(1/4) => P[2]*P[4]
m4 = unity(4, all=T)
K = 2
k = K^(1/4) * m4

x.r = k^3 + k^2 + 3*k
r.gr = expand.grid(x.r, x.r)
x = round0((r.gr[,1] - r.gr[,2])^2)
x = x[ x != 0]
-658409500 + 11069440*x - 4889024*x^2 - 349696*x^3 + 16112*x^4 - 224*x^5 + x^6


m3 = unity(3, all=T)
m3diff = cbind((1-m3[-1]), (1-1/m3[-1]))
m3diff = rbind(m3diff, c(m3[2]-m3[3], m3[3]-m3[2]))
p6fromP3.gen = function(s, c, d) {
	det = sqrt(d^2 - c^3 + 0i)
	p.r = (d + det)^(1/3); q.r = (d - det)^(1/3)
	r = p.r * m3 + q.r / m3 - s
	r = c(r, p.r*m3diff[,1] + q.r*m3diff[,2])
	x = r^2
	p = poly.calc(x)
	return(list(x=x, p=p))
}

# TODO: parametric P6

### P[1]*P[3] => P[3]*P[3]
r4 = roots(c(1,2,0,0, 1))
r.g = expand.grid(r4, r4)
x = round0((r.g[,1] - r.g[,2])^2)
x = unique(x[ x != 0 ])
round0.p(poly.calc(x))
-176 + 288*x - 16*x^2 - 88*x^3 + 56*x^4 - 12*x^5 + x^6


### (x - s) * P3 => P[3]*P[3]
d = 2
c = 1
s = 2 # shift / P[1]
#
p = p6fromP3.gen(s,c,d)
p
x = p$x
-1296 + 33696*x + 2745*x^2 - 3028*x^3 + 510*x^4 - 36*x^5 + x^6



### (x - s) * P3 => P[3]*P[3]
d = 2
c = 1
s = -1 # shift / P[1]
#
p = p6fromP3.gen(s,c,d)
p
x = p$x
-1296 - 4212*x - 3816*x^2 - 193*x^3 + 231*x^4 - 27*x^5 + x^6
(-4 - 12*x - 9*x^2 + x^3) * (324 + 81*x - 18*x^2 + x^3)




##########################
##########################

##########################
### Permutations of P3 ###


library(polynom)
library(pracma)

# generate polynomial based on base roots
# - added also r.g[,2];
fromP3.gen = function(coeff, s) {
	# coeffs: in descending order (+ leading coeff)
	r = roots(coeff) # base P3
	r.g = expand.grid(r, r)
	r.g = r.g[ r.g[,1] != r.g[,2] , ]
	if(nrow(r.g) == 0) return(list(x=NA, p=NA))
	# new roots
	r.g = cbind(r.g[,1], r.g[,2], r.g[,1]*r.g[,2], r.g[,2]^2, r.g[,1]*r.g[,2]^2)
	r.g = r.g[ , 1:length(s)]
	x = sapply(1:nrow(r.g), function(id) round0(sum(s * r.g[id,])))
	p = round0.p(poly.calc(x))
	return(list(x=x, p=p))
}


### from P3 roots: some P6
s = c(1, 0, 1, 0, 0)
coeff = c(1,0,1,1)
p = fromP3.gen(coeff, s)
x = p$x
p
1 - x + 6*x^2 + 3*x^3 - 2*x^5 + x^6

###
s = c(1, 0, 1, 0, -1)
coeff = c(1,0,1,1)
p = fromP3.gen(coeff, s)
x = p$x
p
1 + 4*x + 34*x^2 + 27*x^3 + 4*x^4 + x^5 + x^6

###
s = c(2, 0, -1, 0, 1)
coeff = c(1,-1, 1, 1)
p = fromP3.gen(coeff, s)
x = p$x
p
4 + 4*x - 24*x^3 + 22*x^4 - 6*x^5 + x^6


###
p = sapply(-6:6, function(b) print(fromP3.gen(c(1,b,1,1), c(1,0,1,-1))$p))

p = sapply(-6:6, function(b) print(fromP3.gen(c(1,b,1,1), c(1,0,1,-b))$p))

p = sapply(-6:6, function(b) print(fromP3.gen(c(1,b,1,1), c(1,-1,1,0))$p))

p = sapply(-6:6, function(b) print(fromP3.gen(c(1,0,b,b), c(1,-1,1,0))$p))


##########################
##########################

#################
### Section C ###
#################

### P6 Derivations from P3;

### Variant 1:
# - let B[3] be an initial/baseline polynomial of order 3 with roots r[b];
# - let D[3] be another polynomial:
#  -- the roots of D(x) - D(r) will be the roots of a P[9];
#  -- P[9] = B[3] * P[6];

### Variant 2:
# [viewing Polynomials with the eyes of Lagrange and the skills of Gauss]
# - we can use the structure of B[3] and incorporate it into D[3];


### Example 1:
K = 2
coeff = c(1,1,1,K)
r = roots(coeff)
r.der = -(coeff[3] * r + K); # r^3 + r^2
x = sapply(1:3, function(id) roots(c(1,0, r[id], -r.der[id])))
p1 = round0.p(poly.calc(x))
p = p1 / round0.p(poly.calc(r))
p1; p;
err = 2 - 4*x - x^2 + 5*x^3 - x^4 - x^5 + x^6
round0(err)


### Example 2:
K = 1
coeff = c(1,1,2,K)
r = roots(coeff)
r.der = -(coeff[3] * r + K); # r^3 + r^2
x = sapply(1:3, function(id) roots(c(1,0, r[id], -r.der[id])))
p1 = round0.p(poly.calc(x))
p = p1 / round0.p(poly.calc(r))
p1; p;
err = -1 - 3*x + 3*x^2 + 4*x^3 - 2*x^4 - x^5 + x^6
round0(err)



### Example 3:
K = 1
coeff = c(1,1,2,K)
r = roots(coeff)
r.der = -(coeff[3] * r + K); # r^3 + r^2
x = sapply(1:3, function(id) roots(c(1, r[id], 2*r[id], -2 * r.der[id])))
p1 = round0.p(poly.calc(x))
p = p1 / round0.p(poly.calc(r))
p1; p;
err = -8 - 24*x + 4*x^2 + 12*x^3 - 2*x^5 + x^6
round0(err)


p3der_test.gen = function(coeff) {
	r = roots(coeff)
	r.der = -(coeff[3] * r + coeff[4]); # r^3 + r^2
	x = sapply(1:3, function(id) roots(c(1, r[id], 2*r[id], -2 * r.der[id])))
	p1 = round0.p(poly.calc(x))
	p = p1 / round0.p(poly.calc(r))
	return(list(x=x, p=round0.p(p), p1=p1))
}
###
K = 1
p = sapply(-6:6, function(s) print(p3der_test.gen(c(1,1,s,K))$p))
#
s = -1
p = p3der_test.gen(c(1,1,s,K))
x = p$x
err = 16 - 8*x^2 - 2*x^5 + x^6
round0(err)
x



### simple Examples / Introductory Examples
### x^3 = x + 1

p3derived.gen = function(K, coeff=c(1)) {
	r = roots(c(1,0,coeff,K))
	x.r = -(r * coeff[1] + K)
	x = sapply(x.r, function(r) roots(c(1,0,0,-r)))
	p1 = round0.p(poly.calc(x))
	p = round0.p(p1 / round0.p(poly.calc(r)))
	return(list(x=x, p=p))
}

###
p = p3derived.gen(1)
p
x = p$x
round0(1 - x + x^2 + 2*x^3 - x^4 + x^6)


###
p = sapply(-6:6, function(s) print(p3derived.gen(s)$p))
p = sapply(-6:6, function(s) print(p3derived.gen(1, coeff=c(s))$p))
p = sapply(-6:6, function(s) print(p3derived.gen(s, coeff=c(s))$p))
#
s = -4
p = p3derived.gen(1, coeff=c(s)); x = p$x; p
round0(1 + 4*x + 16*x^2 + 2*x^3 + 4*x^4 + x^6)

###
K = 4
x = p3derived.gen(K, coeff=c(K))$x
round0( K^2 - K^2*x + K^2*x^2 + 2*K*x^3 - K*x^4 + x^6 )
round0( (x^3 - K*x + K)^2 + K*x^4 + K^2*x )



