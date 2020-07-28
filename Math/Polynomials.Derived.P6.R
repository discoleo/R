
########################
###
### Leonard Mada
### [the one and only]
###
### P6 Polynomials
### Derived from Special Factorizations
###
### draft v.0.3-pre-alpha


### Factorization of the P6 Polynomials
### into "conjugated" P2/P3 Polynomials
### with NON-Rational Coefficients


### Introduction

# P6 polynomials can be factorized into:
# 1.) P3 * P3[*]
#  - where P3 & P3[*] are conjugate polynomials;
#  Variants:
#  - Coefficients are based on:
#    (a +/- b*sqrt(n)), cos(2*pi/5),
#    or (m, m^2) where m^3 = 1;
# 2.) P2 * P2[*] * P2[**]
#  - where P2, P2[*], P2[**] are conjugate polynomials;
# - these are trivial or "almost"-trivial factorizations,
#   but they comprise a large sub-family of P6 polynomials;
# 3.) Cubic => sqrt(roots) => Coeff of Quadratic
#  - generates a large subfamily of special/symmetric P6 polynomials;
#  - based on transformation of the roots of a base cubic;

# Formulas are also provided to generate
# a large sub-family of such P6 polynomials
# with b0 = 1 and remaining b[j] = integers;

# Note:
# Even the Class 1 polynomials of order 6
# can be factorized into 2 conjugate P3's:
# P3 * P3[*], where the coefficients are of the form:
# b[j] = b[j,0] +/- b[j,1] * sqrt(K);
# e.g. roots = k^3 - k^2 - k, with k=K^(1/6):
# 62 + 84*x + 6*x^2 - 20*x^3 - 6*x^4 + x^6
#  (8-sqrt(2) + (6-3*sqrt(2))*x - 3*sqrt(2)*x^2 + x^3)*
#  (8+sqrt(2) + (6+3*sqrt(2))*x + 3*sqrt(2)*x^2 + x^3)


#####################

### History
# draft v.0.3-pre-alpha:
# - initial workout of the cubic => Sqrt technique;
# - large sub-family of special P6 polynomials;
# - entire families like:
#   1 - x - b*x^3 - x^5 + x^6 = 0;
#   1 - x + b2*x^2 - b3*x^3 + b2*x^4 - x^5 + x^6 = 0;
# draft v.0.2f-g:
# - Solution & mechanism for:
#   1 + x - x^2 - x^4 + x^5 + x^6 = 0;
#   1 - 3*x - 3*x^5 + x^6 = 0;
#   1 - 3*x + 2*x^3 - 3*x^5 + x^6 = 0;
#   [at the end of this document + more examples]
# draft v.0.2c-e:
# - polynomials having b0 = 1;
# - formula for sqrt entanglement:
#   b0 = (n +/- sqrt(n^2 - 4)) / 2;
# - special case for n = 1:
#   b0 = (1 +/- 1i*sqrt(3))/2;
# - formula for generic sqrt(n) entanglement:
#   b0 = 1 => prod(b0) = 1;
# - various tests/examples for the sqrt formulas;
# draft v.0.2b:
# - 1st formula for sqrt entanglement:
#   b0 = (n +/- sqrt(n^2 - 1));
# draft v.0.2:
# - formulas for the cos(2*pi/5) entanglement;
# draft v.0.1b - v.0.1f:
# - entanglements with:
#   m[3] (m^3 = 1), cos(2*pi/5), cos(2*pi/7);
# - explicit formula for the cos(2*pi/5) entanglement (in v.0.1f);
# - new variants:
#   3 - 3*x + x^6 = 0;
#     factored as 2 cubics:
#     (x^3 + (m-m^2)*x^2 +(m+2*m^2)*x - 2*m-m^2)*
#     (x^3 - (m-m^2)*x^2 +(m^2+2*m)*x - 2*m^2-m);
#   
#   many other, including:
#   [the nice ones tend to be actually trivial]
#   1 + 5*x - 5*x^2 + 3*x^3 + x^6 = 0; # both +/- 5*x: -5*(x+1/2)^2 + (x^3+3/2)^2
#   1 - 5*x - 5*x^2 + 3*x^3 + x^6 = 0; # also as cos(2*pi/5) entanglements
#   1 + 5*x + 5*x^2 - 3*x^3 + x^6 = 0;
#   1 + 5*x + 5*x^2 + 7*x^3 + x^6 = 0;
#   1 + 2*x^2 + 2*x^3 + x^6 = 0 # trivial, to complete the series;
#   1 - 2*x^2 + 2*x^3 + x^6 = 0 # see parametric variants: b2 = -2 + 0i; b3 = 2;
#   -1 - x + 2*x^3 - 2*x^5 + x^6 = 0;
#   1 + 3*x^2 + b*x^3 + 3*x^4 + x^6 = 0; # trivial, just as completion;
# - additional practice with polynomials;
# draft v.0.1
# - moved P6 from Polynomials.Derived.R
#   in a separate file;
# - old version (up to draft v.0.3f):
#   https://github.com/discoleo/R/blob/master/Math/Polynomials.Derived.R
# v.0.3.c - v.0.3.f [in the initial file]:
# - more awesome polynomials with complete solutions:
#   1 + x - x^4 - x^5 + x^6 = 0;
#   1 - x + x^2 + x^3 + x^4 + x^6 = 0;
#   1 + x + x^2 - x^3 - 2*x^4 + x^6 = 0;
#   1 + 2*x + x^2 - x^3 - x^4 + x^6 = 0;
#   & many more;


################

##############
### Theory ###

# TODO:

# P3 or P2 polynomials can be entangled
# to generate very diverse P6 outputs;

### Entanglements
# a.) with Roots of Unity:
# e.g. (x^3 - m*x^2 - x + 1)*(x^3 - m^2*x^2 - x + 1)
# x^6 + x^5 - x^4 + x^3 + 2*x^2 - 2*x + 1
#
# b.) with 2*cos(2*pi/5):
# (x^3 - 2*cos(2*pi/5)*x^2 + 1)*(x^3 - 2*cos(4*pi/5)*x^2 + 1)
# 1 + x - x^2 + 2*x^3 + x^4 + x^6
# Example 2:
# (x^3 + 2*cos(2*pi/5)*x^2 - 1)*(x^3 + 2*cos(4*pi/5)*x^2 - 1)
# 1 + x^2 - 2*x^3 - x^4 - x^5 + x^6
#
# c.) with 2*cos(2*pi/7):
# c3 = 2*cos( 1:3 * 2*pi/7)
# (c3[1]*x^2 - x + 1)*(c3[2]*x^2 + x + 1)*(c3[3]*x^2 - x + 1)
# x^6 + 2*x^5 - 3*x^4 + x^3 + 2*x^2 - 3*x + 1
#
# d.) Roots of Cubic with +/- sqrt(r):
# Cubic => {sqrt(r), -sqrt(r)}
# then: x^2 - r.sq * x + b0 = 0
# generates a large sub-family of P6 polynomials:
sapply(-6:6, function(b) print(p6sq.gen(c(1,5, 4, b))$p))
1 - x - x^2 - 4*x^3 - x^4 - x^5 + x^6
1 - x - x^2 - 3*x^3 - x^4 - x^5 + x^6
1 - x - x^2 - 2*x^3 - x^4 - x^5 + x^6
1 - x - x^2 - x^3 - x^4 - x^5 + x^6
1 - x - x^2 - 0 - x^4 - x^5 + x^6
1 - x - x^2 + x^3 - x^4 - x^5 + x^6
1 - x - x^2 + 2*x^3 - x^4 - x^5 + x^6
1 - x - x^2 + 3*x^3 - x^4 - x^5 + x^6
1 - x - x^2 + 4*x^3 - x^4 - x^5 + x^6
# p6sq.gen(c(1,5, 5, b)): generates the 0*x^2 + 0*x^4 variants;
# and the + x^2 variants:
sapply(-6:6, function(b) print(p6sq.gen(c(1,5, 6, b))$p))
1 - x + x^2 - 4*x^3 + x^4 - x^5 + x^6
1 - x + x^2 - 3*x^3 + x^4 - x^5 + x^6
1 - x + x^2 - 2*x^3 + x^4 - x^5 + x^6
1 - x + x^2 - x^3 + x^4 - x^5 + x^6
1 - x + x^2 - 0 + x^4 - x^5 + x^6
1 - x + x^2 + x^3 + x^4 - x^5 + x^6
1 - x + x^2 + 2*x^3 + x^4 - x^5 + x^6
1 - x + x^2 + 3*x^3 + x^4 - x^5 + x^6
1 - x + x^2 + 4*x^3 + x^4 - x^5 + x^6


####################

### TODO:
# - use exact P3 solver for P6, when applicable;
library(polynom)
library(pracma)

m = unity(3, all=FALSE) # load first helper function unity()!

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
### Generator: Convolved P3
# with b0 = 1 (for the 5nnn and 2xnnn variants);
solve.p6 = function(coeff, type=110) {
	m = unity(3, all=FALSE)
	a = coeff[1]; b = coeff[2]; c = coeff[3]
	# TODO: all permutations
	# TODO: automate permutations;
	if(type == 110) {
		# (x^3 + a*m*x^2 + b*m*x + c)*(x^3 + a*m^2*x^2 + b*m^2*x + c)
		x = c(roots(c(1,a*m,b*m,c)), roots(c(1,a*m^2,b*m^2,c)))
		coeff = c(1, - a, (a^2-b), 2*(a*b+c), (b^2-a*c), - b*c, c^2)
		err = x^6 - a*x^5 + (a^2-b)*x^4 + 2*(a*b+c)*x^3 + (b^2-a*c)*x^2 - b*c*x + c^2
	} else if(type == 111) {
		x = c(roots(c(1,a*m,b*m,c*m)), roots(c(1,a*m^2,b*m^2,c*m^2)))
		coeff = c(1, - a, (a^2-b), (2*a*b-c), (b^2+2*a*c), 2*b*c, c^2)
		err = x^6 - a*x^5 + (a^2-b)*x^4 + (2*a*b-c)*x^3 + (b^2+2*a*c)*x^2 + 2*b*c*x + c^2
	} else if(type == 1) {
		x = c(roots(c(1,a,b,c*m)), roots(c(1,a,b,c*m^2)))
		coeff = c(1, 2*a, (a^2+2*b), (2*a*b-c), (b^2-a*c), -b*c, c^2)
		err = x^6 + 2*a*x^5 + (a^2+2*b)*x^4 + (2*a*b-c)*x^3 + (b^2-a*c)*x^2 - b*c*x + c^2
	} else if(type == 222) {
		# the full version: (c_1 * m + c_2 * m^2); (c0 is not needed)
		a1 = coeff[1]; a2 = coeff[2]
		b1 = coeff[3]; b2 = coeff[4]
		c1 = coeff[5]; c2 = coeff[6]
		x = c(
			roots(c(1, coeff[1]*m+coeff[2]*m^2, coeff[3]*m+coeff[4]*m^2, coeff[5]*m+coeff[6]*m^2)),
			roots(c(1, coeff[1]*m^2+coeff[2]*m, coeff[3]*m^2+coeff[4]*m, coeff[5]*m^2+coeff[6]*m)))
		coeff = c(1, - (a1+a2), (a1^2+a2^2-a1*a2-b1-b2), (2*a1*b1+2*a2*b2-a1*b2-a2*b1-c1-c2),
			(b1^2+b2^2-b1*b2+2*a1*c1+2*a2*c2-a1*c2-a2*c1),
			(2*b1*c1+2*b2*c2-b1*c2-b2*c1), c1^2+c2^2-c1*c2)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == 5229) {
		# TODO: rename to 5221
		# cos(2*pi/5) & b0 = 1 - cos()
		c = 2*cos(2*pi/5 * 1:2)
		a1 = coeff[1]; a2 = coeff[2]
		b1 = coeff[3]; b2 = coeff[4]
		c1 = 1 - c[1]; c2 = 1 - c[2]
		x = c(
			roots(c(1, coeff[1]+coeff[2]*c[1], coeff[3]+coeff[4]*c[1], 1 - c[1])),
			roots(c(1, coeff[1]+coeff[2]*c[2], coeff[3]+coeff[4]*c[2], 1 - c[2])))
		coeff = c(1, 2*a1-a2, (a1^2-a2^2-a1*a2+2*b1-b2), (2*a1*b1 - 2*a2*b2 - a1*b2 - a2*b1 + 3),
			(b1^2-b2^2-b1*b2+3*a1+a2), (3*b1+b2), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == 5220) {
		# cos(2*pi/5) & b0 = 1
		c = 2*cos(2*pi/5 * 1:2)
		a1 = coeff[1]; a2 = coeff[2]
		b1 = coeff[3]; b2 = coeff[4]
		c1 = 1; c2 = 1
		x = c(
			roots(c(1, coeff[1]+coeff[2]*c[1], coeff[3]+coeff[4]*c[1], 1)),
			roots(c(1, coeff[1]+coeff[2]*c[2], coeff[3]+coeff[4]*c[2], 1)))
		coeff = c(1, 2*a1-a2, (a1^2-a2^2-a1*a2+2*b1-b2), (2*a1*b1 - 2*a2*b2 - a1*b2 - a2*b1 + 2),
			(b1^2-b2^2-b1*b2+2*a1-a2), (2*b1-b2), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == -5220) {
		# cos(2*pi/5) & b0 = -1
		# equivalent to 5220 with x = -x and a = -a;
		c = 2*cos(2*pi/5 * 1:2)
		a1 = coeff[1]; a2 = coeff[2]
		b1 = coeff[3]; b2 = coeff[4]
		c1 = -1; c2 = -1
		x = c(
			roots(c(1, coeff[1]+coeff[2]*c[1], coeff[3]+coeff[4]*c[1], -1)),
			roots(c(1, coeff[1]+coeff[2]*c[2], coeff[3]+coeff[4]*c[2], -1)))
		coeff = c(1, 2*a1-a2, (a1^2-a2^2-a1*a2+2*b1-b2), (2*a1*b1 - 2*a2*b2 - a1*b2 - a2*b1 - 2),
			(b1^2-b2^2-b1*b2-2*a1+a2), -(2*b1-b2), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == 21229) {
		# b0 = n +/- sqrt(n^2 - 1)
		n = coeff[1] # passed as 1st coefficient
		a1 = coeff[2]; a2 = coeff[3]
		b1 = coeff[4]; b2 = coeff[5]
		n2 = n^2 - 1; sq = sqrt(n2)
		c1 = n + sq; c2 = n - sq;
		x = c(
			roots(c(1, a1+a2*sq, b1+b2*sq, c1)),
			roots(c(1, a1-a2*sq, b1-b2*sq, c2)))
		coeff = c(1, 2*a1, (2*b1 + a1^2 - a2^2*n2), (2*n + 2*a1*b1 - 2*a2*b2*n2),
			(2*n*a1 - 2*a2*n2 + b1^2 - b2^2*n2), (2*n*b1 - 2*b2*n2), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == 24229) {
		# b0 = (n +/- sqrt(n^2 - 4)) / 2;
		# the remaining coeffs are liniar combinations of the radical:
		# b[j] => b[j, 0] + b[j, 1] * sqrt(n^2 - 4);
		n = coeff[1] # passed as 1st coefficient
		a1 = coeff[2]; a2 = coeff[3]
		b1 = coeff[4]; b2 = coeff[5]
		n2 = n^2 - 4; sq = sqrt(n2)
		c1 = (n + sq)/2; c2 = (n - sq)/2;
		x = c(
			roots(c(1, a1+a2*sq, b1+b2*sq, c1)),
			roots(c(1, a1-a2*sq, b1-b2*sq, c2)))
		coeff = c(1, 2*a1, (2*b1 + a1^2 - a2^2*n2), (n + 2*a1*b1 - 2*a2*b2*n2),
			(n*a1 - a2*n2 + b1^2 - b2^2*n2), (n*b1 - b2*n2), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else if(type == 20221) {
		# b0 = +1; # or - 1;
		# the remaining coeffs are liniar combinations of the radical:
		# b[j] => b[j, 0] + b[j, 1] * sqrt(n);
		n = coeff[1] # passed as 1st coefficient
		a1 = coeff[2]; a2 = coeff[3]
		b1 = coeff[4]; b2 = coeff[5]
		sq = sqrt(n)
		c1 = 1; c2 = 1;
		x = c(
			roots(c(1, a1+a2*sq, b1+b2*sq, c1)),
			roots(c(1, a1-a2*sq, b1-b2*sq, c2)))
		coeff = c(1, 2*a1, (a1^2 - n*a2^2 + 2*b1), 2*(a1*b1 - n*a2*b2+1),
			(2*a1 + b1^2 - n*b2^2), (2*b1), 1)
		err = sapply(x, function(x) sum(coeff*x^(6:0)) )
	} else {
		print("NOT yet implemented!")
	}
	err = round0(err)
	return(list(x=x, coeff=coeff, err=err))
}
### Solver
rootn = function(x, n=3) {
	if(Im(x) != 0 || Re(x) >= 0) return(x^(1/n))
	return( - (-x)^(1/n) )
}
solve.p3 = function(b.coeff, n=3) {
	# coeffs in DESCENDING order
	# TODO: correct implementation when c = m3;
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	m = m^(0:(n-1))
	c = - b.coeff[1]/n
	d = - b.coeff[2]/2
	det = d^2 - c^n
	det = if(Im(det) != 0 || Re(det) >= 0) sqrt(det) else complex(re=0, im=sqrt(-det))
	x = rootn(d + det, n) * m + rootn(d - det, n) / m
	return(x)
}
# generate special polynomials
# derived from cubic => sqrt => quadratic;
p6sq.gen = function(p3.coeff, mult=1, b0=-1, asSq=TRUE) {
	r = roots(p3.coeff)
	x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
	x = sapply(x1, function(r) roots(c(1, mult*r, b0)))
	p = poly.calc(x)
	for(i in 1:length(p)) p[[i]] = round0(p[[i]])
	if(asSq) {
		x = x[1,]^2 # the squares => P6
		p = polynomial(p[c(1,3,5,7,9,11,13)])
	}
	return(list(x=x, p=p))
}
### Other
toPoly = function(coeff, desc=TRUE, digits=5) {
	# TODO: test thoroughly!!!
	if(desc) { coeff = rev(coeff); }
	coeff = round(coeff, digits=digits)
	b0 = coeff[1]
	coeff = coeff[-1]
	isNotZero = (coeff != 0)
	countNotZero = sum(isNotZero)
	countPow = countNotZero; if(coeff[1] != 0) countPow = countPow - 1;
	isOne = (coeff[isNotZero] == 1 | coeff[isNotZero] == -1)
	coeff.txt = as.character(coeff[isNotZero])
	coeff.txt[isOne] = "" # TODO: if coeff = -1
	coeff.txt[isOne & coeff[isNotZero] < 0] = "- "
	op.txt = ifelse( isOne, "", "*")
	oppow.txt = rep("^", countPow)
	pow.txt = as.character(1:(length(coeff)));
	# x^1
	if(coeff[1] != 0) {
		pow.txt = c("", pow.txt[isNotZero][-1])
		oppow.txt = c("", oppow.txt)
	} else pow.txt = pow.txt[isNotZero]
	p.txt = paste(coeff.txt, op.txt, "x", oppow.txt, pow.txt, sep="", collapse=" + ")
	p.txt = paste(b0, p.txt, sep=" + ", collapse="")
	return(p.txt)
}


#######################

####################
### Convolved P3 ###
### Polynomials  ###
####################


#########################
### Roots of Unity M3 ###

### Type = (x^3, m, m, 1)
# Coeffs: x^3 + a*x^2 + b*x + c;
coeff = c(1, 1, -1)
solve.p6(coeff, type=110)

###
coeff = c(-1, 2, -1)
solve.p6(coeff, type=110)


########
### Type = (x^3, m, m, m)
coeff = c(1, 1, -1)
solve.p6(coeff, type=111)

###
coeff = c(1, 1, 1)
solve.p6(coeff, type=111)


###
coeff = c(2, 1, 1)
solve.p6(coeff, type=111)


########
### Type = (x^3, 1, 1, m)
coeff = c(1, 1, -1)
solve.p6(coeff, type=1)

###
coeff = c(2, 1, -1)
solve.p6(coeff, type=1)

###
coeff = c(1, 2, -1)
solve.p6(coeff, type=1)

###
coeff = c(1, -1, 2)
solve.p6(coeff, type=1)


###################
### cos(2*pi/5) ###

# Type 5229:
# - b1, b2: 4 coefficients of type b[j,0] + b[j,1]*cos();
# - b0 = 1 - cos() => prod(b0) = 1;

# the nice variants tend to be trivial;

###
coeff = c(1, 2, 2, -1)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 + 5*x + 10*x^2 + 8*x^3 + x^6

### 1 + ...*x^2 + ...*x^3 + x^6
b = -1
coeff = c(b, 2*b, b^2, -3*b^2)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 - 10*x^2 - 12*x^3 + x^6

###
b = 2
coeff = c(b, 2*b, b^2, -3*b^2)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 - 70*x^2 + 123*x^3 + x^6


### Nice, but trivial
coeff = c(0, 0, -1, -2)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 - 5*x - 5*x^2 + 3*x^3 + x^6
-5*(x + 1/2)^2 + (x^3 + 3/2)^2


### Nice, but trivial
coeff = c(0, 0, 1, 2)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 + 5*x - 5*x^2 + 3*x^3 + x^6
x = -x
1 - 5*x - 5*x^2 - 3*x^3 + x^6
-5*(x + 1/2)^2 + (x^3 - 3/2)^2


### Nice, but trivial
coeff = c(0, 0, 1/5, 2/5)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 + x - 1/5*x^2 + 3*x^3 + x^6
-5*(1/5*x - 1/2)^2 + (x^3 + 3/2)^2


### many more possible
coeff = c(0, 0, -3, -6)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 - 15*x - 45*x^2 + 3*x^3 + x^6


###
coeff = c(0, 0, -3/2, -3)
sol = solve.p6(coeff, type=5229)
sol


###
coeff = c(5, 10, 46, -33)
sol = solve.p6(coeff, type=5229)
sol

x = sol$x
1 + 105*x + 2570*x^2 + 828*x^3 + x^6


################
# Type 5220:
# - b1, b2: 4 coefficients of type b[j,0] + b[j,1]*cos();
# - b0 = 1 => prod(b0) = 1;

###
coeff = c(1, 2, 2, -1)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 + 5*x + 5*x^2 + 7*x^3 + x^6

### unfortunately b3 is complex
coeff = c(1i, 2i, -2, 1)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 - 5*x + 5*x^2 + (2-5i)*x^3 + x^6

### Nice: the -5*x^2 version is in the previous section;
coeff = c(-1, -2, 2, -1)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 + 5*x + 5*x^2 - 3*x^3 + x^6


###
coeff = c(0, -1, -1, -3)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 + x - 10*x^2 - 5*x^3 + x^5 + x^6


###
coeff = c(1, 1, -1, -3)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 + x - 10*x^2 + 10*x^3 + x^5 + x^6



# many trivial variants also generated
### Trivial
coeff = c(0, 0, -1, -2)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 - 5*x^2 + 2*x^3 + x^6

### Trivial
b = 3
coeff = c(0, 0, b, 2*b)
sol = solve.p6(coeff, type=5220)
sol

x = sol$x
1 - 5*b^2*x^2 + 2*x^3 + x^6


#####################
### sqrt(n^2 - 1) ###

# Type 21229:
# - b1, b2: 4 coefficients of type b[j,0] + b[j,1]*sqrt(n^2 - 1);
# - b0 = n +/- sqrt(n^2 - 1) => prod(b0) = 1;
# - function parameters: {n} is first coefficient;

### Test
coeff = c(3, 0, 1, 1, 2)
sol = solve.p6(coeff, type=21229)
sol

x = sol$x
1 - 26*x - 47*x^2 - 26*x^3 - 6*x^4 + x^6


### Test
coeff = c(3, 1, 1, 3, 1)
sol = solve.p6(coeff, type=21229)
sol

x = sol$x
1 + 2*x - 9*x^2 - 4*x^3 - x^4 + 2*x^5 + x^6


#####################
### sqrt(n^2 - 1) ###

# Type 24229:
# - b1, b2: 4 coefficients of type b[j,0] + b[j,1]*sqrt(n^2 - 4);
# - b0 = (n +/- sqrt(n^2 - 4)) / 2 => prod(b0) = 1;
# - function parameters: {n} is first coefficient;

### Test
coeff = c(3, -1, 1, 1, -1)
sol = solve.p6(coeff, type=24229)
sol

x = sol$x
1 + 8*x - 12*x^2 + 11*x^3 - 2*x^4 - 2*x^5 + x^6


### Test
coeff = c(3, -1, 2, -2, 1)
sol = solve.p6(coeff, type=24229)
sol

x = sol$x
1 - 11*x - 14*x^2 - 13*x^3 - 23*x^4 - 2*x^5 + x^6


### Test
coeff = c(23, -1, 1, 0, -1)
sol = solve.p6(coeff, type=24229)
sol

x = sol$x
1 + 525*x - 1073*x^2 + 1073*x^3 - 524*x^4 - 2*x^5 + x^6


###############
### sqrt(n) ###

# Type 20221:
# - b1, b2: 4 coefficients of type b[j,0] + b[j,1]*sqrt(n);
# - b0 = 1 => prod(b0) = 1;
# - function parameters: {n} is first coefficient;

# gsub("\\+ \\-","- ", toPoly(sol$coeff))

###
coeff = c(2, -1, 1, 1, -1)
sol = solve.p6(coeff, type=20221)
sol

x = sol$x
1 + 2*x - 3*x^2 + 4*x^3 + x^4 -2 *x^5 + x^6


###
coeff = c(2, -1, 1, 1, -2)
sol = solve.p6(coeff, type=20221)
sol

x = sol$x
1 + 2*x - 9*x^2 + 8*x^3 + x^4 - 2*x^5 + x^6


###
coeff = c(2, -1, 2, 1, -2)
sol = solve.p6(coeff, type=20221)
sol

x = sol$x
1 + 2*x - 9*x^2 + 16*x^3 - 5*x^4 - 2*x^5 + x^6

#############

#############
### Tests ###

### Large-Scale Tests

###
n.coeff = 2:20
sapply(n.coeff, function(x) toPoly(solve.p6(c(x, 0,0, 1,1), type=24229)$coeff))

###
n.coeff = 2:20
sapply(n.coeff, function(x) toPoly(solve.p6(c(x, 0,-1, 0,0), type=24229)$coeff))

sapply(n.coeff, function(x) toPoly(solve.p6(c(x, 0,-1, 1,0), type=24229)$coeff))

sapply(n.coeff, function(x) toPoly(solve.p6(c(x, 0,-1, 0,1), type=24229)$coeff))

### cos(2*pi/5) variants
n.coeff = (-10):10
sapply(n.coeff, function(x) toPoly(solve.p6(c(0,0, 1, x), type=5229)$coeff))


#########################


#########################
### P2 Decompositions ###
#########################

########################
### Radicals Order 3 ###
########################

# m = unity(3, all=TRUE)
# gsub("\\+ \\-","- ", toPoly(sol$coeff))

polyConv3.gen = function(K, b2, b1) {
	# b = descending order: s2*k^2 + s1*k + s0
	m = unity(3, all=TRUE)
	k = if(K < 0) - (-K)^(1/3) else K^(1/3)

	coeff.gen = function(id, b, pow=coeff.len(b)) sum(b*(k*m[id])^pow)
	coeff.len = function(b) {
		len = length(b) - 1
		return(len:0)
	}
	b2.x = sapply(1:3, coeff.gen, b=b2)
	b1.x = sapply(1:3, coeff.gen, b=b1)

	# TODO: (b0, b1, 1)
	coeff = lapply(1:3, function(id) c(1, b1.x[id], b2.x[id]))
	p = mult.p(mult.p(coeff[[1]], coeff[[2]]), coeff[[3]])
	p = round0(p)

	p.l = list(p=p, p.str = toPoly(Re(p)))
	return(p.l)
}

K = 2
polyConv3.gen(K, c(1, -1), c(1,-2,0))


K = 9
polyConv3.gen(K, c(1, -2), c(1,-2,0))


K = 28
polyConv3.gen(K, c(1, -3), c(1,-2,0))


sapply(1:10, function(id) polyConv3.gen( 2^3+1, c(1, -2), c(1,-2,id))$p.str)


sapply(1:10, function(id) polyConv3.gen( id^3+1, c(1, -id), c(1,-2,0))$p.str)

b = 2
sapply(1:10, function(id) polyConv3.gen( (id^3+1)/b^3, c(b, -id), c(1,-2,0))$p.str)

b = 3 # (8^3+1) / 27 == 19
sapply(1:10, function(id) polyConv3.gen( (id^3+1)/b^3, c(b, -id), c(1,-2,0))$p.str)



###########################
###########################

### Variants

### TODO:
# - clean; re-organize;

### n = 3
n = 3
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
# m = m^(0:(n-1))
### cos(2*pi/7) & cos(2*pi/5) entanglement
c3 = 2*cos(2*pi/7 * (1:3))
c2 = 2*cos(2*pi/5 * (1:2))
### other
r = 1 + c(1i,-1i)*sqrt(2)


### Example 1:
# coeffs only +/-1, no b5
x = c(solve.p3(c(-1, -m)) * m, solve.p3(c(-1, -m^2)) * m^2)
x
round0( 1 - x + x^2 + x^3 + x^4 + x^6 )
x = -x
round0( 1 + x + x^2 - x^3 + x^4 + x^6 ) # -x

###
x = c(roots(c(m,m^2,-m,1)), roots(c(m^2,m,-m^2,1)) )
x
round0( 1 + x - x^4 - x^5 + x^6 )
x = -1/x
round0( 1 + x - x^2 - x^5 + x^6 )


### Example 1x: trivial, for completion;
x = c(solve.p3(c(1, -m)), solve.p3(c(1, -m^2))) # actually {-1, -1, ...}
1 + x + x^2 + x^3 - x^4 + x^6
x = m^(1:2)
1 + x + x^2 + x^4 + x^5 + x^6 # (x^2 + x + 1)*(x^4 + 1)
1 + x - x^3 + x^5 + x^6 # (x^2 + x + 1)*(x^4 - x^2 + 1)
x = - m^(1:2)
1 - x + x^2 + x^3 + x^6 # (x^2 - x + 1)*(1 + x^3 + x^4)
x = c(1i, -1i)
1 + x^3 + x^5 + x^6
1 - x + x^5 + x^6
1 + x - x^5 + x^6
1 + x + x^3 + x^6
1 + x + x^2 + x^4 - x^5 + x^6
x = c(1, unity(5))
1 - x - x^5 + x^6
1 + x + x^5 + x^6 # -x
x = c(-1, unity(5))
-1 - x + x^5 + x^6
# partly
1 - x - x^3 + x^4 - x^5 + x^6 # (x-1)*(x^5 - x^3 - 1)


###############
### Example 2a: 2*x^4
x = c(solve.p3(c(-1, m)), solve.p3(c(-1, m^2)))
1 + x + x^2 - x^3 - 2*x^4 + x^6
1 - x + x^2 + x^3 - 2*x^4 + x^6 # -x
x = 1/x
1 - 2*x^2 - x^3 + x^4 + x^5 + x^6
x = c(roots(c(1,m^2,1i*sqrt(3),1)), roots(c(1,m,-1i*sqrt(3),1)))
1 + 2*x^2 - x^3 + x^4 - x^5 + x^6 # (x^2 + x - 1)*...

x = c(roots(c(1, c2[1], c2[2], -1)), roots(c(1, c2[2], c2[1], -1)))
1 + x + x^3 - 2*x^4 - x^5 + x^6
x = -1/x
1 + x - 2*x^2- x^3 - x^5 + x^6
x = unlist(lapply(c2, function(x) roots(c(x,-1,x-1,x)) ))
1 + x - 2*x^2 - x^3 - x^5 + x^6

### Example 2b: trivial: (x^2 + x + 1) * ...
x = c(solve.p3(c(m, m)) * m, solve.p3(c(m^2, m^2)) * m^2)
1 - x + x^2 - x^3 + 2*x^4 + x^6

# needs root shift or library(pracma)
x = c(roots(c(1,m^2,-1,1)) * m^2, roots(c(1,m,-1,1)) * m)
1 + x + 2*x^4 - x^5 + x^6
1 - x + 2*x^4 + x^5 + x^6 # -x

x = unlist(lapply(c2, function(x) roots(c(1,x,0,x)) ))
-1 - 2*x^2 - x^3 - x^4 - x^5 + x^6
x = unlist(lapply(c2, function(x) roots(c(x,-1,x,x-1)) ))
-1 + x - 2*x^2 + x^4 - x^5 + x^6


##############
### Example 3: 2*x^3
x = c(solve.p3(c(-m, 1)) * m, solve.p3(c(-m^2, 1)) * m^2)
1 + x + x^2 + 2*x^3 + x^4 + x^6
1 - x + x^2 - 2*x^3 + x^4 + x^6 # -x
x = 1/x
1 + x^2 + 2*x^3 + x^4 + x^5 + x^6
x = unlist(lapply(c2, function(x) roots(c(1,x,0,1)) ))
1 - x^2 + 2*x^3 - x^4 - x^5 + x^6

x = unlist(lapply(c2, function(x) roots(c(1,0,x,1)) ))
1 - x - x^2 + 2*x^3 - x^4 + x^6
x = -1/x
1 - x^2 - 2*x^3 - x^4 + x^5 + x^6
x = unlist(lapply(c2, function(x) roots(c(x,0,1,-x)) ))
1 - x - x^2 - 2*x^3 + x^4 + x^6

x = c(solve.p3(c(1, 1)) * m, solve.p3(c(1, 1)) * m^2)
1 - x + x^2 + 2*x^3 - x^4 + x^6
x = 1/x
1 - x^2 + 2*x^3 + x^4 - x^5 + x^6

# needs root shift or library(pracma)
x = c(roots(c(1,m^2,0,-1)) * m^2, roots(c(1,m,0,-1)) * m)
1 + x^2 - 2*x^3 + x^4 - x^5 + x^6

# trivial:
x = c(1i, -1i)
1 - x + x^2 - 2*x^3 + x^4 - x^5 + x^6
(x^2 + 1)*(x^2 + m*x + m^2)*(x^2 + m^2*x + m)



### Example 4a:
x = unlist(lapply(c2, function(x) roots(c(1,x,1,1)) ))
1 + 2*x + x^3 + x^4 - x^5 + x^6

x = c(solve.p3(c(-m^2, m^2)) * m, solve.p3(c(-m, m)) * m^2)
1 - 2*x + x^2 - x^3 + x^4 + x^6

### Example 4b:
x = c(solve.p3(c(-m, -m)) * m^2, solve.p3(c(-m^2, -m^2)) * m)
1 + 2*x + x^2 + x^3 + x^4 + x^6 # trivial from above

### Example 4c: (x^2 + x + 1) * ...
x = c(solve.p3(c(1, m)), solve.p3(c(1, m^2)))
1 + 2*x + x^2 - x^3 - x^4 + x^6

x = unlist(lapply(c2, function(x) roots(c(1,0,-x,x)) ))
-1 + 2*x - x^2 - x^3 + x^4 + x^6
x = unlist(lapply(c2, function(x) roots(c(1,-1,0,1/x)) ))
-1 - x^2 + x^3 + x^4 - 2*x^5 + x^6
x = unlist(lapply(c2, function(x) roots(c(1,1,x-1/x,1/x)) ))
-1 - x + 2*x^2 - x^3 - x^4 + 2*x^5 + x^6

##############
### Example 5: multiple coeffs of 2
x = c(solve.p3(c(m, 1)) * m^2, solve.p3(c(m^2, 1)) * m^2)
1 + 2*x + x^2 + 2*x^3 + 2*x^4 + x^6

x = c(solve.p3(c(-m^2, -1)), solve.p3(c(-m, -1)))
1 + 2*x + x^2 - 2*x^3 - 2*x^4 + x^6
x = -x
1 - 2*x + x^2 + 2*x^3 - 2*x^4 + x^6

# library(pracma)
x = c(roots(c(1,m^2,-m,-m)) , roots(c(1,m,-m^2,-m^2)) )
1 + 2*x + 2*x^2 + 2*x^3 + 2*x^4 - x^5 + x^6
x = -1/x
1 + x + 2*x^2 - 2*x^3 + 2*x^4 - 2*x^5 + x^6

x = c(roots(c(1i,1,0,1)) * 1i, roots(c(-1i,1,0,1)) * -1i)
1 - 2*x^2 - 2*x^3 + x^4 + 2*x^5 + x^6
x = c(roots(c(m,m^2,-1i*sqrt(3),1))*1, roots(c(m^2,m,1i*sqrt(3),1))*1)
1 + 2*x^2 + 2*x^3 - 2*x^4 - x^5 + x^6

x = c(roots(c(1,-m^2,1,m^2)) * m, roots(c(1,-m,1,m)) * m^2)
1 + 2*x + 2*x^2 - 2*x^5 + x^6
x = -1/x
1 + 2*x + 2*x^4 - 2*x^5 + x^6

### 2x2
x = c(roots(c(1,m^2,-m,-m)) * m^2, roots(c(1,m,-m^2,-m^2)) * m)
1 - x - x^2 + 2*x^3 + 2*x^4 - x^5 + x^6
x = 1/x
1 - x + 2*x^2 + 2*x^3 - x^4 - x^5 + x^6

x = c(roots(c(1,m,-1,m^2)) * m^2, roots(c(1,m^2,-1,m)) * m)
1 + x + 2*x^4 + 2*x^5 + x^6
x = 1/x
1 + 2*x + 2*x^2 + x^5 + x^6 # 1/x

x = c(roots(c(-m^2,-m^2,-m,1)) * 1, roots(c(-m,-m,-m^2,1)) * 1)
1 + x + 2*x^2 + 2*x^5 + x^6

x = unlist(lapply(c3, function(x) roots(c(1,1/x,-1/x+x)) ))
-1 - x + 2*x^3 - 2*x^5 + x^6

x = unlist(lapply(c2, function(x) roots(c(x,-1,-1,x)) ))
1 - x - 2*x^2 - 2*x^4 - x^5 + x^6 # trivial

x = c(roots(c(m,1+m^2,-m^2,1)) * m^2, roots(c(m^2,1+m,-m,1)) * m )
1 - 2*x - x^2 + x^3 + 2*x^4 + x^5 + x^6
x = 1/x
1 + x + 2*x^2 + x^3 - x^4 - 2*x^5 + x^6

x = c(roots(c(1,-m^2,m^2,-m^2)) * m, roots(c(1,-m,m,-m)) * m^2)
1 + x + 2*x^3 - 2*x^5 + x^6
x = c(roots(c(1,-m^2,m,-m^2)) * m^2, roots(c(1,-m,m^2,-m)) * m)
1 - 2*x + 2*x^3 + x^5 + x^6 # also 1/x

x = c(roots(c(m,-m,-m^2,1)) * m, roots(c(m^2,-m^2,-m,1)) * m^2 )
1 + x + 2*x^2 - 2*x^3 - x^4 + x^5 + x^6

# unusual
b = (1 + c(1,-1)*1i*sqrt(3))/2
x = c(roots(c(b[1], sqrt(3)*1i, 0,1)), roots(c(b[2], -sqrt(3)*1i, 0,1)))
1 + x^3 + 3*x^4 + 3*x^5 + x^6
x = c(roots(c(b[1], sqrt(3)*1i, 0,1)) * m, roots(c(b[2], -sqrt(3)*1i, 0,1))*m^2)
1 - 3*x^2 + x^3 + 3*x^4 - 3*x^5 + x^6
x = c(roots(c(b[1], sqrt(3)*1i, 0,1)) * m^2, roots(c(b[2], -sqrt(3)*1i, 0,1))*m)
1 + 3*x^2 + x^3 + 3*x^4 + x^6
x = c(roots(c(1,0,1i*sqrt(3),1))*m, roots(c(1,0,-1i*sqrt(3),1))*m^2)
1 + 3*x + 3*x^2 + 2*x^3 + 3*x^4 + x^6
#
x = c(roots(c(1i,1,0,1)) * m, roots(c(-1i,1,0,1)) * m^2)
1 - x^2 + x^4 + sqrt(3)*x^5 + x^6
x = c(roots(c(b[1], sqrt(3)*1i, 0,1)) * -1i, roots(c(b[2], -sqrt(3)*1i, 0,1))*1i)
1 + sqrt(3)*x^3 + 3*x^4 + sqrt(3)*x^5 + x^6
x = c(roots(c(1,1i*m^2,0,1)), roots(c(1,-1i*m,0,1)))
1 + sqrt(3)*x^2 + 2*x^3 + x^4 + sqrt(3)*x^5 + x^6
#
x = c(roots(c(1i/sqrt(3),1,0,1)) * m, roots(c(-1i/sqrt(3),1,0,1)) * m^2)
3 - 3*x^2 + 3*x^4 + 3*x^5 + x^6
x = c(roots(c(-1i/sqrt(3),1,m^2,1)) * m, roots(c(1i/sqrt(3),1,m,1)) * m^2)
3 - 3*x + 6*x^3 - 3*x^5 + x^6
x = c(roots(c(1/2-1i/2/sqrt(3),1,m^2,1)) * m, roots(c(1/2+1i/2/sqrt(3),1,m,1)) * m^2)
3 - 3*x + 9*x^3 - 3*x^5 + x^6
x = c(roots(c(1i,m^2,-1i,1)) * -1i, roots(c(-1i,m,1i,1)) * 1i)
1 + 2*x + 2*x^2 + 3*x^3 + 3*x^4 + x^5 + x^6
x = c(roots(c(1,1i*m^2*sqrt(3),0,1)), roots(c(1,-1i*m*sqrt(3),0,1)))
1 + 3*x^2 + 2*x^3 + 3*x^4 + 3*x^5 + x^6
x = c(roots(c(1,-1i*m^2*sqrt(3),1i*m/sqrt(3),1)), roots(c(1,1i*m*sqrt(3),-1i*m^2/sqrt(3),1)))
1 - x - 8/3*x^2 + 3*x^3 + 2*x^4 - 3*x^5 + x^6
x = c(roots(c(1,-1i*m^2*sqrt(3),1i/sqrt(3),1)), roots(c(1,1i*m*sqrt(3),-1i/sqrt(3),1)))
1 - 8/3*x^2 + 3*x^3 + 3*x^4 - 3*x^5 + x^6

### cos(2*pi/7) entanglement
c3 = 2*cos(2*pi/7 * (1:3))
x = c(roots(c(1,-1,1/c3[1])), roots(c(1,-1,1/c3[2])),roots(c(1,-1,1/c3[3])))
1 + x - 3*x^2 + 3*x^3 + x^4 - 3*x^5 + x^6

x = c(roots(c(1,-c3[1],c3[1])), roots(c(1,-c3[2],c3[2])),roots(c(1,-c3[3],c3[3])))
1 - 3*x + x^2 + 3*x^3 - 3*x^4 + x^5 + x^6

x = c(roots(c(1,-1,c3[1])), roots(c(1,-1,c3[2])),roots(c(1,-1,c3[3])))
1 + 2*x - 3*x^2 + x^3 + 2*x^4 - 3*x^5 + x^6
x = unlist(lapply(c3, function(x) roots(c(x,-1,x)) ))
1 + 2*x + 2*x^2 + 3*x^3 + 2*x^4 + 2*x^5 + x^6

x = unlist(lapply(c2, function(x) roots(c(x,0,x-1,1)) ))
-1 + 3*x - x^2 + x^3 + x^4 + x^6

x = c(roots(c(1,-1/c3[1],1/c3[1])), roots(c(1,-1/c3[2],1/c3[2])),roots(c(1,-1/c3[3],1/c3[3])))
1 - 3*x + 2*x^2 + x^3 - 3*x^4 + 2*x^5 + x^6

x = unlist(lapply(c3, function(x) roots(c(1,x,1/x+1)) ))
-1 + 3*x + 2*x^2 - 2*x^3 - x^4 - x^5 + x^6

x = unlist(lapply(c2, function(x) roots(c(1,-x,1,x)) ))
-1 - x + 3*x^2 + x^4 + x^5 + x^6

###
x = unlist(lapply(c3, function(x) roots(c(1,1/x,1/x+1)) ))
-1 - x - 5*x^3 - 2*x^5 + x^6


### Symmetric
x = sapply(1:3, function(id) roots(c(1, 3 - c3[id] - 2*c3[id]^2, 1)))
1 - 4*x^2 + 7*x^3 - 4*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 10+10*c3[id]-8*c3[id]^2-5*c3[id]^3, 1)))
1 - 18*x^2 - 7*x^3 - 18*x^4 + x^6
# r[2]
x = sapply(r, function(r) roots(c(1, -3+3*r, -3+3*r, -1)))
1 + 18*x^2 + 34*x^3 + 18*x^4 + x^6
x = sapply(r, function(r) roots(c(1, -3+3*r, -3+3*r, 1)))
1 + 18*x^2 + 38*x^3 + 18*x^4 + x^6
x = sapply(r, function(r) roots(c(1, -2+2*r, -2+2*r, 1)))
1 + 8*x^2 + 18*x^3 + 8*x^4 + x^6

### Asymmetric
x = sapply(1:3, function(id) roots(c(1, -10*c3[id]+2*c3[id]^2+5*c3[id]^3, 3-c3[id]-c3[id]^2)))
1 - x^2 - 14*x^3 - 16*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 2+c3[id]-c3[id]^2, c3[id])))
1 - 9*x^2 - 14*x^3 - 8*x^4 + x^6
# r[2]
x = sapply(r, function(r) roots(c(1, -1+r, -3+3*r, -1)))
1 + 2*x^2 + 10*x^3 + 18*x^4 + x^6

###
x = unlist(lapply(c3, function(x) roots(c(1,-1/x+2*x,x^2)) ))
1 + 7*x + 6*x^2 - 2*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 1 + 9*c3[id] - 2*c3[id]^2 - 4*c3[id]^3, 3-c3[id]-c3[id]^2)))
1 + 7*x + 6*x^2 - 2*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 6 + 5*c3[id] - 5*c3[id]^2-3*c3[id]^3, 1-c3[id]-c3[id]^2)))
1 + 7*x + 12*x^2 - 8*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 8 + 9*c3[id] - 7*c3[id]^2-5*c3[id]^3, 1-c3[id]^2)))
1 + 7*x + 13*x^2 - 9*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 1 + 1*c3[id] - 2*c3[id]^2-2*c3[id]^3, 1-c3[id]^2)))
1 + 7*x - 15*x^2 - 23*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1,-2 + 3*c3[id] + 1*c3[id]^2-1*c3[id]^3, -2-3*c3[id]-c3[id]^2)))
1 + 14*x + 40*x^2 - 15*x^4 + x^6


### Special
b = 4
x = sapply(r, function(r) roots(c(1, -1+r, (-3+b)+(1+r), -(1+r))))
6 - 4*b*x + (b^2-2*b-1)*x^2 + 2*b*x^4 + x^6
#
x = sapply(r, function(r) roots(c(1, -1+r, -2+(1+r), -(1+r))))
6 - 4*x - 2*x^2 + 2*x^4 + x^6
x = sapply(r, function(r) roots(c(1, -1+r, -4+(1+r), -(1+r))))
6 + 4*x + 2*x^2 - 2*x^4 + x^6
x = sapply(r, function(r) roots(c(1, -1+r,  0+(1+r), -(1+r))))
6 - 12*x + 2*x^2 + 6*x^4 + x^6
#
x = sapply(r, function(r) roots(c(1, -1+r, -3+2*r, -2*(1+r))))
24 - 8*x + x^2 + x^6
#
x = sapply(r, function(r) roots(c(1, 0, 1*(-1+r), -1+r)))
2 + 4*x + 2*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, -1+r)))
2 + 8*x + 8*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 0, -3+3*r, -1+r)))
2 + 12*x + 18*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, 2*(-1+r))))
8 + 16*x + 8*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, 3*(-1+r))))
18 + 24*x + 8*x^2 + x^6
#
x = sapply(1:3, function(id) roots(c(1, 2+c3[id]-c3[id]^2, 2-c3[id])))
7 + 14*x + 7*x^2 + x^6


b = 3
x = sapply(r, function(r) roots(c(1, -2+2*r, -5+r, -b)))
b^2 + 8*b*x + 18*x^2 + (8-2*b)*x^3 + x^6
x = sapply(r, function(r) roots(c(1, -2+2*r, -5+r, -1)))
1 + 8*x + 18*x^2 + 6*x^3 + x^6
x = sapply(r, function(r) roots(c(1, -2+2*r, -5+r, -4)))
16 + 32*x + 18*x^2 + x^6
x = sapply(r, function(r) roots(c(1, b*(-1+r), -b^2-1+r, -2*b)))
poly.calc(x)
4 + 4*x + 3*x^2 + x^6 # for b = 1
#
b = 2
sq = 3
x = sapply(sqrt(sq+0i)*c(1i,-1i), function(r) roots(c(1, b*r, -sq*b^2/2+r, -sq*b)))
poly.calc(x)

x = sapply(r, function(r) roots(c(1, 1-r, -2+r, -3+5*r)))
54 + 16*x - 17*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 1-r, -3+2*r, -1+5*r)))
66 + 32*x - 11*x^2 + x^6
x = sapply(r, function(r) roots(c(1, 1-r, -4+3*r, 1+5*r)))
86 + 48*x - x^2 + x^6
#
b = -1
x = sapply(r, function(r) roots(c(1, 1-r, -1 + b*(-1+r), 2*b + 5*(-1+r))))
poly.calc(x)
# e.g. 54 - 16*x - 17*x^2 + x^6 for b = -1


b3 = 2
b2 = -1-0i # alternative: b2 = 2
#
x = sapply(r, function(r) roots(c(1, 0, -sqrt(b2/2)+sqrt(b2/2)*r, b3/2)))
(b3/2)^2 + b2*x^2 + b3*x^3 + x^6
#
x = sapply(r, function(r) roots(c(1, 0, -sqrt(1/2)+sqrt(1/2)*r, 1)))
1 + x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -1+r, 1)))
1 + 2*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -sqrt(3/2)+sqrt(3/2)*r, 1)))
1 + 3*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -sqrt(2)+sqrt(2)*r, 1)))
1 + 4*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -sqrt(3)+sqrt(3)*r, 1)))
1 + 6*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, 1)))
1 + 8*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -sqrt(5)+sqrt(5)*r, 1)))
1 + 10*x^2 + 2*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -4+4*r, 1)))
1 + 32*x^2 + 2*x^3 + x^6
#
x = sapply(r, function(r) roots(c(1, 0, -1+r, 2)))
4 + 2*x^2 + 4*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, 2)))
4 + 8*x^2 + 4*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -4+4*r, 2)))
4 + 32*x^2 + 4*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -1+r, 4)))
16 + 2*x^2 + 8*x^3 + x^6
x = sapply(r, function(r) roots(c(1, 0, -2+2*r, 4)))
16 + 8*x^2 + 8*x^3 + x^6
#
x = sapply(r, function(r) roots(c(1, -2+2*r, -7+3*r, -1-2*r)))
17 + 18*x^2 + 18*x^3 + x^6
b = 1
x = sapply(r, function(r) roots(c(1, b*(-1+r), -b^2+(b+1)*(-1+r), -3+b*(1-r))))
poly.calc(x)
# free term in root: (2+2/b) + b*(1-r), but int only for b = {...}
x = sapply(r, function(r) roots(c(1, b*(-1+r), -b^2+(b+1)*(-1+r), -2-2/b+b*(1-r))))
poly.calc(x)

b = -1
x = sapply(r, function(r) roots(c(1, b*(-1+r), -b^2+(b+1)*(-1+r), 3+b*(-1+r))))
11 - 6*x + 5*x^2 + 6*x^3 + x^6 # for b = -1
17 - 16*x + 34*x^2 + 14*x^3 + x^6 # for b = -2

b = -3
x = sapply(r, function(r) roots(c(1, -2+2*r, -4+(2+b)*(-1+r), b+2*r)))
(b^2+4*b+12) + 2*(20+4*b+b^2)*x^2 + 10*(b+2)*x^3 + x^6
9 + 34*x^2 - 10*x^3 + x^6 # for b = -3
8 + 32*x^2 + x^6 # trivial for b = -2
9 + 34*x^2 + 10*x^3 + x^6 # for b = -1
x = sapply(r, function(r) roots(c(1, -2+2*r, -6+2*r, 2*r)))
12 + 40*x^2 + 20*x^3 + x^6
x = sapply(r, function(r) roots(c(1, -2+2*r, -7+3*r, 1+2*r)))
17 + 50*x^2 + 30*x^3 + x^6
b1 = 3
b = -3
x = sapply(r, function(r) roots(c(1, b1*(-1+r), -2*b1+(b+b1)*(-1+r), b+b1*r)))
poly.calc(x)
# eliminates x & x^5 terms
# TODO: more work

x = sapply(r, function(r) roots(c(1, -2+2*r, -4, 3-2*r)))
9 - 8*x + 2*x^3 + x^6


### x^6 + b1*x + b0
# see also below
# ("Entanglement with Roots of unity")
# e.g.: 3 - 3*x + x^6

x = sapply(r, function(r) roots(c(1, 2-2*r, -6+2*r, 5+3*r)))
82 - 40*x + x^6



################################

### Entanglement with Roots of unity

### nice:
coeff = c(1,-1, 1,2, -2,-1)
sol = solve.p6(coeff, type=222)
sol

x = sol$x
err = 3 - 3*x + x^6
round0(err)


###
coeff = c(1,-1, 1,1, 1,1)
sol = solve.p6(coeff, type=222)
sol

x = sol$x
err = 1 + 2*x + x^2 - 2*x^3 + x^4 + x^6
round0(err)


### Derivation
# Coeffs      | c1*m^2+c2*m      | b1*m^2+b2*m      | a1*m^2+a2*m      | 1
# c1*m+c2*m^2 | c1^2+c2^2-c1*c2  | b1*c1+b2*c2+     | a1*c1+a2*c2+     | c1*m+c2*m^2
#             |                  | b1*c2*m+b2*c1*m^2| a1*c2*m+a2*c1*m^2|
# b1*m+b2*m^2 | b1*c1+b2*c2+     | b1^2+b2^2-b1*b2  | a1*b1+a2*b2+     | b1*m+b2*m^2
#             | b1*c2*m^2+b2*c1*m|                  | a1*b2*m+a2*b1*m^2|
# a1*m+a2*m^2 | a1*c1+a2*c2+     | a1*b1+a2*b2+     | a1^2+a2^2-a1*a2  | a1*m+a2*m^2
#             | a1*c2*m^2+a2*c1*m| a1*b2*m^2+a2*b1*m|                  |
# 1           | c1*m^2+c2*m      | b1*m^2+b2*m      | a1*m^2+a2*m      | 1


x^6 - (a1+a2)*x^5 + (a1^2+a2^2-a1*a2-b1-b2)*x^4 + (2*a1*b1+2*a2*b2-a1*b2-a2*b1-c1-c2)*x^3 +
(b1^2+b2^2-b1*b2+2*a1*c1+2*a2*c2-a1*c2-a2*c1)*x^2 + (2*b1*c1+2*b2*c2-b1*c2-b2*c1)*x + c1^2+c2^2-c1*c2

coeff = c(1, - (a1+a2), (a1^2+a2^2-a1*a2-b1-b2), (2*a1*b1+2*a2*b2-a1*b2-a2*b1-c1-c2),
(b1^2+b2^2-b1*b2+2*a1*c1+2*a2*c2-a1*c2-a2*c1), (2*b1*c1+2*b2*c2-b1*c2-b2*c1), c1^2+c2^2-c1*c2)



########################

### Brute-Force Approach

n = 3
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
# m = m^(0:(n-1))
### cos(2*pi/7) & cos(2*pi/5) entanglement
c3 = 2*cos(2*pi/7 * (1:3))
c2 = 2*cos(2*pi/5 * (1:2))

# faster version (but still relatively slow)
for(s1 in (-1):10) {
for(s2 in (-10):10) {
	for(s3 in (-10):10) {
	for(s0 in (-2):12) {
		x.p = lapply(1:3, function(id) c(1, s0+s1*c3[id]+s2*c3[id]^2+s3*c3[id]^3, -2-1*c3[id]+2*c3[id]^2))
		p = round0(mult.p(x.p[[1]], x.p[[2]]))
		p = round0(mult.p(p, x.p[[3]]))
		sum0 = sum(p == 0)
		# cat(sum0); cat(", ")
		if(sum0 > 1 && any(p[c(2,4,6)] != 0) && max(abs(p)) <= 200) {
			print(p)
			if(sum0 > 2) cat("==> !! ")
			print(c(s0, s1, s2, s3))
		}
	}
}
}
}

# very slow version
for(s1 in (-10):10) {
for(s2 in (-10):10) {
	for(s3 in (-10):10) {
	for(s0 in (-2):10) {
		x = sapply(1:3, function(id) roots(c(1, s0+s1*c3[id]+s2*c3[id]^2+s3*c3[id]^3, -2-3*c3[id]-1*c3[id]^2)))
		p = round(poly.calc(x))
		p1 = as.vector(p)
		sum0 = sum(p1 == 0)
		# cat(sum0); cat(", ")
		if(sum0 > 1 && any(p1[c(2,4,6)] != 0) && max(abs(p1)) <= 200) {
			print(p)
			if(sum0 > 2) cat("==> !! ")
			print(c(s0, s1, s2, s3))
		}
	}
}
}
}


# c(1,2,-2), c(1,2,-1), c(1,2,-2),
# c(1,3,1), c(3,-1,-1), c(4,6,-7), c(1,3,2), c(0,2,1), c(1,0,-1)
# c(-2,-3,-1), c(-2,-1,2)
for(s0 in (-10):10) {
for(s1 in (-10):10) {
	for(s2 in (-10):10) {
		pr = round(prod(s0 + s1*c3 + s2*c3^2))
		# cat(sum0); cat(", ")
		if(pr == 1) {
			print(c(s0, s1, s2))
		}
	}
}
}

#######

r = 1 + c(1i,-1i)*sqrt(2)
conj.roots = function(r) {
	r = rev(r)
	c(1, s20+s21*r, s10 + s11*r, s00+s01*r)
}

for(s11 in (-6):6) {
for(s10 in (-7):7) {
for(s01 in (-7):7) {
for(s00 in (-7):7) {
	for(s21 in (-3):5) {
	for(s20 in (-3):5) {
		x.p = lapply(r, conj.roots)
		p = round0(mult.p(x.p[[1]], x.p[[2]]))
		if(p[7] == 0) next
		sum0 = sum(p == 0)
		# cat(sum0); cat(", ")
		if(sum0 > 2 && any(p[c(2,4,6)] != 0) && max(abs(p)) <= 200) {
			print(p)
			if(sum0 > 2) cat("==> !! ")
			print(c(s20, s21, s10, s11, s00, s01))
		}
		p.rm = p[-7]
		if((p[7] == 1 && sum0 > 1) || all(p.rm == 1 | p.rm == -1 | p.rm == 0)) {
			cat("\n\n !!! ==>\n")
			print(p)
			print(c(s20, s21, s10, s11, s00, s01))
		}
	}
}}}}
}


###
p.c = t(sapply(1:3, function(id) c(c3[id]^2-c3[id], c3[id]^3-1, 1)))
mult.p(mult.p(p.c[1,], p.c[2,]), p.c[3,])
x = roots(rev(p.c[1,]))
1 - 2*x^2 + 6*x^4 - 7*x^5 + x^6

####################

p6sq.gen = function(p3.coeff, mult=1, b0=-1, asSq=TRUE) {
	r = roots(p3.coeff)
	x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
	x = sapply(x1, function(r) roots(c(1, mult*r, b0)))
	if(asSq) {
		x = x[1,]^2
	}
	p = poly.calc(x)
	for(i in 1:length(p)) p[[i]] = round0(p[[i]])
	return(list(x=x, p=p))
}

### Examples

sapply(-6:6, function(b) print(p6sq.gen(c(1, 5, 5, b - 0))$p))
# the values of the roots must be extracted explicitly from p6sq.gen();
1 - x - 6*x^3 - x^5 + x^6 
1 - x - 5*x^3 - x^5 + x^6 
1 - x - 4*x^3 - x^5 + x^6 
1 - x - 3*x^3 - x^5 + x^6 
1 - x - 2*x^3 - x^5 + x^6 
1 - x - x^3 - x^5 + x^6 
1 - x - x^5 + x^6 
1 - x + x^3 - x^5 + x^6 
1 - x + 2*x^3 - x^5 + x^6 
1 - x + 3*x^3 - x^5 + x^6 
1 - x + 4*x^3 - x^5 + x^6 
1 - x + 5*x^3 - x^5 + x^6 
1 - x + 6*x^3 - x^5 + x^6

sapply(-13:3, function(b) print(p6sq.gen(c(1,0, -15,b))$p))
sapply(-13:3, function(b) print(p6sq.gen(c(1,1, -11,b))$p))
sapply(-5:5, function(b) print(p6sq.gen(c(1,5, b, -2))$p))
sapply(-6:6, function(b) print(p6sq.gen(c(1,5, b, -3))$p))
sapply(-6:6, function(b) print(p6sq.gen(c(1,-2, -22, b - 8))$p))
sapply(-6:6, function(b) print(p6sq.gen(c(1,-3, -26, b - 14))$p))
sapply(-6:6, function(b) print(p6sq.gen(c(1,-3, -28, b - 16))$p))

### 1 - x - x^2 - x^4 - x^5 + x^6
p = p6sq.gen(c(1,5, 4,-2))
p
x = p$x[1,]^2
1 - x - x^2 - x^4 - x^5 + x^6


### 1 + x - x^2 - x^4 + x^5 + x^6 [same as previous]
# the initial workout of the solution;
x = solve(polynomial(c(1,1,-1,0,-1,1,1)))
1 + x - x^2 - x^4 + x^5 + x^6

poly.calc(x-1/x)
x.r = x - 1/x
-4 - 4*x.r^2 + 3*x.r^4 + x.r^6
poly.calc((x-1/x)[c(1,3,5)]^2 + 1)
2 - 7*x + x^3
# Solution starts from: 2 - 7*x + x^3;


# Analysis [old]
p = poly.calc(x[c(3:6)])
p
# p has the coefficients as roots of the following polynomials:
x.r = p[[2]]
x.r^3 - 2*x.r^2 - 3*x.r + 2
x.r = p[[3]]
x.r^3 - 2*x.r^2 - 6*x.r + 8

p = poly.calc(x[c(1:2)])
p
# p has the coefficient as roots of the following polynomial:
x.r = p[[2]]
x.r
x.r^3 - x.r^2 - 4*x.r + 2
2 - 4*x.r - x.r^2 + x.r^3


##########


### 1 - 3*x - 3*x^5 + x^6
r = roots(c(1,3,-3,-4))
x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
x = sapply(x1, function(r) roots(c(1,r,-1)))
poly.calc(x)
1 - 3*x^2 - 3*x^10 + x^12
x = x[1,]^2
1 - 3*x - 3*x^5 + x^6


### 1 - 3*x + 2*x^3 - 3*x^5 + x^6
r = roots(c(1,3,-3,-2))
x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
x = sapply(x1, function(r) roots(c(1,r,-1)))
poly.calc(x)
1 - 3*x^2 + 2*x^6 - 3*x^10 + x^12
x = x[1,]^2
1 - 3*x + 2*x^3 - 3*x^5 + x^6


### 1 - 6*x^2 + 8*x^3 - 6*x^4 + x^6
r = roots(c(1,6,+3,-2))
x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
x = sapply(x1, function(r) roots(c(1,r,-1)))
poly.calc(x)
1 - 6*x^4 + 8*x^6 - 6*x^8 + x^12
x = x[1,]^2
1 - 6*x^2 + 8*x^3 - 6*x^4 + x^6


### 1 - 10*x^2 + 16*x^3 - 10*x^4 + x^6
r = roots(c(1,6,-1,-2))
x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
x = sapply(x1, function(r) roots(c(1,r,-1)))
poly.calc(x)
1 - 10*x^4 + 16*x^6 - 10*x^8 + x^12
x = x[1,]^2
1 - 10*x^2 + 16*x^3 - 10*x^4 + x^6


### 1 - 3*x + 4*x^2 - 6*x^3 + 4*x^4 - 3*x^5 + x^6
r = roots(c(1,3,1,-2))
x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
x = sapply(x1, function(r) roots(c(1,r,-1)))
poly.calc(x)
1 - 3*x^2 + 4*x^4 - 6*x^6 + 4*x^8 - 3*x^10 + x^12
x = x[1,]^2
1 - 3*x + 4*x^2 - 6*x^3 + 4*x^4 - 3*x^5 + x^6


###
sapply(-6:6, function(b) print(p6sq.gen(c(1, 5, 6, b))$p))

sapply(-6:6, function(b) print(p6sq.gen(c(1,5, b,-2))$p))

sapply(-10:10, function(b) print(p6sq.gen(c(1, 5, 8, b))$p))


# the full series:
sapply(-6:6, function(b) print(p6sq.gen(c(1,5, 5, b))$p))
# 1 - x - b*x^3 - x^5 + x^6 # b from -6 to +6 (fully generalized);

### 1 - x - x^3 - x^5 + x^6
p = p6sq.gen(c(1,5, 5, -1))
x = p$x[1,]^2
1 - x - x^3 - x^5 + x^6

### 1 - x + x^3 - x^5 + x^6
# factorization: (x^2-x+1), but still part of a series;
p = p6sq.gen(c(1,5, 5, 1))
x = p$x[1,]^2
1 - x + x^3 - x^5 + x^6


### 1 - x - 4*x^3 - x^5 + x^6
p = p6sq.gen(c(1,5, 5,-4))
p
x = p$x[1,]^2
1 - x - 4*x^3 - x^5 + x^6

### trivial, but still interesting

### 1 - x - 3*x^3 - x^5 + x^6
# factorization: (x^2+x+1)*...
p = p6sq.gen(c(1,5, 5,-3))
p
x = p$x[1,]^2
1 - x - 3*x^3 - x^5 + x^6

### 1 - x - 2*x^3 - x^5 + x^6
# roots: include +i, -i
p = p6sq.gen(c(1,5, 5,-2))
p
x = p$x[1,]^2
1 - x - 2*x^3 - x^5 + x^6


###
p6sq.gen(c(1,5,-1,-2))
1 - x^2 - 6*x^4 + 10*x^6 - 6*x^8 - x^10 + x^12

