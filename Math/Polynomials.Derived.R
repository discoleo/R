########################
###
### Leonard Mada
### [the one and only]
###
### Derived Polynomials
### v.0.4c-coeff-x1

### Note:
# This is the 1st part towards:
# - introducing various derived polynomials,
#   e.g. based on polynomials of Class 1;
# - brief introduction of polynomials of Class 1;
# - introducing interesting properties of polynomial roots;
# - includes a different approach to polynomials.

### History
# v.0.4b - v.0.4c:
# - [started] parametric P[5] derived from:
#   x^5 - x + K = 0;
# v.0.4.a:
# - added various derived polynomials starting from:
#   x^n = x^j + K => roots of: x^j = r^n - K;
#   Example: x^4 + x^3 + 1 = 0
#   => 1 + 2*x^3 - x^4 + x^6 - x^7 + x^8 = 0;
# v.0.3.x:
# - moved P6 polynomials to:
#   Polynomials.Derived.P6.R
# v.0.3.c - v.0.3.f:
# - more awesome P6 polynomials with complete solutions;
# - moved now to separate file (see v.0.3.x);
#   1 + x - x^4 - x^5 + x^6 = 0;
#   1 - x + x^2 + x^3 + x^4 + x^6 = 0;
#   1 + x + x^2 - x^3 - 2*x^4 + x^6 = 0;
#   1 + 2*x + x^2 - x^3 - x^4 + x^6 = 0;
#   & many more;
# v.0.3.b:
# - a nice polynomial: 1 - x + x^2 + x^5 + x^6 + x^10
#   [solution based on roots of: x^5 - x - 1 = 0]
# v.0.3a: added a Class 3 Polynomial example (roots based on cos());
# v.0.2b: a slight generalization;
# v.0.2a: new technique to construct interesting polynomials;
# v.0.1: first drafts;

### Theory:
# - let P[n] be a polynomial of order n with integer(/rational) coefficients;
# - let r be the n roots of this polynomial;
# - let f be a polynomial function with integer coefficients;
# - then f(r) are the roots of a polynomial of order n with integer(/rational) coefficients;


library(polynom)
library(pracma) # needed if the coefficients are complex

# needed to get the roots of the base polynomial;
# the derived polynomial is currently derived "empirically";
# [it is relatively easy to derive parametrically as well]

######################

### helper functions

# - see file: Polynomials.Helper.R;
# - TODO: clean section;

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
# safe sqrt
sqrt.c = function(x) {
	x = as.complex(x)
	x.sqrt = sqrt(x)
	return(x.sqrt)
}
# safe n-th root
nroot.c = function(x, n=3) {
	if(Im(x) == 0) {
		x = Re(x)
		if(x >= 0) {
			x = complex(re = x^(1/n), im=0)
		} else if(n %% 2 == 0) {
			x = complex(re = 0, im=(-x)^(1/n))
		} else {
			x = - (-x)^(1/n)
		}
	} else {
		x = x^(1/n)
	}
	return(x)
}
### Solve P3
solve.P3Base <- function(c, d, n=3, all=TRUE) {
	det = (d^2 - c^n)
	det = sqrt.c(det)
	
	p = nroot.c(d + det, n)
	q = if(Re(d) >= Re(det)) { nroot.c(d - det, n); } else { - nroot.c(-d + det, n); }
	if(all) {
		m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
		id = 0:(n-1)
		m = m^id
		return(p*m + q/m)
	}
	return(p + q)
}
shift.poly = function(b) {
	# assumes b = c(b[n-1], ..., b0)
	len = length(b)
	shift = - b[1] / len
	# TODO: len > 3
	b.shifted = b
	b = c(1, b)
	for(i.base in 1:len) {
		for(id in i.base:len) {
			pow = id + 1 - i.base
			b.shifted[id] = b.shifted[id] + b[i.base] * choose(len + 1 - i.base, pow) * shift^pow
		}
		# print(b.shifted)
	}
	b.shifted[1] = 0
	return(list("b"=b.shifted, "shift"=shift))
}
solve.P3 = function(b) {
	# Coefficients in DESC powers
	len = length(b)
	if(len == 3 && b[1] != 0) {
		b.shifted = shift.poly(b)
		b = b
	} else if(len == 4 && b[2] != 0) {
		b.shifted = shift.poly(b[-1])
	} else {
		if(len < 3) { b = c(0, b); }
		else if(len > 3) { b = b[(len-2):len]; }
		b.shifted = list("b"=c(0, b[c(len-1, len)]), "shift"=0)
	}
	b = b.shifted$b
	print(b.shifted$b)
	return (solve.P3Base(-b[2]/3, -b[3]/2) + b.shifted$shift)
}

##############

#############
### Examples:

### Base polynomial:
# x^5 - x - 1
p = polynomial(c(-1,-1,0,0,0,1))
p

x0 = solve(p)
x0

### Derived polynomials

# lets create various derived polynomials of order 5 with integer coefficients;


###
x = x0^3 + x0^2 + x0
x

poly.calc(x)
err = -4 - 14*x - 21*x^2 - 11*x^3 + x^5
round0(err)

# Test
# x0 = 2 # other values to Test composition & decomposition
x = x0^3 + x0^2 + x0
err = -4 - 14*x - 21*x^2 - 11*x^3 + x^5
round0(err)
x = x0
# p.calc() o (x^3 + x^2 + x) =
err = - 4 - 14*x - 35*x^2 - 67*x^3 - 96*x^4 - 107*x^5 - 93*x^6 - 51*x^7 - 3*x^8 + 34*x^9 + 51*x^10 + 45*x^11 + 30*x^12 + 15*x^13 + 5*x^14 + x^15
round0(err)
# ==
err = (4 + 10*x + 25*x^2 + 42*x^3 + 54*x^4 + 57*x^5 + 46*x^6 + 30*x^7 + 15*x^8 + 5*x^9 + x^10) *
(x^5 - x - 1)
round0(err)
# == 0
# x = x0^3 + x0^2 + x0 is indeed a root of the above polynomial!


##############

###
x = x0^4 - x0^2 + x0^3 - x0
x

poly.calc(x)
err = -1 + 3*x - 16*x^2 + 18*x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^3 - x0
x

poly.calc(x)
err = -1 + 5*x - 8*x^2 + 4*x^3 + x^5
round0(err)

###
x = x0^3 + x0
x

poly.calc(x)
err = -1 - 5*x - 8*x^2 - 4*x^3 + x^5
round0(err)

###
x = x0^4 + x0^3 + x0^2 + x0
x

poly.calc(x)
err = -1 - 5*x - 10*x^2 - 10*x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^4 - x0^3 + x0^2 - x0
x

poly.calc(x)
err = -1 + 5*x - 10*x^2 + 10*x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^4 + x0
x

poly.calc(x)
err = -4 - 2*x + 7*x^2 + x^3 - 4*x^4 + x^5
round0(err)

###
x = x0^4 - x0
x

poly.calc(x)
err = -4 + 12*x - 15*x^2 + 11*x^3 - 4*x^4 + x^5
round0(err)

###
# works because free term of the initial polynomial is +/- 1!
x = x0 + 1/x0
x

poly.calc(x)
err = -1 + 4*x - 4*x^2 - 5*x^3 + x^4 + x^5
round0(err)

###
# works because free term of the initial polynomial is +/- 1!
x = 1/x0 - x0
x

poly.calc(x)
err = 1 + 4*x + 4*x^2 + 5*x^3 + x^4 + x^5
round0(err)


### x^5 - 3*x^2 - x - 1
x = 1/(x0^2 - x0 + 1 - 2/(x0^2+x0+1))
# (x0^2 + x0 + 1)/(x0^4 + x0^2 - 1)
# (x0^2 - x0 - 1)/(x0^4 - x0^2 - 1)
# == x0^3
x

poly.calc(x)
err = x^5 - 3*x^2 - x - 1
round0(err)
#
x = x0
x^15 - 3*x^6 - x^3 - 1
# *nice* polynomial:
m3 = complex(re=cos(2*pi/3), im=sin(2*pi/3))
x = c(x0*m3, x0*m3^2)
err = 1 - x + x^2 + x^5 + x^6 + x^10
round0(err)

x = c(1/x0 * m3, 1/x0*m3^2, x0*m3, x0*m3^2)
poly.calc(x)
1 - x + x^2 + x^4 + x^5 + x^6 + x^7 + x^8 - x^9 + 6*x^10 +
- x^11 + x^12 + x^13 + x^14 + x^15 + x^16 + x^18 - x^19 + x^20


##########################

##################
### Parametric ###
##################

### Derivation:
p = roots.derived(5)
p = mult.lpm(p)
p = sort.pm(p, c(4,3), "x")
p = p[, c("x", paste0("r", 5:1), paste0("s", 4:1), "coeff")]
rownames(p) = seq(nrow(p))

# helper functions:
# - moved to Polynomials.Helper.R & Polynomials.Helper.EP.R;


K = -1
# x^5 - x + K
x0 = roots(c(1,0,0,0,-1, K))

### Examples:
s4=1; s3=-5; s2=0; s1=3;
r = sapply(seq(5), function(id) sum(x0[id]^seq(4) * c(s1,s2,s3,s4)));
round0.p(poly.calc(r))
#
eval.pm(p[p$x == 1,], c(1, x0, s4,s3,s2,s1))

x^5 - 4*s4*x^4 + (6*s4^2 + 5*s1*s4*K + 5*K*s3*s2 - 4*s3*s1 - 2*s2^2)*x^3 +
	- (5*s4*s3^2*K^2 + 5*s4^2*s2*K^2 - 5*s2^2*s1*K - 5*s3*s1^2*K - 3*s3^3*K + 2*s4*s3*s2*K +
		+ 11*s4^2*s1*K + 4*s2*s1^2 + 4*s3^2*s2 - 4*s4*s2^2 - 8*s4*s3*s1 + 4*s4^3)*x^2 +
	+ (5*s4^3*s3*K^3 + 5*s3^2*s2^2*K^2 - 5*s4*s2^3*K^2 - 5*s3^3*s1*K^2 - 5*s4*s3*s2*s1*K^2 + 5*s4^2*s1^2*K^2 +
		- s4^2*s3^2*K^2 + 6*s4^3*s2*K^2 + 5*s2*s1^3*K - s3*s2^3*K + 7*s3^2*s2*s1*K - 3*s4*s2^2*s1*K +
		- 13*s4*s3*s1^2*K + s4*s3^3*K - 3*s4^2*s3*s2*K + 7*s4^3*s1*K - s1^4 + s2^4 - s3^4 + s4^4 +
		- 4*s3*s2^2*s1 + 2*s3^2*s1^2 + 4*s4*s2*s1^2 + 4*s4*s3^2*s2 - 2*s4^2*s2^2 - 4*s4^2*s3*s1)*x +
	- s4^5*K^4 + (s3^5 - s4^4*s3 - 5*s4*s3^3*s2 + 5*s4^2*s3*s2^2 + 5*s4^2*s3^2*s1 - 5*s4^3*s2*s1)*K^3 +
		- s2^5*K^2 + 5*s3*s2^3*s1*K^2 - 5*s3^2*s2*s1^2*K^2 - 5*s4*s2^2*s1^2*K^2 + 5*s4*s3*s1^3*K^2 +
		+ s3^4*s2*K^2 - 4*s4*s3^2*s2^2*K^2 + 2*s4^2*s2^3*K^2 - s4*s3^3*s1*K^2 + 7*s4^2*s3*s2*s1*K^2 +
		- 3*s4^3*s1^2*K^2 - s4^4*s2*K^2 +
		+ s1^5*K - s2^4*s1*K + 4*s3*s2^2*s1^2*K - 2*s3^2*s1^3*K - 4*s4*s2*s1^3*K + s3^4*s1*K +
		- 4*s4*s3^2*s2*s1*K + 2*s4^2*s2^2*s1*K + 4*s4^2*s3*s1^2*K - s4^4*s1*K;


### Derivation:
# p[p$x == 2,]
# sort.rpm(p[p$x == 2, ])
# unique.rpm(p[p$x == 2, ])
replace.rpm(p[p$x == 1, ])
print.p(replace.pm(replace.rpm(p[p$x == 0, ])[-1], data.frame(K=1, coeff=-1), "E5"), "K")


# x^5 + (s4*S4 + s3*S3 + s2*S2 + s1*S1) +
#	+ (s4^2*E2_44 + s3*s4*E2_43 + s2*s4*E2_42 + s1*s4*E2_41 +
#		+ s3^2*E2_33 + s3*s2*E2_32 + s3*s1*E2_31 +)
#		+ s2^2*E2_22 + s2*s1*E2_21)*x^3 +
#	- (s4*3*E3_444 + s4^2*s3*E3_443 + s4*s3^2*E3_433 +
#		+ s4^2*s2*E3_442 + s4*s3*s2*E3_432 + s4*s2^2*E4_422 +
#		+ s4^2*s1*E3_441 + s4*s3*s1*E3_431 + s4*s2*s1*E3_421 + s4*s1^2*E3_411 +
#		+ s3^3*E3_333 + s3^2*s2*E3_332 + s3*s2^2*E3_322 +
#		+ s3^2*s1*E3_331 + s3*s2*s1*E3_321 + s3*s1^2*E3_311 +
#		+ s2^3*E3_222 + s2^2*s1*E3_221 + s2*s1^2*E3_211 + s1^3*E3_111)*x^2 + ...


### Sn = sum(r^n)
# Epoly.gen(n, 5, 1)
S4 = -4; # S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4
S3 = 0; # (S^3 - 3*E2*S + 3*E3 = 0;)
S2 = 0; S1 = 0;
### E2_nm
# Epoly.distinct(c(n,m), 5)
# print.p(Epoly.distinct(c(4,4), 5), "S")
# - Note: S = 0 (and omitted); E2 = 0; E3 = 0; E4 = -1;
E2_44 = 6; # E2^4 + 4*E2*E3^2 - 4*E2^2*E4 + 6*E4^2 - 4*E3*E5
E2_43 = 0; E2_42 = 0;
E2_41 = 5*K; # (5*E2*E3 - 5*E5 = 5*K;)
E2_33 = 0; # E2^3 + 3*E3^2 - 3*E2*E4
E2_32 = 5*K; # - E2*E3 - 5*E5
E2_31 = -4; # - 2*E2^2 + 4*E4
E2_22 = -2; # E2^2 + 2*E4
E2_21 = 0;


##########################
##########################

### From: x^5 - x^2 + K

### Derivation:
p = roots.derived(5)
p = mult.lpm(p)
p = sort.pm(p, c(4,3), "x")
p = p[, c("x", paste0("r", 5:1), paste0("s", 4:1), "coeff")]
rownames(p) = seq(nrow(p))

K = -1
# x^5 - x^2 + K
x0 = roots(c(1,0,0,-1, 0, K))

### Examples:
s4=1; s3=-5; s2=0; s1=3;
r = sapply(seq(5), function(id) sum(x0[id]^seq(4) * c(s1,s2,s3,s4)));
round0.p(poly.calc(r))
#
eval.pm(p[p$x == 1,], c(1, x0, s4,s3,s2,s1))

x^5 - 3*s3*x^4 +
	+ (4*s4^2*K + 5*s3*s2*K + 5*s4*s1*K + 3*s3^2 - 3*s4*s2 - 3*s2*s1)*x^3 +
	+ (- 5*s4*s3^2*K^2 - 5*s4^2*s2*K^2 - s4^2*s3*K - 7*s3^2*s2*K + 8*s4*s2^2*K + s4*s3*s1*K + 5*s2^2*s1*K +
		+ 5*s3*s1^2*K - s4^3 - s3^3 - s2^3 - s1^3 + 3*s4*s3*s2 - 3*s4^2*s1 + 3*s3*s2*s1 - 3*s4*s1^2)*x^2 +
	+ (5*s4^3*s3*K^3 + 2*s4^4*K^2 + 2*s4*s3^3*K^2 - 5*s4*s2^3*K^2 + 7*s4^3*s1*K^2 - 5*s3^3*s1*K^2 +
		- 4*s4^2*s3*s2*K^2 + 5*s3^2*s2^2*K^2 - 5*s4*s3*s2*s1*K^2 + 5*s4^2*s1^2*K^2 + 2*s2^4*K +
		+ 2*s4^3*s2*K + 2*s3^3*s2*K + 5*s2*s1^3*K - 6*s4*s3*s2^2*K + 6*s4^2*s2*s1*K - 6*s3*s2^2*s1*K +
		- 3*s3^2*s1^2*K + 9*s4*s2*s1^2*K)*x +
	+ (- s4^5*K^4 + s3^5*K^3 - 2*s4^4*s2*K^3 + s4^3*s3^2*K^3 - 5*s4*s3^3*s2*K^3 + 5*s4^2*s3*s2^2*K^3 +
		+ 5*s4^2*s3^2*s1*K^3 - 5*s4^3*s2*s1*K^3 - s2^5*K^2 + 2*s3^4*s1*K^2 - s4^3*s2^2*K^2 - s3^3*s2^2*K^2 +
		+ 3*s4*s3*s2^3*K^2 + 2*s4^3*s3*s1*K^2 + 5*s3*s2^3*s1*K^2 + 5*s4*s3*s1^3*K^2 - 6*s4*s3^2*s2*s1*K^2 +
		- 3*s4^2*s2^2*s1*K^2 + 7*s4^2*s3*s1^2*K^2 - 5*s3^2*s2*s1^2*K^2 - 5*s4*s2^2*s1^2*K^2 + s1^5*K +
		+ 3*s4*s1^4*K + s4^3*s1^2*K + s3^3*s1^2*K + s2^3*s1^2*K + 3*s4^2*s1^3*K - 3*s3*s2*s1^3*K - 3*s4*s3*s2*s1^2*K)


### Derivation:
rpl = list(
	"E5" = data.frame(K=1, coeff=-1),
	"E3" = data.frame(coeff=1)
);

pK = coef.rpm(p, flt=c("S", "E2", "E4"), rpl=rpl)
sapply(5:0, function(pow) cat("(", print.p(pK[pK$x == pow, - match("x", names(pK))], "K"), ")*x^", pow, " +\n", sep=""))



##########################

##########################
### Decomposition of   ###
### derived polynomial ###

### Base polynomial:
# x^5 - x - 1
p = polynomial(c(-1,-1,0,0,0,1))
p

x0 = solve(p)
x0

###
x = x0^2 + x0
x

poly.calc(x)
err = -1 - 5*x - 9*x^2 - 2*x^3 + x^5
round0(err)
# unfortunately only a shift of the polynom:
det = sqrt(1 + 4*(x0^2 + x0)) / 2
x = c(-1/2 + det, -1/2 - det)
err = 1 + 4*x + 10*x^2 + 10*x^3 + 5*x^4 + x^5
round0(err)
x = x + 1
err = x^5 - x + 1
round0(err)


###
x = (x0^2 + x0)^2 + (x0^2 + x0)
x

poly.calc(x)
err = -4 - 32*x - 77*x^2 - 35*x^3 - 4*x^4 + x^5
round0(err)
#
det = sqrt(1 + 4*(x0^2 + x0)^2 + 4*(x0^2 + x0)) / 2
x = c(-1/2 + det, -1/2 - det)
err = 4 + 12*x + 13*x^2 + 8*x^3 + 5*x^4 + x^5
round0(err)
x = x + 1
x[c(3,4,6,7,10)]
err = 1 - 5*x + 9*x^2 - 2*x^3 + x^5
round0(err)


################
################

# a nice technique to create "awesome" polynomials

n = 5
m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
m = m^(0:(n-1))

# K = 2 is currently fixed value for the first P5!
K = 2 # the 2 derived polynomials work with other values as well;
#
k = K^(1/n)
k = k*m # all 5 (n) roots;
r = k^3 + k^2 + k # currently fixed as well!
# the "derived" polynomial:
x = r
err = - 14 - 10*x - 20*x^2 - 10*x^3 + x^5
round0(err) # all 5 roots are valid;
### Deriving another polynomial based on the previous P5
# [a derivation of the derived polynomial]
# new roots: solutions of x^3 + x^2 + x - r = 0
r = k^3 + k^2 + k
x = sapply(r, function(r) solve.P3(c(1, 1, -r)) )
x = t(x)
x
# only 10 of the 15 values are valid roots!
# [but we have *all* roots!]
err = K^2 + K + 1 + 5*x + 15*x^2 + 30*x^3 + 45*x^4 + (51+K)*x^5 + 45*x^6 + 30*x^7 + 15*x^8 + 5*x^9 + x^10
err = (x^2 + x + 1)^n + K*x^n + K^2 + K
round0(err)

# for n = 4:
# (x^2+x+1)^4 + K*(x+1)^4 - 2*K*x^2 + K^2

# new roots: solutions of x^3 - x^2 + x - r = 0
r = k^3 - k^2 + k
x = sapply(r, function(r) solve.P3(c(-1, 1, -r)) )
x = t(x)
x
# only 10 of the 15 values are valid roots!
# [but we have *all* roots!]
err = K^2 - K + 1 - 5*x + 15*x^2 - 30*x^3 + 45*x^4 - (51-K)*x^5 + 45*x^6 - 30*x^7 + 15*x^8 - 5*x^9 + x^10
err = (x^2 - x + 1)^n + K*x^n + K^2 - K
round0(err)

# for n = 4:
# (x^2-x+1)^4 + K*(x-1)^4 - 2*K*x^2 + K^2


###############

id = (0:4)*2 + 1
m = 2*cos(id * 2 * pi/11)


# Base Polynomial: a Class 3 Polynomial
x = m^2 + m
err = x^5 - 8*x^4 + 19*x^3 - 15*x^2 + x + 1
round0(err)

# the derived's polynomial Derived Polynomial:
det = sqrt(1 + 4*m^2 + 4*m) / 2
x = c(-1/2 + det, -1/2 - det)
x
x = x[c(3,4,6,7,10)]
# 5 of the roots are the real roots
err = 1 - 2*x - 5*x^2 + 2*x^3 + 4*x^4 + x^5
round0(err)
poly.calc(1/(m-1/m)) # same polynomial as the derived!
err = 1 - 2*x - 5*x^2 + 2*x^3 + 4*x^4 + x^5
round0(err)


# Base Polynomial: a Class 3 Polynomial
x = m^2 - m
err = 1 + 3*x - 25*x^2 + 29*x^3 - 10*x^4 + x^5
round0(err)

# the derived's polynomial Derived Polynomial:
det = sqrt(1 + 4*m^2 - 4*m + 0i) / 2
x = c(1/2 + det, 1/2 - det)
x
x = x[c(2,3,4,6,10)]
x
# 5 of the roots are the real roots
err = 1 - 6*x - x^2 + 10*x^3 - 6*x^4 + x^5
round0(err)


###################

library(polynom)
library(pracma)

m5 = complex(re=cos(2*pi/5), im=sin(2*pi/5))
m5 = m5^(1:4)


###
sol = sapply(m5, function(x) roots(c(1, x^3, x)) )
sol

x = sol # *NOT* roots of -1!
1 - x + x^3 - x^4 + x^5 - x^7 + x^8


###
sol = sapply(m5, function(x) roots(c(1, x^4-x^2, x)) )
sol

poly.calc(x)
x = sol
1 - x^2 + 10*x^3 + 11*x^4 - 10*x^5 - x^6 + x^8


###
sol = sapply(m5, function(x) roots(c(1, x^2, -x)) )
sol

x = sol
1 + x + 2*x^2 - 2*x^3 - 3*x^5 + 2*x^6 - x^7 + x^8


##########################
##########################

######################
### P4 Derivations ###

# the P6 variants:
# - are in file Polynomials.Derived.P6.fromP4.R;
#  [the file covers both from P3 & P4]

### other variants

library(polynom)
library(pracma)

### P4 variants

p4der.gen = function(K, coeff=c(1)) {
	r = roots(c(1,coeff,0,0,K))
	if(coeff == 0) {return(list(x=NA, p=NA));}
	x.r = -(r^4 + K) / coeff
	x = sapply(x.r, function(r) roots(c(1,0,0,-r)))
	p1 = round0.p(poly.calc(x))
	p = round0.p(p1 / round0.p(poly.calc(r)))
	return(list(x=x, p=p))
}
p4der2.gen = function(K, coeff=c(1, 0, 1)) {
	r = roots(c(1,coeff,K))
	if(coeff[1] == 0 && coeff[2] == 0) {return(list(x=NA, p=NA));}
	x.r = -(r^4 + coeff[2]*r^2 + K) / coeff[1]
	x = sapply(x.r, function(r) roots(c(1,0,coeff[3],-r)))
	p1 = round0.p(poly.calc(x))
	p = round0.p(p1 / round0.p(poly.calc(r)))
	return(list(x=x, p=p))
}

###
K = 1
p = p4der2.gen(K)
x = p$x
p
round0(predict(p$p, p$x))
err = 4 + 4*x + 5*x^2 + 4*x^3 + 7*x^4 - x^5 + 5*x^6 - x^7 + x^8
round0(err)

###
K = 1
p = p4der.gen(K)
x = p$x
p
err = 1 + 2*x^3 - x^4 + x^6 - x^7 + x^8
round0(err)

###
K = 2
p = p4der.gen(K)
x = p$x
p
err = 4 + 4*x^3 - 2*x^4 + x^6 - x^7 + x^8
round0(err)

###
K = 2
p = p4der.gen(K, coeff=3)
x = p$x
p
err = 4 + 12*x^3 - 2*x^4 + 9*x^6 - 3*x^7 + x^8
round0(err)
# (x^2-x+1) * P6
err = 4 + 4*x + 8*x^3 + 6*x^4 - 2*x^5 + x^6
round0(err)

###
K = 1
p = sapply(-6:6, function(s) print(p4der.gen(K, coeff=s)$p))
p = sapply(-6:6, function(s) print(p4der.gen(s, coeff=s)$p))
#
s = -1
p = p4der.gen(s, coeff=s)
x = p$x
round0( 1 + 2*x^3 + x^4 + x^6 + x^7 + x^8 )


# old [very simple]
r = roots(c(1,0,0,1,1))
x.r = r^4 + r^2
x1 = sapply(x.r, function(r) roots(c(1,0,1,0,-r)))
x = (x1^2)[c(1,3),]
p1 = round0.p(poly.calc(x1))
p2 = round0.p(poly.calc(x))
p1 # = (x^4+x+1)*(x^4-x+1)*(5 + 9*x^2 + 8*x^4 + 4*x^6 + x^8)
p2 # = (x^4 + 2*x^2 - x + 1)*(5 + 9*x + 8*x^2 + 4*x^3 + x^4)


x.r = r^2*(r+1)
x = sapply(x.r, function(r) roots(c(1,1,0,-r)))
p1 = round0.p(poly.calc(x))
p1
round0.p(p1 / round0.p(poly.calc(r)))
1 - x - 4*x^2 - x^3 + 5*x^4 + 6*x^5 + 6*x^6 + 4*x^7 + x^8


x.r = -(r+1)
x = sapply(x.r, function(r) roots(c(1,0,0,0,-r)))
p1 = round0.p(poly.calc(x))
p1
round0.p(p1 / round0.p(poly.calc(r)))
1 - x + x^2 - x^3 + 3*x^4 - 2*x^5 + x^6 + 3*x^8 - x^9 + x^12



### P5
d = 2
c = 1
r = roots(c(1,0,-5*c,0,5*c^2, -2*d))
x.r = (r^5 + 5*c^2*r - 2*d)/5/c
x = sapply(x.r, function(r) roots(c(1,0,0,-r)))
p1 = round0.p(poly.calc(x))
p1
round0.p(p1 / round0.p(poly.calc(r)))
err = 16 + 20*x + 25*x^2 + 40*x^3 + 25*x^4 + 4*x^5 + 20*x^6 + 5*x^8 + x^10
round0(err)


### P5: x^5 = x+1
b1 = -1
b0 = -1
r = roots(c(1,0,0,0, b1, b0))
x.r = -(b1*r + b0)
x = sapply(x.r, function(r) roots(c(1,0,0,0,0,-r)))
p1 = round0.p(poly.calc(x))
p1
round0.p(p1 / round0.p(poly.calc(r)))
err = 1 - x + x^2 - x^3 + x^4 - 4*x^5 + 3*x^6 - 2*x^7 + x^8 + 6*x^10 - 3*x^11 + x^12 - 4*x^15 + x^16 + x^20
round0(err)


######################
######################

### Analysis

# - various experimental analysis;

### x^5 - x^2 + k = 0
# - 3 real roots;
curve(x^5, from=-1.2, to=1.2)
curve(x^2 - 3/5 * (4/25)^(1/3), add=T, col="green")
curve(-x^2 + 3/5 * (4/25)^(1/3), add=T, col="red")
curve(x^2, add=T, col="orange")


#########################
#########################

### Experimental

### x^5 - x = R

### roots: r1, Conj(r1) =>
(a1+b1*1i)^5 - (a1+b1*1i) - R # = 0
(a1-b1*1i)^5 - (a1-b1*1i) - R # = 0

### Diff =>
5*a1^4 - 10*a1^2*b1^2 + b1^4 - 1 # = 0
5*a2^4 - 10*a2^2*b2^2 + b2^4 - 1 # = 0

### Sum =>
a1^5 - 10*a1^3*b1^2 + 5*a1*b1^4 - a1 - R # = 0
a2^5 - 10*a2^3*b2^2 + 5*a2*b2^4 - a2 - R # = 0

### Diff Eq(r1) - Eq(r2):
(a1+b1*1i)^4 + (a2+b2*1i)^4 + (a1+b1*1i)*(a2+b2*1i)*(a1^2+a2^2-b1^2-b2^2 + 2*(a1*b1+a2*b2)*1i) +
	+ (a1+b1*1i)^2*(a2+b2*1i)^2 - 1 # = 0
a1^4 + a2^4 + b1^4 + b2^4 - 6*(a1^2*b1^2 + a2^2*b2^2) + 4*(a1^3*b1 + a2^3*b2 - a1*b1^3 - a2*b2^3)*1i +
	+ (a1*a2 - b1*b2 + (a1*b2+a2*b1)*1i)*(a1^2+a2^2-b1^2-b2^2 + 2*(a1*b1+a2*b2)*1i) +
	+ (a1*a2 - b1*b2 + (a1*b2+a2*b1)*1i)^2 - 1 # = 0
# Re =>
a1^4 + a2^4 + b1^4 + b2^4 - 6*(a1^2*b1^2 + a2^2*b2^2) +
	+ (a1*a2 - b1*b2)*(a1^2+a2^2 - (b1^2+b2^2)) - 2*(a1*b2 + a2*b1)*(a1*b1 + a2*b2) +
	+ (a1*a2 - b1*b2)^2 - (a1*b2 + a2*b1)^2 - 1 # = 0
# TODO + Im()

# [alternative] =>
a1^5 + 5*a1^4*b1*1i - 10*a1^3*b1^2 - 10*a1^2*b1^3*1i + 5*a1*b1^4 + b1^5*1i - (a1+b1*1i) +
	- a2^5 - 5*a2^4*b2*1i + 10*a2^3*b2^2 + 10*a2^2*b2^3*1i - 5*a2*b2^4 - b2^5*1i + (a2+b2*1i) # = 0
a1^5 - a2^5 - 10*a1^3*b1^2 + 10*a2^3*b2^2 + 5*a1*b1^4 - 5*a2*b2^4 - (a1-a2) # = 0


### Debug:
R = 1
x = roots(c(1,0,0,0,-1,-R));
a1 = Re(x[2]); b1 = Im(x[2]);
a2 = Re(x[4]); b2 = Im(x[4]);

