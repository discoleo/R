
########################
###
### Leonard Mada
### [the one and only]
###
### Derived Polynomials
### v.0.3x

### Note:
# This is the 1st part towards:
# - introducing various derived polynomials,
#   e.g. based on polynomials of Class 1;
# - brief introduction of polynomials of Class 1;
# - introducing interesting properties of polynomial roots;
# - includes a different approach to polynomials.

### History
# v.0.3.x:
# - moved P6 poynomials to:
#   https://github.com/discoleo/R/blob/master/Math/Polynomials.Derived.P6.R
# v.0.3.c - v.0.3.f:
# - more awesome polynomials with complete solutions:
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

# needed to get the roots of the base polynomial;
# the derived polynomial is currently derived "empirically";
# [it is relatively easy to derive parametrically as well]

######################

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


###############
### Decomposing
### derived polynomial

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

# a nice technique to create "awsome" polynomials

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

