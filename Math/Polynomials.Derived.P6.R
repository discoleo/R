
########################
###
### Leonard Mada
### [the one and only]
###
### Derived Polynomials: P6
### P6 Polynomials
###
### draft v.0.1d



### History
# draft v.0.1b - v.0.1d:
# - entanglements with:
#   m[3] (m^3 = 1), cos(2*pi/5), cos(2*pi/7);
# - new variants:
#   3 - 3*x + x^6 = 0;
#     factored as 2 cubics:
#     (x^3 + (m-m^2)*x^2 +(m+2*m^2)*x - 2*m-m^2)*
#     (x^3 - (m-m^2)*x^2 +(m^2+2*m)*x - 2*m^2-m)
#   7 + 14*x + 7*x^2 + x^6 = 0;
#   many other, including:
#   -1 - x + 2*x^3 - 2*x^5 + x^6 = 0;
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


####################

### TODO:
# - use exact P3 solver for P6, when applicable;
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
solve.p6 = function(coeff, type=110) {
	m = unity(3, all=FALSE)
	a = coeff[1]; b = coeff[2]; c = coeff[3]
	# TODO: all permutations
	# TODO: automate permutations;
	if(type == 110) {
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


#######################

####################
### Convolved P3 ###
### Polynomials  ###
####################


########
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


#######################
#######################

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

x = unlist(lapply(c3, function(x) roots(c(1,-1/x+2*x,x^2)) ))
1 + 7*x + 6*x^2 - 2*x^4 + x^6

x = sapply(1:3, function(id) roots(c(1, 2+c3[id]-c3[id]^2, c3[id])))
1 - 9*x^2 - 14*x^3 - 8*x^4 + x^6

### Symmetric
x = sapply(1:3, function(id) roots(c(1, 3 - c3[id] - 2*c3[id]^2, 1)))
1 - 4*x^2 + 7*x^3 - 4*x^4 + x^6
x = sapply(1:3, function(id) roots(c(1, 10+10*c3[id]-8*c3[id]^2-5*c3[id]^3, 1)))
1 - 18*x^2 - 7*x^3 - 18*x^4 + x^6

### Asymmetric
x = sapply(1:3, function(id) roots(c(1, -10*c3[id]+2*c3[id]^2+5*c3[id]^3, 3-c3[id]-c3[id]^2)))
1 - x^2 - 14*x^3 - 16*x^4 + x^6

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
x = sapply(1:3, function(id) roots(c(1, 2+c3[id]-c3[id]^2, 2-c3[id])))
7 + 14*x + 7*x^2 + x^6


############

            | c1*m^2+c2*m      | b1*m^2+b2*m      | a1*m^2+a2*m      | 1
c1*m+c2*m^2 | c1^2+c2^2-c1*c2  | b1*c1+b2*c2+     | a1*c1+a2*c2+     | c1*m+c2*m^2
            |                  | b1*c2*m+b2*c1*m^2| a1*c2*m+a2*c1*m^2|
b1*m+b2*m^2 | b1*c1+b2*c2+     | b1^2+b2^2-b1*b2  | a1*b1+a2*b2+     | b1*m+b2*m^2
            | b1*c2*m^2+b2*c1*m|                  | a1*b2*m+a2*b1*m^2|
a1*m+a2*m^2 | a1*c1+a2*c2+     | a1*b1+a2*b2+     | a1^2+a2^2-a1*a2  | a1*m+a2*m^2
            | a1*c2*m^2+a2*c1*m| a1*b2*m^2+a2*b1*m|                  |
1           | c1*m^2+c2*m      | b1*m^2+b2*m      | a1*m^2+a2*m      | 1


x^6 - (a1+a2)*x^5 + (a1^2+a2^2-a1*a2-b1-b2)*x^4 + (2*a1*b1+2*a2*b2-a1*b2-a2*b1-c1-c2)*x^3 +
(b1^2+b2^2-b1*b2+2*a1*c1+2*a2*c2-a1*c2-a2*c1)*x^2 + (2*b1*c1+2*b2*c2-b1*c2-b2*c1)*x + c1^2+c2^2-c1*c2

coeff = c(1, - (a1+a2), (a1^2+a2^2-a1*a2-b1-b2), (2*a1*b1+2*a2*b2-a1*b2-a2*b1-c1-c2),
(b1^2+b2^2-b1*b2+2*a1*c1+2*a2*c2-a1*c2-a2*c1), (2*b1*c1+2*b2*c2-b1*c2-b2*c1), c1^2+c2^2-c1*c2)


###
coeff = c(1,-1, 1,1, 1,1)
sol = solve.p6(coeff, type=222)
sol

x = sol$x
err = 1 + 2*x + x^2 - 2*x^3 + x^4 + x^6
round0(err)


###
coeff = c(1,-1, 1,2, -2,-1)
sol = solve.p6(coeff, type=222)
sol

x = sol$x
err = 3 - 3*x + x^6
round0(err)


################


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


# c(1,2,-2), c(1,3,1), c(3,-1,-1), c(4,6,-7), c(1,3,2), c(0,2,1), c(1,0,-1)
# c(-2,-3,-1)
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
