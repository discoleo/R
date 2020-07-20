
########################
###
### Leonard Mada
### [the one and only]
###
### Derived Polynomials: P6
### P6 Polynomial
###
### draft v.0.1


### History
# draft v.0.1
# - moved P6 from Polynomials.Derived.R
#   in a separate file;
# - old version (up to draft v.0.3f):
#   https://github.com/discoleo/R/blob/master/Math/Polynomials.Derived.R
# v.0.3.c - v.0.3.f [initial file]:
# - more awesome polynomials with complete solutions:
#   1 + x - x^4 - x^5 + x^6 = 0;
#   1 - x + x^2 + x^3 + x^4 + x^6 = 0;
#   1 + x + x^2 - x^3 - 2*x^4 + x^6 = 0;
#   1 + 2*x + x^2 - x^3 - x^4 + x^6 = 0;
#   & many more;


####################

### TODO:
# - use exact P3 solver for P6, when applicable;
library(pracma)
m = unity(3, all=FALSE)

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
coeff = c(1, 1, -1)
solve.p6(coeff, type=110)

###
coeff = c(-1, 2, -1)
solve.p6(coeff, type=110)


########
### Type = (x^3, m, m, m)
coeff = c(1, 1, -1)
solve.p6(coeff, type=111)

### type = (x^3, m, m, m)
coeff = c(1, 1, 1)
solve.p6(coeff, type=111)


### type = (x^3, m, m, m)
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


### Example 1: coeffs only +/-1, no b5
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

### Example 2b: trivial: (x^2 + x + 1) * ...
x = c(solve.p3(c(m, m)) * m, solve.p3(c(m^2, m^2)) * m^2)
1 - x + x^2 - x^3 + 2*x^4 + x^6

# needs root shift or library(pracma)
x = c(roots(c(1,m^2,-1,1)) * m^2, roots(c(1,m,-1,1)) * m)
1 + x + 2*x^4 - x^5 + x^6
1 - x + 2*x^4 + x^5 + x^6 # -x


##############
### Example 3: 2*x^3
x = c(solve.p3(c(-m, 1)) * m, solve.p3(c(-m^2, 1)) * m^2)
1 + x + x^2 + 2*x^3 + x^4 + x^6
1 - x + x^2 - 2*x^3 + x^4 + x^6 # -x
x = 1/x
1 + x^2 + 2*x^3 + x^4 + x^5 + x^6

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
x = c(solve.p3(c(-m^2, m^2)) * m, solve.p3(c(-m, m)) * m^2)
1 - 2*x + x^2 - x^3 + x^4 + x^6

### Example 4b:
x = c(solve.p3(c(-m, -m)) * m^2, solve.p3(c(-m^2, -m^2)) * m)
1 + 2*x + x^2 + x^3 + x^4 + x^6 # trivial from above

### Example 4c: (x^2 + x + 1) * ...
x = c(solve.p3(c(1, m)), solve.p3(c(1, m^2)))
1 + 2*x + x^2 - x^3 - x^4 + x^6


##############
### Example 5: multiple coeffs of 2
x = c(solve.p3(c(m, 1)) * m^2, solve.p3(c(m^2, 1)) * m^2)
1 + 2*x + x^2 + 2*x^3 + 2*x^4 + x^6

x = c(solve.p3(c(-m^2, -1)), solve.p3(c(-m, -1)))
1 + 2*x + x^2 - 2*x^3 - 2*x^4 + x^6

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


