########################
##
## Leonard Mada
## [the one and only]
##
## P6 Polynomials:
## Symmetric Polynomials
##
## draft v.0.1d


### Decomposition of Symmetric Polynomials of order [2n]
### into Polynomials of order [n]

# - every *strictly* symmetric polynomial of order 2*n
#   can be "decomposed" into polynomials of order n;
# - strictly symmetric:
#   1 + b1*x + b2*x^2 + ... + b2*x^(n-2) + b1*x^(n-1) + x^n;
#   b[j] = b[n-j] ***and*** b0 == 1;
# - "strictly" symmetric polynomials of order 2*n+1
#   have the trivial root x = -1 and can be factored into:
#   (x+1)*P[2n], where P[2n] is strictly symmetric;


###############

###############
### History ###

# draft v.0.1c:
# - generalized quasi-symmetric:
#   x^6 + b[1]*x^5 + b[2]*x^4 + b[3]*x^3 + b[2]*R*x^2 + b[1]*R^2*x + R^3 = 0;
# draft v.0.1b:
# - added the decompositions for P8 & P10;
# - the P[5] (for P10) is solvable numerically;
# - TODO: the "minus" versions as well;
# draft v.0.1a:
# - new file, moved from Polynomials.Derived.P6.R
#   (during the last sub-version of v.0.3b);



################

library(polynom)
library(pracma)

### Helper Functions

# Roots of Unity
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
### Solver: simplified
solve.p6sym = function(b, type="symmetric") {
	# b = symmetric coefficients
	# 1 + b1*x + b2*x^2 + b3*x^3 + b2*x^4 + b1*x^5 + x^6
	len = length(b)
	if(len == 3) {
		# P[6]: only Strictly Symmetric;
		if(type == "minus") {
			r = roots(c(1, b[1], -b[2]+3, b[3] + 2*b[1]))
			b0 = -1;
		} else {
			r = roots(c(1, b[1], b[2]-3, b[3] - 2*b[1]))
			b0 = 1;
		}
		x = sapply(r, function(r) roots(c(1,-r,b0)))
		p = round0.p(poly.calc(x))
		b = c(b0, b, rev(b[-len]), 1)
		if(type == "minus") {
			b[5] = -b[5]; # b4
		}
		id = 0:(2*len)
		err = sapply(x, function(x) round0(sum(b*x^id)))
		return(list(x=x, p=p, err=err))
	} else if(len == 4) {
		if(type == "minus") {
			# TODO
		} else {
			r = roots(c(1, b[1], b[2]-4, b[3]-3*b[1], b[4]-2*b[2]+2))
			b0 = 1;
		}
	} else if(len == 5) {
		if(type == "minus") {
			# TODO
		} else {
			r = roots(c(1, b[1], b[2]-5, b[3] - 4*b[1], b[4] - 3*b[2] + 5, b[5] + 2*b[1] - 2*b[3]))
			b0 = 1;
		}
	} else {
		print("Not yet implemented!")
	}
	# x + 1/x = r;
	x = sapply(r, function(r) roots(c(1,-r,b0)));
	p = round0.p(poly.calc(x))
	b = c(b0, b, rev(b[-len]), 1)
	if(type == "minus") {
		b[5] = -b[5]; # b4??? # TODO
	}
	id = 0:(2*len)
	err = sapply(x, function(x) round0(sum(b*x^id)))
	return(list(x=x, p=p, err=err))
}

################

### Strictly Symmetric P[6]

### Examples:

###
b = c(-2, 0, -1)
p = solve.p6sym(b)
p; x = p$x;
1 - 2*x - x^3 - 2*x^5 + x^6

###
b = c(-2, 3, 2)
p = solve.p6sym(b)
p; x = p$x;
1 - 2*x + 3*x^2 + 2*x^3 + 3*x^4 - 2*x^5 + x^6

###
b = c(1, 3, 3)
p = solve.p6sym(b)
p; x = p$x;
1 + x + 3*x^2 + 3*x^3 + 3*x^4 + x^5 + x^6

###
b = c(1, 0, 0)
p = solve.p6sym(b, type = "minus")
p; x = p$x;
-1 + x + x^5 + x^6

##############

##############
### Series ###

### 1 + b1*x + b1*x^5 + x^6
p = sapply(-5:5, \(b) print(solve.p6sym(c(b, 0, 0))$p))
1 - 5*x - 5*x^5 + x^6 
1 - 4*x - 4*x^5 + x^6 
1 - 3*x - 3*x^5 + x^6 
1 - 2*x - 2*x^5 + x^6 
1 - x - x^5 + x^6 
1 + x^6 
1 + x + x^5 + x^6 
1 + 2*x + 2*x^5 + x^6 
1 + 3*x + 3*x^5 + x^6 
1 + 4*x + 4*x^5 + x^6 
1 + 5*x + 5*x^5 + x^6

### e.g. - x^2 - x^4
p = sapply(-5:5, \(b) print(solve.p6sym(c(-1, -1, b))$p))
1 - x - x^2 - 5*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 - 4*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 - 3*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 - 2*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 - x^3 - x^4 - x^5 + x^6 
1 - x - x^2 - 0 - x^4 - x^5 + x^6 
1 - x - x^2 + x^3 - x^4 - x^5 + x^6 
1 - x - x^2 + 2*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 + 3*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 + 4*x^3 - x^4 - x^5 + x^6 
1 - x - x^2 + 5*x^3 - x^4 - x^5 + x^6
### e.g. 0*(x^2 + x^4)
p = sapply(-5:5, function(b) print(solve.p6sym(c(-1, 0, b))$p))
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
### e.g. +(x^2 + x^4)
p = sapply(-5:5, function(b) print(solve.p6sym(c(-1, 1, b))$p))
1 - x + x^2 - 5*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 - 4*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 - 3*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 - 2*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 - x^3 + x^4 - x^5 + x^6 
1 - x + x^2 - 0 + x^4 - x^5 + x^6 
1 - x + x^2 + x^3 + x^4 - x^5 + x^6 
1 - x + x^2 + 2*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 + 3*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 + 4*x^3 + x^4 - x^5 + x^6 
1 - x + x^2 + 5*x^3 + x^4 - x^5 + x^6


### Minus derivation of *Strictly* Symmetric
# -1 + b1*x + b2*x^2 + b3*x^3 - b2*x^4 + b1*x^5 + x^6

# [in the old version] [deprecated]
# b1 => b1*1i, b3 => b3*-1i; x = -1i * x;
# p6sq.gen(c(1, b1 + 6, 4*b1 + 9 + b2, 2*b1 + 2 + 2*b2 + b3)


### "Minus"-series: +b2*x^2 - b2*x^4 & b0 = -1;
p = sapply(-5:5, function(b) print(solve.p6sym(c(-1, 1, b), type="minus")$p))
-1 - x + x^2 - 5*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 - 4*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 - 3*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 - 2*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 - x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 - 0 - x^4 - x^5 + x^6 
-1 - x + x^2 + x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 + 2*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 + 3*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 + 4*x^3 - x^4 - x^5 + x^6 
-1 - x + x^2 + 5*x^3 - x^4 - x^5 + x^6


###
p = solve.p6sym(c(1,2,2), type="minus")
x = p$x
-1 + x + 2*x^2 + 2*x^3 - 2*x^4 + x^5 + x^6

###
p = solve.p6sym(c(1,-2,2), type="minus")
x = p$x
-1 + x - 2*x^2 + 2*x^3 + 2*x^4 + x^5 + x^6

###
p = solve.p6sym(c(1,3,3), type="minus")
x = p$x
-1 + x + 3*x^2 + 3*x^3 - 3*x^4 + x^5 + x^6

###
p = solve.p6sym(c(1,-3,3), type="minus")
x = p$x
-1 + x - 3*x^2 + 3*x^3 + 3*x^4 + x^5 + x^6

###
p = solve.p6sym(c(-2,-3,5), type="minus")
x = p$x
-1 - 2*x - 3*x^2 + 5*x^3 + 3*x^4 - 2*x^5 + x^6


##########################
##########################

### x^6 + b1*x^5 + b2*x^4 + b3*x^3 + b2*R*x^2 + b1*R^2*x + R^3 = 0

### Solution:
x^6 + b1*x^5 + b2*x^4 + b3*x^3 + b2*R*x^2 + b1*R^2*x + R^3 # = 0 / x^3
x^3 + (R/x)^3 + b1*(x^2 + (R/x)^2) + b2*(x + R/x) + b3 # = 0
### R/x = y =>
x^3 + y^3 + b1*(x^2 + y^2) + b2*(x + y) + b3 # = 0
### S = x + y =>
S^3 - 3*x*y*S + b1*(S^2 - 2*x*y) + b2*S + b3
S^3 - 3*R*S + b1*(S^2 - 2*R) + b2*S + b3
S^3 + b1*S^2 + (b2 - 3*R)*S - 2*b1*R + b3

### Step 2:
# solve:
# x + y = S
# x*y = R

solve.P6Sym = function(R, b) {
	coeff = c(1, b[1], (b[2] - 3*R[1]), - 2*b[1]*R[1] + b[3])
	S = roots(coeff)
	xy.d = sqrt(S^2 - 4*R[1] + 0i)
	x = (S + xy.d) / 2;
	y = (S - xy.d) / 2;
	return(c(x, y))
}

### Examples:
R = 2
b = c(1,-1,0)
x = solve.P6Sym(R, b)

err = x^6 + b[1]*x^5 + b[2]*x^4 + b[3]*x^3 + b[2]*R*x^2 + b[1]*R^2*x + R^3
round0(err)


### Ex 2:
R = 2
b = c(1,-3,1)
x = solve.P6Sym(R, b)

err = x^6 + b[1]*x^5 + b[2]*x^4 + b[3]*x^3 + b[2]*R*x^2 + b[1]*R^2*x + R^3
round0(err)


### Ex 3:
R = 5
b = c(1, 0, 0)
x = solve.P6Sym(R, b)

err = x^6 + b[1]*x^5 + b[2]*x^4 + b[3]*x^3 + b[2]*R*x^2 + b[1]*R^2*x + R^3
round0(err)

### Special
m = complex(re=cos(pi/3), im=sin(pi/3))

R = m
b = c(m^2, -m, -1)
x = solve.P6Sym(R, b) / m^2

err = -1 + x - x^2 - x^3 + x^4 + x^5 + x^6
round0(err)


### Ex 2: the same
R = 1
b = c(1i, -1, 1i)
x = solve.P6Sym(R, b) / 1i

err = -1 + x - x^2 - x^3 + x^4 + x^5 + x^6
round0(err)


########################
########################

### P2 o P3 Entanglement
### "Skewed"

# but quasi-trivial;

(x^2 + b0/x + b1)^2 + c1*(x^2 + b0/x + b1) + c0 # = 0
# omitted: b2*x (for simplicity);
x^4 + b0^2/x^2 + b1^2 + 2*b1*x^2 + 2*b0*x + 2*b0*b1/x + c1*x^2 + b0*c1/x + b1*c1 + c0 # = 0
x^4 + (2*b1 + c1)*x^2 + 2*b0*x + b1*c1 + b1^2 + c0 + 2*b0*b1/x + b0*c1/x + b0^2/x^2 # = 0
x^6 + (2*b1 + c1)*x^4 + 2*b0*x^3 + (b1*c1 + b1^2 + c0)*x^2 + b0*(2*b1+c1)*x + b0^2 # = 0

### Example:
# c0 = - b1*(b1 + c1)
# but is trivial!
x^6 + (2*b1 + c1)*x^4 + 2*b0*x^3 + b0*(2*b1+c1)*x + b0^2 # = 0

### Ex 1:
b0 = 3;
b1 = 1; c1 = -1;
c0 = - b1*(b1 + c1)
#
r = roots(c(1, c1, c0))
x = sapply(r, function(r) roots(c(1, 0, b1-r, b0)))

x^6 + x^4 + 2*b0*x^3 + b0*x + b0^2 # = 0


####################
####################


####################
### Symmetric P8 ###

### Generalized

### Example:
b = c(-1, 3, 2, -2); b1=b[1]; b2=b[2]; b3=b[3]; b4=b[4];
c = sqrt(3);
x = roots(c(1,b1,b2,b3,b4,c*b3,c^2*b2,c^3*b1, c^4));

x^8 + b1*x^7 + b2*x^6 + b3*x^5 + b4*x^4 + b3*c*x^3 + b2*c^2*x^2 + b1*c^3*x + c^4 # = 0

# Transformed Polynomial:
z = x + c/x;
z^4 - 4*c*(z^2 - 2*c) - 6*c^2 + b1*(z^3 - 3*c*z) + b2*(z^2 - 2*c) + b3*z + b4 # = 0
z^4 + b1*z^3 + (b2 - 4*c)*z^2 + (b3 - 3*c*b1)*z + b4 - 2*c*b2 + 2*c^2 # = 0


### Old Examples

###
b = c(1, 0, 0, 1)
# add the function parameter to the list: e.g. b[2] + bx;
p = sapply(-5:5, function(bx) print(solve.p6sym(c(b[1], b[2] + bx, b[3], b[4]))$p))
1 + x - 5*x^2 + x^4 - 5*x^6 + x^7 + x^8 
1 + x - 4*x^2 + x^4 - 4*x^6 + x^7 + x^8 
1 + x - 3*x^2 + x^4 - 3*x^6 + x^7 + x^8 
1 + x - 2*x^2 + x^4 - 2*x^6 + x^7 + x^8 
1 + x - x^2 + x^4 - x^6 + x^7 + x^8 
1 + x - 0 + x^4 - 0 + x^7 + x^8 
1 + x + x^2 + x^4 + x^6 + x^7 + x^8 
1 + x + 2*x^2 + x^4 + 2*x^6 + x^7 + x^8 
1 + x + 3*x^2 + x^4 + 3*x^6 + x^7 + x^8 
1 + x + 4*x^2 + x^4 + 4*x^6 + x^7 + x^8 
1 + x + 5*x^2 + x^4 + 5*x^6 + x^7 + x^8


### more Examples
b = c(1, 1, -1, -1)
p = solve.p6sym(b)
p
x = p$x
1 + b[1]*x + b[2]*x^2 + b[3]*x^3 + b[4]*x^4 + b[3]*x^5 + b[2]*x^6 + b[1]*x^7 + x^8


###
b = c(1, 1, 2, -1)
p = solve.p6sym(b)
p
x = p$x
1 + b[1]*x + b[2]*x^2 + b[3]*x^3 + b[4]*x^4 + b[3]*x^5 + b[2]*x^6 + b[1]*x^7 + x^8



#####################
#####################

#####################
### Symmetric P10 ###

# Solution: the P5 is solved numerically!
# coeffs of P5:
# c(1, b[1], b[2]-5, b[3] - 4*b[1], b[4] - 3*b[2] + 5, b[5] + 2*b[1] - 2*b[3])


### more Examples
b = c(10, 1, -1, -1, 0)
p = solve.p6sym(b)
p
x = p$x
1 + b[1]*x + b[2]*x^2 + b[3]*x^3 + b[4]*x^4 + b[5]*x^5 + b[4]*x^6 + b[3]*x^7 + b[2]*x^8 + b[1]*x^9 + x^10


###
b = c(1, 1, 2, -1, -1)
p = solve.p6sym(b)
p
x = p$x
1 + b[1]*x + b[2]*x^2 + b[3]*x^3 + b[4]*x^4 + b[5]*x^5 + b[4]*x^6 + b[3]*x^7 + b[2]*x^8 + b[1]*x^9 + x^10


b = c(1, 0, 0, 0, 0)
# add the function parameter to the list: e.g. b[2] + bx;
p = sapply(-5:5, function(bx) print(solve.p6sym(c(b[1], b[2] + bx, b[3], b[4], b[5]))$p))
1 + x - 5*x^2 - 5*x^8 + x^9 + x^10
1 + x - 4*x^2 - 4*x^8 + x^9 + x^10
1 + x - 3*x^2 - 3*x^8 + x^9 + x^10
1 + x - 2*x^2 - 2*x^8 + x^9 + x^10
1 + x - x^2 - x^8 + x^9 + x^10
1 + x - 0 - 0 + x^9 + x^10
1 + x + x^2 + x^8 + x^9 + x^10
1 + x + 2*x^2 + 2*x^8 + x^9 + x^10
1 + x + 3*x^2 + 3*x^8 + x^9 + x^10
1 + x + 4*x^2 + 4*x^8 + x^9 + x^10
1 + x + 5*x^2 + 5*x^8 + x^9 + x^10


# if( b[c(1, 3, 5)] == 0) => actually trivial: r = +/- 1i;
# 1 + b2*x^2 + b4*x^4 + b4*x^6 + b2*x^8 + x^10 = 0;
# P5: c(1, 0, b[2]-5, 0, b[4] - 3*b[2] + 5, 0)


#########################
#########################


# experimental: power 3
p6sq3.old.gen = function(p3.coeff, mult=1, b0=-1, asSq=TRUE) {
	m = unity(3, all=F)
	r = c(roots(p3.coeff))
	r = ifelse(Im(r) == 0 & Re(r) < 0, - (-Re(r))^(1/3), r^(1/3))
	x1 = c(r, m*r, m^2*r)
	# x1 = c(sqrt(r+0i), -(sqrt(r+0i)))
	x = sapply(x1, function(r) roots(c(1, mult*r, 0, b0)))
	p = poly.calc(x)
	for(i in 1:length(p)) p[[i]] = round0(p[[i]])
	if(asSq) {
		x = x[1,]^2 # the squares => P6
		len = length(p3.coeff)*4 - 3
		p = polynomial(p[seq(from=1, to=len, by=2)])
	}
	return(list(x=x, p=p))
}
# new version
p6sq3.gen = function(p3.coeff, mult=1, b0=-1, asSq=T) {
	# m = unity(3, all=F)
	r = c(roots(p3.coeff))
	x = sapply(r, function(r) roots(c(1, -3*b0 + mult*r, 3*b0^2 + mult*r, -b0^3)))
	p = poly.calc(x)
	for(i in 1:length(p)) p[[i]] = round0(p[[i]])
	return(list(x=x, p=p))
}

# p6sq.gen(c(1, b1 + 6, 4*b1 + 9 + b2, 2*b1 + 2 + 2*b2 + b3)

### Test
p = p6sq3.gen(c(1,-1,1,-1), asSq=F)
p

p = p6sq3.gen(c(1,0,1,0,-1), asSq=F)
polynomial(p$p) / polynomial(c(1,1))^4
1 + 8*x + 29*x^2 + 60*x^3 + 75*x^4 + 60*x^5 + 29*x^6 + 8*x^7 + x^8

###
p = p6sq3.gen(c(1,-15,30))
x = p$x
1 + 21*x + 105*x^2 + 170*x^3 + 105*x^4 + 21*x^5 + x^6

########################
########################


### Old
### TODO: cleanup


### Strictly Symmetric P6


#######################
### *Strictly* Symmetic
# b0 == 1
# 1 + b1*x + b2*x^2 + b3*x^3 + b2*x^4 + b1*x^5 + x^6
# new simplified version: solve.p6sym(c(b1, b2, b3))
# old variant:
# p6sq.gen(c(1, b1 + 6, 4*b1 + 9 + b2, 2*b1 + 2 + 2*b2 + b3)


### Examples

p = sapply(-6:6, function(b) print(solve.p6sym(c(b, 0, 0))$p))
# the values of the roots must be extracted explicitly from solve.p6sym();
1 - 5*x - 5*x^5 + x^6
1 - 4*x - 4*x^5 + x^6
1 - 3*x - 3*x^5 + x^6
1 - 2*x - 2*x^5 + x^6
1 - x - x^5 + x^6
1 + x^6
1 + x + x^5 + x^6 # basically -x;
1 + 2*x + 2*x^5 + x^6
1 + 3*x + 3*x^5 + x^6
1 + 4*x + 4*x^5 + x^6
1 + 5*x + 5*x^5 + x^6
p = sapply(-6:6, function(b) print(solve.p6sym(c(-1, 0, b))$p))
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

series.id = -5:5
p = sapply(series.id, function(b) print(solve.p6sym(c(0,-1, b))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(1, 1, b))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(5, b, -2))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(5, b, -3))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(-2, -2, b))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(-3, -3, b))$p))
p = sapply(series.id, function(b) print(solve.p6sym(c(-3, -4, b))$p))


### 1 - x - x^2 - x^4 - x^5 + x^6
p = solve.p6sym(c(-1,-1, 0))
p
x = p$x
1 - x - x^2 - x^4 - x^5 + x^6


### [old workout]
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


#################


### TODO: cleanup

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


###########################
###########################

#######################
### Symmetrical P12 ###

m1 = complex(re=cos(pi/3), im=sin(pi/3)) # -1!
m3 = complex(re=cos(2*pi/3), im=sin(2*pi/3))
c2 = 2*cos(2*pi/5 * 1:2)
c3 = 2*cos(2*pi/7 * 1:3)

###
r = roots(c(1,3, 0,0,0, 3,1))
x = c(r * m1, r/m1)
poly.calc(x)
1 + 3*x + 9*x^2 + 3*x^5 - 7*x^6 + 3*x^7 + 9*x^10 + 3*x^11 + x^12

###
r1 = roots(c(1,3*m1, 0,0,0, 3*m1,1))
r2 = roots(c(1,3/m1, 0,0,0, 3/m1,1))
x = c(r1, r2)
poly.calc(x)
1 + 3*x + 9*x^2 + 3*x^5 + 20*x^6 + 3*x^7 + 9*x^10 + 3*x^11 + x^12


### Mixing m1 with m-1
r1 = roots(c(1,3*m1, 2*m3,0,2*m3, 3*m1,1))
r2 = roots(c(1,3/m1, 2/m3,0,2/m3, 3/m1,1))
x = c(r1, r2)
poly.calc(x)
1 + 3*x + 7*x^2 + 6*x^3 + 2*x^4 + 9*x^5 + 28*x^6 + 9*x^7 + 2*x^8 + 6*x^9 + 7*x^10 + 3*x^11 + x^12

### Mixing m1 with m-1
r1 = roots(c(1,3*m1, 2/m3,0,2/m3, 3*m1,1))
r2 = roots(c(1,3/m1, 2*m3,0,2*m3, 3/m1,1))
x = c(r1, r2)
poly.calc(x)
1 + 3*x + 7*x^2 - 12*x^3 + 2*x^4 - 9*x^5 + 28*x^6 - 9*x^7 + 2*x^8 - 12*x^9 + 7*x^10 + 3*x^11 + x^12

### Only m-1
r1 = roots(c(1,3*m1, 2/m1,0,2/m1, 3*m1,1))
r2 = roots(c(1,3/m1, 2*m1,0,2*m1, 3/m1,1))
x = c(r1, r2)
poly.calc(x)
1 + 3*x + 11*x^2 - 6*x^3 + 6*x^4 - 3*x^5 + 28*x^6 - 3*x^7 + 6*x^8 - 6*x^9 + 11*x^10 + 3*x^11 + x^12


### Sqrt()
r1 = roots(c(1, 3*sqrt(2), 2-sqrt(2),0,2-sqrt(2), 3*sqrt(2),1))
r2 = roots(c(1,-3*sqrt(2), 2+sqrt(2),0,2+sqrt(2),-3*sqrt(2),1))
x = c(r1, r2)
poly.calc(x)
1 - 14*x^2 + 12*x^3 + 6*x^4 + 12*x^5 - 30*x^6 + 12*x^7 + 6*x^8 + 12*x^9 - 14*x^10 + x^12


### Sqrt()
r1 = roots(c(1, 1i,-1i, 0,-1i, 1i,1))
r2 = roots(c(1,-1i, 1i, 0, 1i,-1i,1))
x = c(r1, r2)
poly.calc(x)
1 + x^2 - 2*x^3 + x^4 - 2*x^5 + 6*x^6 - 2*x^7 + x^8 - 2*x^9 + x^10 + x^12


### cos(2*pi/5)
r1 = roots(c(1, c2[1], 0,0,0, c2[1],1))
r2 = roots(c(1, c2[2], 0,0,0, c2[2],1))
x = c(r1, r2)
poly.calc(x)
1 - x - x^2 - x^5 - x^7 - x^10 - x^11 + x^12


### cos(2*pi/7)
r1 = roots(c(1, c3[1], 0,0,0, c3[1],1))
r2 = roots(c(1, c3[2], 0,0,0, c3[2],1))
r3 = roots(c(1, c3[3], 0,0,0, c3[3],1))
x = c(r1, r2, r3)
poly.calc(x)
1 - x - 2*x^2 + x^3 - x^5 - x^6 + x^7 - 2*x^8 - 2*x^10 +  
+ x^11 - x^12 - x^13 + x^15 - 2*x^16 - x^17 + x^18


### cos(2*pi/7)
r1 = roots(c(1, c3[1]^2 - c3[1], 0,0,0, c3[1]^2 - c3[1],1))
r2 = roots(c(1, c3[2]^2 - c3[2], 0,0,0, c3[2]^2 - c3[2],1))
r3 = roots(c(1, c3[3]^2 - c3[3], 0,0,0, c3[3]^2 - c3[3],1))
x = c(r1, r2, r3)
poly.calc(x)
1 + 6*x + 5*x^2 + x^3 + 6*x^5 + 13*x^6 + 15*x^7 + 5*x^8 + 5*x^10 +
+ 15*x^11 + 13*x^12 + 6*x^13 + x^15 + 5*x^16 + 6*x^17 + x^18


### (...)^(1/3)
k = 2^(1/3) * c(1, m3, m3^2)
r1 = roots(c(1, k[1], 0,0,0, k[1],1))
r2 = roots(c(1, k[2], 0,0,0, k[2],1))
r3 = roots(c(1, k[3], 0,0,0, k[3],1))
x = c(r1, r2, r3)
poly.calc(x)
1 + 2*x^3 + 3*x^6 + 6*x^7 + 6*x^11 + 3*x^12 + 2*x^15 + x^18
#
r1 = roots(c(1, k[1] * 1i, 0,0,0, k[1] * 1i,1)) / 1i
r2 = roots(c(1, k[2] * 1i, 0,0,0, k[2] * 1i,1)) / 1i
r3 = roots(c(1, k[3] * 1i, 0,0,0, k[3] * 1i,1)) / 1i
x = c(r1, r2, r3)
poly.calc(x)
-1 + 2*x^3 + 3*x^6 + 6*x^7 + 6*x^11 - 3*x^12 + 2*x^15 + x^18


### 2x entangled (...)^(1/3)
m3.grid = expand.grid(c(1, m3, m3^2), c(1, m3, m3^2))
k = 3^(1/3) * m3.grid[,1] - 2^(1/3) * m3.grid[,2]
r = sapply(k, function(k) roots(c(k, 1, 0,0,0, 1, k)) )
x = sort(r)
poly.calc(x)
1 + 165*x^3 + 12*x^6 + 495*x^7 + 991*x^9 + 18*x^10 + 495*x^11 + 45*x^12 + 2979*x^13 + 45*x^14 +
+ 2640*x^15 + 54*x^16 + 3006*x^17 + 153*x^18 + 7425*x^19 + 135*x^20 + 4374*x^21 + 99*x^22 + 7425*x^23 +
+ 309*x^24 + 10026*x^25 + 153*x^26 + 4950*x^27 + 153*x^28 + 10026*x^29 + 309*x^30 + 7425*x^31 +
+ 99*x^32 + 4374*x^33 + 135*x^34 + 7425*x^35 + 153*x^36 + 3006*x^37 + 54*x^38 + 2640*x^39 + 45*x^40 +
+ 2979*x^41 + 45*x^42 + 495*x^43 + 18*x^44 + 991*x^45 + 495*x^47 + 12*x^48 + 165*x^51 + x^54
# predict(p, x[54]) # NOT better!
#
r = sapply(k, function(k) roots(c(1, 1/k, 0,0,0, 1/k, 1)) )
x = sort(r)
poly.calc(x)
# the same polynomial


###############

m7 = unity(7, all=F)
m7 = m7^(1:6) # without 1!

###
x = sapply( m7, function(m) roots(c(1, 1-m, m)) )
round0.p(poly.calc(x))
1 - 7*x + 20*x^2 - 28*x^3 + 15*x^4 + 7*x^5 - 8*x^6 - 7*x^7 + 15*x^8 + 28*x^9 + 20*x^10 + 7*x^11 + x^12

### the same as with -m^1
x = sapply( m7, function(m) roots(c(1, 1-m^2, -m^2)) )
round0.p(poly.calc(x))
1 + 7*x + 22*x^2 + 42*x^3 + 57*x^4 + 63*x^5 + 64*x^6 + 63*x^7 + 57*x^8 + 42*x^9 + 22*x^10 + 7*x^11 + x^12
#
x = sapply( m7, function(m) roots(c(1, 1-m^2, m^2)) )
round0.p(poly.calc(x))
1 - 7*x + 20*x^2 - 28*x^3 + 15*x^4 + 7*x^5 - 8*x^6 - 7*x^7 + 15*x^8 + 28*x^9 + 20*x^10 + 7*x^11 + x^12


#################
