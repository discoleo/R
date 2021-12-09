
### Leonard Mada
###
### P3 Polynomial Systems
### Solver: Exact solutions
###
### draft 0.4-pre-alpha4

### P3 Systems
# v.0.4-pre-alpha3 - v.0.4-pre-alpha4:
# - special cases:
#  -- partly asymmetric order 2 with exact solution;
#  -- partly asymmetric order 3 with exact solution;
# - TODO: a lot of cleanup;
# v.0.4-pre-alpha2:
# - partly asymmetric order 2 S3;
# v.0.4-pre-alpha:
# - basic ideas to solve Asymmetric higher order systems;
# - basic example for fully asymmetric order 2 S3;
# v.0.3e: a less simplified version (see v.0.3b);
# v.0.3d:
# - partial extension of asymmetric system to 4 variables;
# - TODO: need the v.0.3.a solver as proper function!
# v.0.3c: added minimally asymmetric P[2];
# v.0.3b:
# - added also a greatly simplified version
#   of the basic asymmetric system;
#   [but with exact solution]
# v.0.3a: basic Asymmetric system;
# v.0.2d: more roots + classical "solution" to simple P3S3 (the P[9] polynomial);
# v.0.2c: Test the Linear decomposition concept;
# v.0.2a: S3 system + linear (x+y+z) terms;
# v.0.1: simple S3 system: the Base System;


######################

### Theory

### Order n = O(n) Simple S3 System
# B1 * (x^n, y^n, z^n) = R1
# B2 * (x*y, x*z, y*z) = R2
# x*y*z = R3

# 1.) Symmetric O(n):
# - B1 = B2 = (1, 1, 1);
# - can be decomposed into 2 entangled systems:
# - solve a P[n]-polynomial: with exact solution for n <= 3;
# - solve a derived basic S3: exact solution;
# 2.) Partly Asymmetric O(n)
# - b12 = b21^n, b13 = b22^n;
# - solve a P[3*n] polynomial:
#  -- exact solution for n=1 or for other special cases;
# 3.) Fully Asymmetric O(n)
# - solve a P[6*n] polynomial;
# - a P[12] example is shown below;

### Order n = O(n) Complex P3 System
# see below for Theory of complex O(3) P3;


#######################

### A.) Base S3 system

### S3 System:
# x^3 + y^3 + z^3 = R1
# x*y + x*z + y*z = R2
# x*y*z = R3


### Solution Steps:

### Subsystem 1:
#
# S^3 - 3*R2*S + 3*R3 - R1 = 0
# solve for S;
# where S = x + y + z;

### Subsystem 2:
#
# x + y + z = S
# x*y + x*z + y*z = R2
# x*y*z = R3
# - solve for x, y, z;
# - exact solution provided;


### B.) Base System: Classical Solution
# - the classical solution to the base-system [A];
# - involves P[9]:
# x*y*z = C => y*z = C/x; # C = P (in part [A])
# x*y + x*z + y*z = B => x*(y+z) = B - C/x
# => y+z = B/x - C/x^2;
# =>
# y^3 + z^3 = (y+z)*((y+z)^2 - 3*yz)
# = (B/x - C/x^2) * ((B/x - C/x^2)^2 - 3*C/x)
# =>
# x^9 - A*x^6 - 3*B*C*x^4 + (B^3+3*C^2)*x^3 - 3*B^2*C*x^2 + 3*B*C^2*x - C^3 # = 0

# Test
x = sqrt(c(2,3,5))
A = sum(x^3)
B = sum(x[1]*x[2], x[1]*x[3], x[2]*x[3])
C = prod(x)
x^9 - A*x^6 - 3*B*C*x^4 + (B^3+3*C^2)*x^3 -3*B^2*C*x^2 + 3*B*C^2*x - C^3


### C.) Complex S3 System
# b31*(x^3+y^3+z^3) + b21*(x^2+y^2+z^2) + b11*(x+y+z) + e21*(xy+xz+yz) + e31*x*y*z = A1
# b32*(x^3+y^3+z^3) + b22*(x^2+y^2+z^2) + b12*(x+y+z) + e22*(xy+xz+yz) + e32*x*y*z = A2
# b33*(x^3+y^3+z^3) + b23*(x^2+y^2+z^2) + b13*(x+y+z) + e23*(xy+xz+yz) + e33*x*y*z = A3

# see below for complet *exact* solution!
# 1.) the system is decomposed using a liniar decomposition;
# 2.) then the "base"-terms are computed;
# 3.) then the S3 system is solved based on [A];


### D.) Perturbations to break Symmetry
# [C] is still symetrical;
# see below for some ideas;



### Asymmetrical Simple Systems

### E.) Basic Assymetrical System:
# b11 * x + b12 * y + b13 * z = R1
# b21 * x*y + b22 * x*z + b23 * y*z = R2
# x*y*z = R3
# see part [E] for solution;
# [but invovles numerical solution to a P[6] polynomial;]

### E.2.) Simplified versions of [E]
# - but with exact solution;
# - a greatly simplified;
# - and the "classic" partly simplified version;


####################

### helper functions

# dir.create(tempdir())

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
### complex/matrix round
round0 = function(m, tol=1E-10) {
	isZero = abs(Im(m)) < tol
	m[isZero] = Re(m[isZero])
	m[isZero & abs(Re(m)) < tol] = 0
	
	isZero = ( ! isZero ) & (abs(Re(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	return(m)
}

### Solution
# Cubic Polynomial
solveP3 <- function(c, d, n=3, all=TRUE) {
	det = (d^2 - c^n)
	det = sqrt.c(det)
	
	p = nroot.c(d + det, n)
	q = nroot.c(d - det, n)
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
# Cubic System
solveP3Outer = function(s, bi, pr, n=3) {
	if(n == 1) {
		d2 = s
	} else if(n == 2) {
		d2 = sqrt(s + 2*bi)
	} else if(n == 3) {
		# y^3 - 3*bi * y + 3*pr - s
		d = -(3*pr - s) / 2
		d2 = solveP3(bi, d)
		# Test
		y = d2
		err = y^3 - 3*bi * y + 3*pr - s
		cat("err Tri = ")
		print(round0(err))
	} else {
		print("Not yet implemented!")
		break
	}
	return(d2)
}
solveP3OuterShifted = function(B, b.shifts, n=3) {
	# X^3 + 3*b21*X^2 - (3*BB + 3*b11 - b31)*X + 3*P - A = 0
	A  = B[1]
	BB = B[2]
	P  = B[3]
	coeff = c(3*b.shifts[2], - 3*(BB + b.shifts[3] - b.shifts[1]/3), 3*P - A)
	if(length(b.shifts) / n == 2) {
		# TODO: X^2
	}
	b = shift.poly(coeff)
	b.b = b$b
	r = solveP3( - b.b[2]/3, - b.b[3]/2)
	r = r + b$shift
	return(r)
}
solveP3S.Base = function(s, bi, pr, d2, n=3, b.shifts=c(0,0,0)) {
	### Solve P3 System:
	# y^3 - 3*B*y^2 + pr*d2*y - pr^2
	# y = shifted by B = bi/3
	# y^3 - 3*B^2*y + pr*d2*y - 2*B^3 + pr*d2*B - pr^2
	B = bi / 3
	d = B^3 - pr*(d2*B - pr)/2
	y = solveP3(B^2 - pr*d2/3, d, all=TRUE)
	cat("\nx_shifted: ")
	print(paste("x =", round(y, 5)))
	y = y + B
	# Test
	err = y^3 - 3*B*y^2 + pr*d2*y - pr^2
	print(paste("err =", round0(err), collapse=", "))
	# x1, x2, x3
	x1 = pr/y
	syz = d2 - x1
	diff = sqrt.c(syz^2 - 4*y)
	x2 = (syz + diff)/2
	x3 = (syz - diff)/2
	p = matrix(c(x1, x2, x3), ncol=3)
	# Test
	m = 1
	test = c(
		apply(p, m, function(p) sum(p^n) + b.shifts[1]*sum(p)),
		apply(p, m, function(p) sum(p[1] * p[-1], p[2]*p[3]) + b.shifts[2]*sum(p)),
		apply(p, m, function(p) prod(p) + b.shifts[3]*sum(p)))
	return(list("x"=p, "test" = test))
}
t3.m = function(m, ncol=3) {
	matrix(m, ncol=ncol, byrow=T)
}
solveP3S = function(s, bi, pr, n=3) {
	d2 = solveP3Outer(s, bi, pr, n)
	l = lapply(d2, function(d2) solveP3S.Base(s, bi, pr, d2, n))
	# print(l) ### DEBUG
	l.x = sapply(l, function(l) t(l$x))
	l.test = sapply(l, function(l) t(l$test))
	return(list(x=t3.m(l.x), test=t3.m(l.test)))
}
solveP3SS = function(B, b, n=3) {
	# X^3 + 3*b21*X^2 - (3*BB + b31 + b11)*X + 3*P - A = 0
	r = solveP3OuterShifted(B, b, n)
	solve.p3 = function(r) {
		B.shifted = B - r * b
		solveP3S.Base(B.shifted[1], B.shifted[2], B.shifted[3], r, n, b.shifts=b)
	}
	l = lapply(r, solve.p3)
	l.x = sapply(l, function(l) t(l$x))
	l.test = sapply(l, function(l) t(l$test))
	return(list(x=t3.m(l.x), test=t3.m(l.test)))
}

############################

### Test

solveP3S(1,2,3)

# Parameters
A = 1
B = 2
C = -1
#
x.all = solveP3S(A, B, C)
x.all
# Test
x = x.all$x
err = x^9 - A*x^6 - 3*B*C*x^4 + (B^3+3*C^2)*x^3 -3*B^2*C*x^2 + 3*B*C^2*x - C^3
round0(err)


# TODO:
# - CORRECTED: transposition of the solutions was needed;
# - still check if all roots are computed correctly & add missing permutations;
# - known: tests are triplicated;
unique(round(x.all$x, 12))


#############################


### P3 System:
# x^3 + y^3 + z^3 + b31 * (x+y+z) = A
# x*y + x*z + y*z + b21 * (x+y+z) = B
# x*y*z + b11 * (x+y+z) = P

### Steps:

### Subsystem 1:
#
# X^3 + 3*b21*X^2 - (3*B + 3*b11 - b31)*X + 3*P - A = 0
# solve for X;
# then solve Subsystem 2;

### Test

B = c(1, 2, -1)
b = c(1, 2, 3)
#
x.all = solveP3SS(B, b)
x.all


B = c(1, 2, -1)
b = c(1, 2, 4)
#
x.all = solveP3SS(B, b)
x.all

###
x = x.all$x[1,]
# ...

###################
###
### "Linear"-System
###
###################

### Terms:
# E[i,j] = (x^3+y^3+z^3) & E[i] on row j (i=2:3);
# B[i, j] = currently only (x+y+z) on row j;
# - TODO: add the (x^2+y^2+z^2) terms;

### "ALL"-Solutions
# TODO: verify if correct!

### Test
E.m = matrix(
c(1, 2, 3,
  1, -1, 2,
 -1, 1, 4), ncol=3, byrow=TRUE)

B.m = matrix(
c(-2, -3, -5), ncol=1, byrow=TRUE)
V.m = matrix(
c(1, 2, -1), ncol=1)

# Solve "Linear"-System
B.r = solve(E.m, cbind(V.m, B.m))

x.all = solveP3SS(B.r[,1], B.r[,2])
x.all

# Test solution
x = x.all$x[2,]
X = sum(x)
XE = c(sum(x^3), sum(x[1]*x[2], x[1]*x[3], x[2]*x[3]), prod(x))
E.m %*% XE + B.m * X
V.m


##########################

### TODO:

### A.) Solve the basic types of systems

### Extension Path:
# A.1.) add X^2 term, where X = sum(x[i]);

# A.2.) solve "liniar" meta-system:
# [X is grouped with the power terms,
# NOT with the elementary polynomials,
# because of the way of the solution]
# X^2, X = "slack" variables;

# p31 * sum(x^3) + e31 * E3 + e21 * E2 + b21 * X^2 + b11 * X = v1
# p32 * sum(x^3) + e32 * E3 + e22 * E2 + b22 * X^2 + b12 * X = v2
# p33 * sum(x^3) + e33 * E3 + e23 * E2 + b23 * X^2 + b13 * X = v3

# [p, e] x [sum(x^3), E] = V - B x [X^2, X]

# A.3.) add sum(x[i]^2) terms instead of the X^2
# - after having the "linear" solution: it is easy to transform and solve;
#   sum(x[i]^2) = X^2 - 2*E2;


#######################
###
### Perturbation Theory

### D.) Perturbations to break Symmetry

# This will be a new branch in mathematics.
# It may be useful to study systematically such perturbations and
# understand the behaviour of the perturbed systems.


# B.1.) Class 1 Perturbations
# B.1.a.) x1 = x1' + s1;
# B.1.b.) x1 = a1*x1' + s1;
# [+ applied on multiple variables]

# B.2.) Class 2 Perturbations
# B.2.1.) Symmetrical:
#         x1 = x1' + x2'; x2 = x2' + x3'; x3 = x3' + x1';
# B.2.2.) x1 = a1*x1' + a2*x2;
# [+ applied on multiple variables]

# B.3.) Class 3 Perturbations: Powers
# when it affects only 1 of the variables:
# so the remaining have consistent exponents lower than the target variable;

# B.3..a.) x1 = x1'^2;
# B.3.b.) x1 = x1'^2 + a2*x2;
# It is easier to solve an order 3 system, than an order 6 system.


########################

########################
### Assymetrical Systems
### E.) Basic System

# (b11, b12, b13) * (x, y, z) = R1
# (b21, b22, b23) * (xy, xz, yz) = R2
# x*y*z = R3

# => y*z = R3 / x
# (b12, b13) * (y, z) = - b11 * x + R1
# (b21, b22) * (y, z) = R2/x - b23 * R3/x^2
# => linear solution for y & z;


library(polynom)
# needed for poynomials with complex coefficients
library(pracma)

solve.P3asym = function(b1, b2, R) {
	# needs package "pracma"
	# Note:
	# - package polynom cannot handle polynomials with complex coefficients;
	# - such polynomials do come up in section [F];
	
	# Solve linear sub-system
	b.yz = rbind(b1[-1], b2[-3])
	# R.r2 = (x, b0, 1/x, 1/x^2)
	R.r2 = rbind( c( -b1[1], R[1], 0, 0), c(0, 0, R[2], -R[3] * b2[3]) )
	yz = solve(b.yz, R.r2)
	# (b11, b12, b13) * (x, y, z) = R1
	# should be == 0;
	# b.coeff = c(b1[1], -R[1], 0, 0) +  b1[2] * yz[1,] + b1[3] * yz[2,]
	# print(b.coeff)

	# Solve polynomial
	# TODO: cleanup, improve;
	# print(yz)
	p1 = polynomial.c(rev(yz[1,]))
	p2 = polynomial.c(rev(yz[2,]))
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	p = p - c(0,0,0,R[3],0,0,0)
	# p = polynomial.c(p)
	print(p)
	# x = solve(p)
	x = roots(rev(p))
	y = sapply(x, function(x) sum(yz[1,] * c(x, 1, 1/x, 1/x^2)) )
	z = sapply(x, function(x) sum(yz[2,] * c(x, 1, 1/x, 1/x^2)) )

	### Solution
	sol = rbind(x, y, z)
	sol
}

### free Parameters
R = c(1, 1, 1)
b1 = c(1, 2, 3)
b2 = c(1, 3, 4)
# Solve linear sub-system
b.yz = rbind(b1[-1], b2[-3])
# R.r2 = (x, b0, 1/x, 1/x^2)
R.r2 = rbind( c( -b1[1], R[1], 0, 0), c(0, 0, R[2], -R[3] * b2[3]) )
yz = solve(b.yz, R.r2)
# (b11, b12, b13) * (x, y, z) = R1
# should be == 0;
b.coeff = c(b1[1], -R[1], 0, 0) +  b1[2] * yz[1,] + b1[3] * yz[2,]
b.coeff

# Solve polynomial
p = polynomial(rev(yz[1,])) * polynomial(rev(yz[2,])) - polynomial(c(0,0,0,R[3]))
p
x = solve(p)
y = sapply(x, function(x) sum(yz[1,] * c(x, 1, 1/x, 1/x^2)) )
z = sapply(x, function(x) sum(yz[2,] * c(x, 1, 1/x, 1/x^2)) )

### Solution
sol = rbind(x, y, z)
sol
### Test
b1 %*% sol - R[1]
#
err = sapply(1:ncol(sol), function(id) sum(b2 * c(sol[1,id]*sol[2,id], sol[1,id]*sol[3,id], sol[2,id]*sol[3,id]))) - R[2]
round0(err)
#
err = sapply(1:ncol(sol), function(id) prod(sol[,id])) - R[3]
round0(err)

###########################
### E.2.) P3 variant System
###   Greatly Simplified
###   (partly) Assymetric

### Exact Solution
# b11*x + b1*(y+z) = R1
# b2*(x*y + x*z) + b23*y*z = R2
# x*y*z = R3

# => y*z = R3 / x
# => y + z = R1/b1 - b11/b1 * x
# b2*x*(R1/b1 - b11/b1 * x) + b23*R3/x = R2
# - b2*b11/b1*x^3 + b2/b1*R1*x^2 - R2*x + b23*R3 = 0
# x^3 - R1/b11 * x^2 + b1*R2/(b2*b11) * x - b1*b23*R3/(b2*b11) = 0

### free Parameters:
b1 = c(1, 3) # b11, b12_b13
b2 = c(2, 5) # b21_b22, b23
R = c(1,1,1)
### Solution
b3 = - b2[1]*b1[1]/b1[2]
b = c(b2[1]/b1[2] * R[1], - R[2], b2[2]*R[3])
b = b / b3
b
b.shifted = shift.poly(b)
b.shifted
# x
x = solveP3(b.shifted$b[2]/-3, b.shifted$b[3]/-2) + b.shifted$shift
x
# y, z
yz = R[3] / x
yz.s = R[1]/b1[2] - b1[1]/b1[2] * x
yz.minus = sqrt(yz.s^2 - 4*yz)
y = (yz.s + yz.minus)/2
z = (yz.s - yz.minus)/2
# complete solution:
sol = cbind(x, y, z)
sol = rbind(sol, cbind(x, z, y))
sol

### Test
sapply(1:nrow(sol), function(id) sum(b1[c(1,2,2)]*sol[id,]))
sapply(1:nrow(sol), function(id) sum(b2[c(1,1,2)]*sol[id,c(2,1,3)]*sol[id,c(1,3,2)]))
sapply(1:nrow(sol), function(id) prod(sol[id,]))


############################
### E.2.b) P3 variant System
### Partly Assymetric
### Less Simplified than E.2.a.

### Exact Solution
# b11*x + b12*y + b13*z = R1
# b2c * x*(b12*y + b13*z) + b23*y*z = R2
# x*y*z = R3

# => y*z = R3 / x
# => b12*y + b13*z = R1 - b11 * x
# b2c*x*(R1 - b11 * x) + b23*R3/x = R2
# - b2c*b11*x^3 + b2c*R1*x^2 - R2*x + b23*R3 = 0

solve.yz = function(yz, yz.s, b1, sign=+1) {
	yz.minus = sign * sqrt(yz.s^2 - 4*b1[2]*b1[3]*yz)
	y = (yz.s + yz.minus)/2 / b1[2]
	z = (yz.s - yz.minus)/2 / b1[3]
	cbind(y, z)
}

### free Parameters:
b1 = c(1, 2, 3)
b2c = 2
b2 = c(b1[-1], 5)
R = c(1,1,1)
### Solution
b3 = - b2c*b1[1]
b = c(b2c * R[1], - R[2], b2[3]*R[3])
b = b / b3
b
b.shifted = shift.poly(b)
b.shifted
# x
x = solveP3(b.shifted$b[2]/-3, b.shifted$b[3]/-2) + b.shifted$shift
x
# y, z
yz = R[3] / x
yz.s = R[1] - b1[1] * x # b12*y + b13*z
sol  = solve.yz(yz, yz.s, b1)
sol2 = solve.yz(yz, yz.s, b1, sign=-1)
# complete solution:
sol = cbind(x, sol)
sol = rbind(sol, cbind(x, sol2))
sol

### Test
sapply(1:nrow(sol), function(id) sum(b1*sol[id,]))
sapply(1:nrow(sol), function(id) sum(c(b2c,b2c,1) * b2 * sol[id,c(2,1,3)]*sol[id,c(1,3,2)]))
sapply(1:nrow(sol), function(id) prod(sol[id,]))


############################

############################
### Asymmetrical Systems ###
############################

#########################
### F.) Higher Orders ###
#########################


#########################
### Partly Asymmetric ###
#########################

# b11*x^2 + b21^2*y^2 + b22^2*z^2 = R1
# x*(b21*y + b22*z) + b23*y*z = R2
# x*y*z = R3

# => y*z = R3/x
# => b21*y + b22*z = R2/x - b23*R3/x^2
#
# b11*x^2 + (b21*y + b22*z)^2 - 2*b21*b22*R3/x = R1
# b11*x^6 + (R2*x - b23*R3)^2 - 2*b21*b22*R3*x^3 = R1*x^4
b11*x^6 - R1*x^4 - 2*b21*b22*R3*x^3 + R2^2*x^2 - 2*b23*R2*R3*x + b23^2*R3^2 # = 0


library(polynom)

b11 = 2
b2  = c(1, 2, 3)
R   = c(1,1,1)
### Solution:
b = c(- R[1], - 2*b2[1]*b2[2]*R[3], R[2]^2, - 2*b2[3]*R[2]*R[3], b2[3]^2*R[3]^2) / b11
b = rev(c(1, 0, b))

p = polynomial(b)
p
x = solve(p)
#
yz = R[3] / x
yz.p = (R[2] - b2[3]*yz)/x
yz.m = sqrt(yz.p^2 - 4*b2[1]*b2[2]*yz)
y = (yz.p + yz.m) / 2 / b2[1]
z = (yz.p - yz.m) / 2 / b2[2]
sol = cbind(x, y, z)
sol

# TODO: all solutions

### Test
sapply(1:nrow(sol), function(id) sum(sol[id,]^2 * c(b11, b2[-3]^2)))
sapply(1:nrow(sol), function(id) sum(b2*sol[id,c(1,1,2)]*sol[id,c(2,3,3)]))
sapply(1:nrow(sol), function(id) prod(sol[id,]))


### Special Cases:
# R2 = 0
b11*x^6 - R1*x^4 - 2*b21*b22*R3*x^3 + b23^2*R3^2 # = 0
# b11 = 1; b23 = {+/-} b21*b22;
x^6 - R1*x^4 - 2*b21*b22*R3*x^3 + (b21*b22)^2*R3^2 # = 0
(x^3 - b21*b22*R3)^2 - R1*x^4 # = 0

library(pracma)

solve.psym.S3P2 = function(R, b) {
	if(R[2] != 0) stop("Can solve only when R[2] == 0!")
	# assumes b23 = b[1]*b[2]!
	x = roots(c(1, 0, -R[1], - 2*b[1]*b[2]*R[3], 0, 0, (b[1]*b[2])^2*R[3]^2))
	yz = R[3] / x;
	yz.s = - b[1]*b[2]*R[3]/x^2;
	yz.d = sqrt(yz.s^2 - 4*b[1]*b[2]*yz + 0i);
	y = (yz.s + yz.d)/2 / b[1];
	z = (yz.s - yz.d)/2 / b[2];
	# x = c(x, x); tmp = y; y = c(y, z); z = c(z, tmp);
	return(cbind(x=x, y=y, z=z))
}

###
R = c(1, 0, 2)
b = c(-1, 2)
#
sol = solve.psym.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + b[1]^2*y^2 + b[2]^2*z^2 # - R1
x*(b[1]*y + b[2]*z) + b[1]*b[2]*y*z # = 0
x*y*z # - R3



#########################
### Partly Asymmetric
### [very simple version]

# x^2 + y^2 + z^2 = R1
# b2*x*(y+z) + b23*y*z = R2
# x*y*z = R3

# Note:
# - the very simple version does NOT simplify additionally the system;

### Solution
# => y*z = R3/x
# => x*(y+z) = R2 - b23*R3/x
# x^2*(y+z)^2 = R2^2 - 2*b23*R2*R3/x + b23^2*R3^2/x^2
# x^2*(y^2 + z^2) = -2*R3*x + R2^2 - 2*b23*R2*R3/x + b23^2*R3^2/x^2
# =>
# x^4 - 2*R3*x + R2^2 - 2*b23*R2*R3/x + b23^2*R3^2/x^2 = R1*x^2
# x^6 - R1*x^4 - 2*R3*x^3 + R2^2*x^2 - 2*b23*R2*R3*x + b23^2*R3^2 = 0


library(polynom)

### free Parameters:
# [reusing code from TRUE Partly Asymmetric]
b11 = 1
b2  = c(1, 1, 3)
R   = c(1, 1, 1)
### Solution:
b = c(- R[1], - 2*b2[1]*b2[2]*R[3], R[2]^2, - 2*b2[3]*R[2]*R[3], b2[3]^2*R[3]^2) / b11
b = rev(c(1, 0, b))

p = polynomial(b)
p
x = solve(p)
#
yz = R[3] / x
yz.p = (R[2] - b2[3]*yz)/x
yz.m = sqrt(yz.p^2 - 4*b2[1]*b2[2]*yz)
y = (yz.p + yz.m) / 2 / b2[1]
z = (yz.p - yz.m) / 2 / b2[2]
sol = cbind(x, y, z)
sol

# TODO: all solutions

### Test
sapply(1:nrow(sol), function(id) sum(sol[id,]^2 * c(b11, b2[-3]^2)))
sapply(1:nrow(sol), function(id) sum(b2*sol[id,c(1,1,2)]*sol[id,c(2,3,3)]))
sapply(1:nrow(sol), function(id) prod(sol[id,]))


#########################


#########################
### Completely Assymetric
###
### using Cardan-type Polynomials

# (b21*y + b22*z)^n = b21^n*y^n + b22^n*z^n + f((b21*y + b22*z), y*z);

### P3 Order 2:
# b11*x^2 + b12*y^2 + b13*z^2 = R1
# b21*x*y + b22*x*z + b23*y*z = R2
# x*y*z = R3

# => y*z = R3/x
# => b21*y + b22*z = R2/x - b23*R3/x^2
# =>
# b12*y^2   + b13*z^2   = - b11*x^2 + R1
# b21^2*y^2 + b22^2*z^2 = (R2/x - b23*R3/x^2)^2 - 2*b21*b22*R3/x
# b21^2*y^2 + b22^2*z^2 = - 2*b21*b22*R3/x + R2^2/x^2 - 2*b23*R2*R3/x^3 + b23^2*R3^2/x^4
# =>
# solve linear system: => y^2, z^2;
# => substitute y^2 & z^2:
# x^2*y^2*z^2 = R3^2
# x^8*y^2*z^2 - R3^2 * x^6 = 0

library(polynom)

### free Parameters
b1 = c(1,2,3)
b2 = c(1,1,2)
R  = c(1,1,1)
### linear solution for sub-system:
yz.m = matrix(c(b1[-1], b2[-3]^2), ncol=2, byrow=TRUE)
yz.R = matrix(
	c(-b1[1], 0, R[1], 0, 0, 0, 0,
	c(0, 0, 0, - 2*b2[1]*b2[2]*R[3], R[2]^2, - 2*b2[3]*R[2]*R[3], b2[3]^2*R[3]^2)), ncol=7, byrow=TRUE)

yz.coeff = solve(yz.m, yz.R)
yz.R
yz.coeff

p = polynomial(rev(yz.coeff[1,])) * polynomial(rev(yz.coeff[2,])) - polynomial(c(rep(0,6),R[3]^2))
p
x = solve(p)
x

# solve y, z:
pow = 2 - (0:6)
solve.yz = function(x, yz.coeff, id, pow = 2 - (0:6)) sqrt(sum(yz.coeff[id,] * x^pow))
y = sapply(x, solve.yz, yz.coeff, 1)
z = sapply(x, solve.yz, yz.coeff, 2)
# all solutions
sol = cbind(x, y, z)
sol = rbind(sol, cbind(x, -y, z), 
	cbind(x, y, -z), cbind(x, -y, -z))

# TODO: solve y & z using different approach;

### Test
err1 = apply(sol, 1, function(sol) sum(b1*sol^2))
err2 = apply(sol, 1, function(sol) sum(b2*sol[c(1,1,2)]*sol[c(2,3,3)]))
err3 = apply(sol, 1, function(sol) prod(sol))
err = cbind(err1, err2, err3)
err = round0(err)
correct = (round0(err[,2] - R[2]) == 0) & (round0(err[3] - R[3]) == 0)
err[correct,]
sol[correct,]


###############
### P3 Order 3:
# b11*x^3 + b12*y^3 + b13*z^3 = R1
# b21*x*y + b22*x*z + b23*y*z = R2
# x*y*z = R3

# => y*z = R3/x
# => b21*y + b22*z = R2/x - b23*R3/x^2
# =>
# b12*y^3   + b13*z^3   = R1 - b11*x^3
# b21^3*y^3 + b22^3*z^3 = (R2/x - b23*R3/x^2)^3 - 3*b21*b22*R3/x * (R2/x - b23*R3/x^2)
# =>
# solve linear system: => y^3, z^3;
# => substitute y^3 & z^3:
# x^3*y^3*z^3 = R3^3




#####################
#####################

#####################
### Other Systems ###

# b11*x^2 + b12*y^2 + b13*z^2 = R1
# b21*x*y + b22*x*z + b23*y*z = R2
# x^2*y*z = R3

# => y*z = R3 / x^2
# => b21*y + b22*z = R2/x - b23*R3/x^3
# b21^2*y^2 + b22^2*z^2 = (R2^2 - 2*b21*b22*R3)/x^2 - b23*R2*R3/x^4 + b23^2*R3^2/x^6
# b12*y^2 + b13*z^2 = - b11*x^2 + R1


###
# b11*x^4 + b12*y^4 + b13*z^4 = R1
# b21*x*y + b22*x*z + b23*y*z = R2
# x^4*y*z = R3  # TODO: x^2*y*z = R3

# => y*z = R3 / x^4
# => b21*y + b22*z = R2/x - b23*R3/x^5
# b21^4*y^4 + b22^4*z^4 = (R2 - b23*R3/x^4)/x^4 - 4*b21*b22*R3/x^4*(...) - 6*b21^2*b22^2*R3^2/x^8
# b12*y^4 + b13*z^4 = - b11*x^4 + R1



####################
### Partly Symmetric
# b11*x^3 + (b12*y + b13*z)^3 = R1
# b12*x*y + b13*x*z + b23*y*z = R2
# x^3*y*z = R3

# b11*x^6 + (R2 - b23*R/x^3)^3 = R1*x^3
# b11*x^15 + (R2*x^3 - b23*R3)^3 = R1*x^12


### Special Case:
# R2 = 0
x^9 - R1*x^6 + 3*b1*b2*b3*R3^2*x^3 - b3^3*R3^3 # = 0

### Example:
R = c(2, 0, 1)
b = c(2, 3, 1)
#
x = roots(c(1,0,0,-R[1],0,0,3*b[1]*b[2]*b[3]*R[3]^2,0,0,- b[3]^3*R[3]^3))
yz = R[3]/x;
yz.s = - b[3] * yz / x;
yz.d = sqrt(yz.s^2 - 4*b[1]*b[2]*yz + 0i);
y = (yz.s + yz.d) / 2 / b[1];
z = (yz.s - yz.d) / 2 / b[2];

### Test
x^3 + b[1]^3*y^3 + b[2]^3*z^3
b[1]*x*y + b[2]*x*z + b[3]*y*z
x*y*z


###################
### 
# b11*x^3 + y^3 + z^3 = R1
# y + z = R2
# x*y*z = R3

# b11*x^3 + R2^3 - 3*R2*R3/x = R1
# b11*x^4 + (R2^3-R1)*x - 3*R2*R3 = 0


### 
# b11*x^3 + y^3 + z^3 = R1
# x*(y + z) = R2
# x^2*y*z = R3

# b11*x^3 + R2^3/x^3 - 3*R2*R3/x^3 = R1
# b11*x^6 - R1*x^3 + R2^3 - 3*R2*R3 = 0


###################
###################


###################
### X.) 4 Variables
### Partly Assymetric

# b11 * x1 + b1*(x2, x3, x4) = R1
# b1*x1*(x2, x3, x4) + b22*(x2*x3, x2*x4, x3*x4) = R2
# b22*x1*(x2*x3, x2*x4, x3*x4) + b32*x2*x3*x4 = R3
# x1*x2*x3*x4 = R4

# => x2*x3*x4 = R4/x1
# => b22*(x2*x3, x2*x4, x3*x4) = R3 / x1 - b32*R4 / x1^2
# => b1*(x2, x3, x4) = R2 / x1 - R3 / x1^2 + b32*R4 / x1^3
# =>
# b11*x1 + R2 / x1 - R3 / x1^2 + b32*R4 / x1^3 = R1
# b11*x^4 - R1*x^3 + R2*x^2 - R3*x + b32*R4 = 0


library(polynom)
library(pracma)


####################

### helper functions
polynomial.c = function(coef) {
    a <- coef
    while ((la <- length(a)) > 1 && a[la] == 0) a <- a[-la]
    structure(a, class = "polynomial")
}
root1.f = function(n=5, keep=TRUE, positive=TRUE) {
	# keep = include 1 (or -1 for odd roots of -1);
	i.pos = ifelse(positive, 2, 1)
	m = complex(re=cos(i.pos * pi/n), im=sin(i.pos * pi/n))
	if(positive) {
		m = m^(1:(n-1))
		if(keep) {
			m = c(1, m)
		}
	} else {
		if(n %% 2 == 1) {
			i.max = (n - 3) %/% 2
		} else {
			i.max = n %/% 2 - 1
		}
		m = m^(2* 0:i.max + 1)
		m = c(m, 1/m)
		if(keep && n %% 2 == 1) {
			m = c(-1, m)
		}
	}
	return(m)
}
# function to compute coefficients of polynomial
elemPoly = function(x, start=0, adjustSign=TRUE, tol=1E-7) {
	coeff = sapply(c(start, seq_along(x)), function(n) {
		if(n > 10) {
			print(paste("Iteration:", n))
			flush.console() # to display in real time
		}
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
		# coeff[abs(coeff) < 1E-10 ] = 0
		coeff = round0(coeff, tol=tol)
		coeff = coeff * adj
	}
	isComplex = (Im(coeff) != 0)
	if(length(isComplex[isComplex]) == 0) {
		coeff = as.numeric(coeff)
	}
	
	n = len - 1
	n_1 = n - 1
	poly.coeff = round.complex(coeff)
	#
	poly.str = toPoly(poly.coeff)
	
	poly.list = list(
		r = x,
		poly.coeff = coeff,
		poly = poly.str,
		n = n
	)
	return(poly.list)
}
toPoly = function(poly.coeff, strVar="x") {
	len = length(poly.coeff)
	n   = len - 1
	n_1 = n - 1
	coeff.sign = rep("+", len)
	coeff.sign[1] = ""
	coeff.sign[ Re(poly.coeff) < 0] = "-"
	poly.coeff.pos = poly.coeff
	### TODO: Im() != 0
	isComplex = (Re(poly.coeff) < 0) & (Im(poly.coeff) != 0)
	isNegativ = (Re(poly.coeff) < 0) & (Im(poly.coeff) == 0)
	poly.coeff.pos[isNegativ] = -1 * poly.coeff.pos[isNegativ]
	if(length(isComplex[isComplex]) > 0) {
		poly.coeff.pos[isComplex] = complex(
			re = - Re(poly.coeff.pos[isComplex]),
			im = Im(poly.coeff.pos[isComplex])  )
	}
	coeff.str = as.character( poly.coeff.pos)
	isCoeff = (poly.coeff.pos[-len] != 1)
	coeff.str[-len][isCoeff] = paste(coeff.str[-len][isCoeff], "*", strVar, sep="")
	coeff.str[-len][ ! isCoeff ] = strVar
	coeff.str = paste(coeff.sign, coeff.str, sep=" ")
	
	poly.str = paste(coeff.str[1:n_1], "^", n:2, sep="")
	poly.str = c(poly.str, coeff.str[n:len])
	poly.str = poly.str[poly.coeff != 0]
	return(paste(poly.str, collapse=" "))
}
round.complex = function(m, digits=0) {
	isComplex = (Im(m) != 0)
	if(length(isComplex[isComplex]) > 0) {
		m[isComplex] = complex(
			re=round(Re(m[isComplex]), digits),
			im=round(Im(m[isComplex]), digits))
	}
	m[ ! isComplex] = round(Re(m[ ! isComplex]), digits)
	return(m)
}

####################

### free Parameters
b1 = c(1, 2,3,4)
b2 = c(b1[-1], 3,3,1)
b3 = c(b2[4:6], 2)
R = c(1,1,1,1)
#
b.coeff = c(b1[1], - R[1], R[2], - R[3], b3[4]*R[4])
p = polynomial(rev(b.coeff))

x = solve(p)
x

### TODO:
# - solve assymetric S3 system;
# - implement solver (see previous sections) as function;

solve.S3subsys = function(x) {
	R1.new = R[1] - b1[1]*x
	R.new = c(R1.new, R[2] - b1[1]*x*R1.new, R[4]/x)
	sol = solve.P3asym(b1[2:4], b2[4:6], R.new)
	matrix(t(sol), ncol=3)
}

sol = sapply(x, solve.S3subsys)
sol = solve.S3subsys(x[1])
sol = cbind(x[1], sol)
colnames(sol) = c(paste("x", 1:4, sep=""))
sol

### Test
apply(sol, 1, function(x) sum(b1*x))
apply(sol, 1, function(x) sum(b2*c(x[1]*x[-1], x[2]*x[-(1:2)], x[3]*x[4])))
apply(sol, 1, function(x) sum(b3*c(x[1]*x[2]*x[3], x[1]*x[2]*x[4], x[1]*x[3]*x[4], x[2]*x[3]*x[4])))
apply(sol, 1, function(x) prod(x))

### Debug
a = c(2,3,1,3)
R = c(sum(b1*a),
	sum(b2*c(a[1]*a[-1], a[2]*a[-(1:2)], a[3]*a[4])),
	sum(b3*c(a[1]*a[2]*a[3], a[1]*a[2]*a[4], a[1]*a[3]*a[4], a[2]*a[3]*a[4])),
	prod(a))
sum((a[2]*a[3])^(6:0) * b.coeff)

