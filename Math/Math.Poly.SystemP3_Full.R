
### Leonard Mada
###
### P3 Polynomial Systems
### Solver: Exact solutions
###
### draft 0.3b

### P3 Systems
# v.0.3.b:
#  - added also a greatly simplified version
#    of the basic assymetric system;
#    [but with exact solution]
# v.0.3a: basic Assymetric system;
# v.0.2d: more roots + classical "solution" to simple PS3 (the P[9] polynomial);
# v.0.2c: Test the Linear decomposition concept;
# v.0.2a: P3 system + linear (x+y+z) terms;
# v.0.1: simple P3 system: the Base System;

#####################

### A.) Base P3 system

### P3 System:
# x^3 + y^3 + z^3 = A
# x*y + x*z + y*z = B
# x*y*z = P


### Solution Steps:

### Subsystem 1:
#
# X^3 - 3*B*X + 3*P - A = 0
# solve for X;
# where X = x + y + z;

### Subsystem 2:
#
# x + y + z = X
# x*y + x*z + y*z = B
# x*y*z = P
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
# x^9 - A*x^6 - 3*B*C*x^4 + (B^3+3*C^2)*x^3 -3*B^2*C*x^2 + 3*B*C^2*x - C^3 # = 0

# Test
x = sqrt(c(2,3,5))
A = sum(x^3)
B = sum(x[1]*x[2], x[1]*x[3], x[2]*x[3])
C = prod(x)
x^9 - A*x^6 - 3*B*C*x^4 + (B^3+3*C^2)*x^3 -3*B^2*C*x^2 + 3*B*C^2*x - C^3


### C.) Complex P3 System
# b31*(x^3+y^3+z^3) + b21*(x^2+y^2+z^2) + b11*(x+y+z) + e21*(xy+xz+yz) + e31*x*y*z = A1
# b32*(x^3+y^3+z^3) + b22*(x^2+y^2+z^2) + b12*(x+y+z) + e22*(xy+xz+yz) + e32*x*y*z = A2
# b33*(x^3+y^3+z^3) + b23*(x^2+y^2+z^2) + b13*(x+y+z) + e23*(xy+xz+yz) + e33*x*y*z = A3

# see below for complet *exact* solution!
# 1.) the system is decomposed using a liniar decomposition;
# 2.) then the "base"-terms are computed;
# 3.) then the P3 system is solved based on [A];


### D.) Perturbations to break Symmetry
# [C] is still symetrical;
# see below for some ideas;


### E.) Basic Assymetrical System:
# b11 * x + b12 * y + b13 * z = R1
# b21 * x*y + b22 * x*z + b23 * y*z = R2
# x*y*z = R3
# see part [E] for solution;
# [but invovles numerical solution to a P[6] polynomial;]

### E.2.) a greatly simplified version of [E]
# but with exact solution;


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


##########################

##########################
### E.) Basic Assymetrical
###     System

# (b11, b12, b13) * (x, y, z) = R1
# (b21, b22, b23) * E2 = R2
# x*y*z = R3

# => y*z = R3 / x
# (b12, b13) * (y, z) = - b11 * x + R1
# (b21, b22) * (y, z) = R2/x - b23 * R3/x^2
# => linear solution for y & z;


library(polynom)


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

