
### Leonard Mada
###
### P3 Polynomial Systems
### Solver: Exact solutions
###
### draft 0.2c

### P3 Systems
# v.02c: Test the Linear solution concept;
# v.02: P3 system + linear (x+y+z) terms;
# v.01: simple P3 system: the Base System;

#####################

### Base P3 system

### P3 System:
# x^3 + y^3 + z^3 = A
# x*y + x*z + y*z = B
# x*y*z = P


### Steps:

### Subsystem 1:
#
# X^3 - 3*B*X + 3*P - A = 0
# solve for X;

### Subsystem 2:
#
# x + y + z = X
# x*y + x*z + y*z = B
# x*y*z = P
# solve for x, y, z;


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
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(sum(isZero) > 0) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	return(m)
}

### Solution
# Cubic Polynomial
solveP3 <- function(c, d, n=3, all=FALSE) {
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
		print(paste("err Tri =", round0(err)))
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
	print(paste("y =", round(y, 5)))
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
solveP3S = function(s, bi, pr, n=3) {
	d2 = solveP3Outer(s, bi, pr, n)
	solveP3S.Base(s, bi, pr, d2, n)
}
solveP3SS = function(B, b, n=3) {
	# X^3 + 3*b21*X^2 - (3*BB + b31 + b11)*X + 3*P - A = 0
	r = solveP3OuterShifted(B, b, n)
	B.shifted = B - r * b
	solveP3S.Base(B.shifted[1], B.shifted[2], B.shifted[3], r, n, b.shifts=b)
}

############################

### Test

solveP3S(1,2,3)

solveP3S(1,2,-1)

# computes 3 roots by now;
# TODO: add all roots


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
x = x.all$x[1,]
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

### B.) Perturbations to break Symmetry

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

