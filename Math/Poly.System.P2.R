
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: P2
### v.0.2b


### History
# draft v.0.2b:
# - some examples of derived polynomials of order 8
#   based on the order 2 system & cos(2*pi/5)-entanglement;
# draft v.0.2:
# - added solutions for order 2 & order 4;
# draft v.0.1:
# - solution to the Symmetric order 3 system;
# - used to solve also:
#   x^6 + 3*x^5 + 3*x^4 + b*x^3 + 3*c*x^2 + 3*c^2*x + c^3 = 0;


######################

### helper functions

# TODO: use exact solver for P3
library(polynom)
# library(pracma) # for polynomials with complex coefficients

solve.ps = function(shift, R, n=3) {
	if(n == 2) {
		# this is only a quadratic
		coeff = c(2*shift^2 - R[1] - 2*R[2], 2*shift, 1)
	} else if(n == 3) {
		coeff = c(2*shift^3 - R[1] - 6*shift*R[2], 3*(shift^2-R[2]), 3*shift, 1)
	} else if(n == 4) {
		coeff = c(2*shift^4 - 12*shift^2*R[2] - R[1] + 2*R[2]^2,
			4*shift*(shift^2-3*R[2]), (6*shift^2-4*R[2]), 4*shift, 1)
	} else {
		print("NOT yet implemented!")
	}
	X = solve(polynomial(coeff))
	# X = roots(rev(coeff)) # pracma
	#
	det = X^2 - 4*R[2]
	x.m = sapply(det, function(det) if(Im(det) != 0 | Re(det) >= 0) sqrt(det) else complex(re=0, im=sqrt(-det)) )
	x = (X + x.m)/2
	y = (X - x.m)/2
	sol.df = data.frame(x=c(x, y), y=c(y, x))
	return(sol.df)
}
# round to 0
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

####################

####################
### Symmetric System

### Root Shift ###

### n = 2
# (x + a)^2 + (y + a)^2 = R1
# x*y = R2

# x^2 + y^2 + 2*a*(x + y) + 2*a^2 - R1 = 0
# X = x + y =>
# X^2 + 2*a*X + 2*a^2 - R1 - 2*R2 = 0
# Step 1: Solve for X
# Step 2: x + y = X

### Parameters
a = 1
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=2)

### Test
(sol$x + a)^2 + (sol$y + a)^2
sol$x * sol$y

###
a = -2
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=2)

### Test
(sol$x + a)^2 + (sol$y + a)^2
sol$x * sol$y
sol


### Equivalent Polynomial:
x = sol$x

err = x^2*(x + a)^2 - R[1]*x^2 + (R[2] + a*x)^2
round0(err)


###
c3 = 2*cos(2*pi/7 * (1:3))
R = c(1,1)
x = vapply(c3, solve.ps, R=R, n=2)


###############

### n = 3
# (x + a)^3 + (y + a)^3 = R1
# x*y = R2

# x^3 + y^3 + 3*a*(x^2 + y^2) + 3*a^2*(x + y) + 2*a^3 - R1 = 0
# X = x + y =>
# X^3 - 3*R2*X + 3*a*X^2 - 6*a*R2 + 3*a^2*X + 2*a^3 - R1 = 0
# X^3 + 3*a*X^2 + 3*(a^2-R2)*X + 2*a^3 - R1 - 6*a*R2 = 0
# Step 1: Solve for X
# Step 2: x + y = X


### Parameters
a = 1
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=3)

### Test
(sol$x + a)^3 + (sol$y + a)^3
sol$x * sol$y


### Equivalent Polynomial:
### Special Cases
b = 2*a^3 - R[1]
x = sol$x
err = x^6 + 3*a*x^5 + 3*a^2*x^4 + b*x^3 + 3*a^2*R[2]*x^2 + 3*a*R[2]^2*x + R[2]^3
round0(err)


# some special cases
### NO x^3
a = 1
c = 1
#
sol = solve.ps(a, c(2*a^3, c), n=3)
x = sol$x
err = x^6 + 3*a*x^5 + 3*a^2*x^4 + 3*a^2*c*x^2 + 3*a*c^2*x + c^3
round0(err)
sol

# some special cases
b = 3
c = 1
#
sol = solve.ps(1, c(2 - b, c), n=3)
x = sol$x
err = x^6 + 3*x^5 + 3*x^4 + b*x^3 + 3*c*x^2 + 3*c^2*x + c^3
round0(err)
sol

###
b = 1
c = 1
#
sol = solve.ps(1, c(2 - b, c), n=3)
x = sol$x
err = x^6 + 3*x^5 + 3*x^4 + b*x^3 + 3*c*x^2 + 3*c^2*x + c^3
round0(err)
sol


##################
##################

### n = 4
# (x + a)^4 + (y + a)^4 = R1
# x*y = R2

### TODO:
# (x + a)^4 + (y + a)^4 - R1 = 0
# x^4 + y^4 + 4*a*(x^3 + y^3) + 6*a^2*(x^2 + y^2) + 4*a^3*(x + y) + 2*a^4 - R1 = 0
# X = x + y =>
# X^4 - 4*R2*(x^2 + y^2) - 6*R2^2 +
#  + 4*a*X^3 - 12*a*R2*X + 6*a^2*X^2 - 12*a^2*R2 + 4*a^3*X + 2*a^4 - R1 = 0
#
# X^4 + 4*a*X^3 + (6*a^2-4*R2)*X^2 + 4*(a^3-3*a*R2)*X + 2*a^4 - R1 + 2*R2^2 - 12*a^2*R2 = 0
# Step 1: Solve for X
# Step 2: x + y = X


### Example 1:
# Parameters
a = 1
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=4)

### Test
(sol$x + a)^4 + (sol$y + a)^4
sol$x * sol$y
sol


### Example 2
# Parameters
a = 2
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=4)

### Test
(sol$x + a)^4 + (sol$y + a)^4
sol$x * sol$y
sol


### Example 3
# Parameters
a = -3
R = c(-2, 2)
### Solution
sol = solve.ps(a, R, n=4)

### Test
(sol$x + a)^4 + (sol$y + a)^4
sol$x * sol$y
sol



#########################
#########################

### Derived Polynomials

c2 = 2*cos(2*pi/5 * 1:2)

###
R1 = c(c2[1]-c2[2], 1)
R2 = c(c2[2]-c2[1], 1)
sol = rbind(
	solve.ps(c2[1] - c2[2], R1, n=2),
	solve.ps(c2[2] - c2[1], R2, n=2) )

# poly.calc(sol$x)
x = sol$x
1 + 20*x^3 + 57*x^4 + 20*x^5 + x^8
# Note: Set {x} == {1/x}


###
R1 = c(c2[2] - c2[1], 1)
R2 = c(c2[1] - c2[2], 1)
sol = rbind(
	solve.ps(c2[1]^3-c2[2]^3, R1, n=2),
	solve.ps(c2[2]^3-c2[1]^3, R2, n=2) )

# poly.calc(sol$x)
x = sol$x
1 - 40*x^3 + 1437*x^4 - 40*x^5 + x^8


###
R1 = c(0, c2[1])
R2 = c(0, c2[2])
sol = rbind(
	solve.ps(c2[1]-c2[2], R1, n=2),
	solve.ps(c2[2]-c2[1], R2, n=2) )

poly.calc(sol$x)
x = sol$x
1 + 10*x + 50*x^2 + 110*x^3 + 123*x^4 + 10*x^5 + x^8


###
R1 = c(c2[1]-c2[2], c2[1])
R2 = c(c2[2]-c2[1], c2[2])
sol = rbind(
	solve.ps(c2[1]-c2[2], R1, n=2),
	solve.ps(c2[2]-c2[1], R2, n=2) )

poly.calc(sol$x)
x = sol$x
1 + 10*x + 45*x^2 + 100*x^3 + 118*x^4 + 30*x^5 + x^8


###
R1 = c(c2[2]-c2[1], c2[1]-c2[2])
R2 = c(c2[1]-c2[2], c2[2]-c2[1])
sol = rbind(
	solve.ps(c2[1] - c2[2], R1, n=2),
	solve.ps(c2[2] - c2[1], R2, n=2) )

poly.calc(sol$x)
x = sol$x
25 + 100*x + 200*x^2 + 200*x^3 + 105*x^4 + x^8


###
R1 = c(0, c2[1]-c2[2])
R2 = c(0, c2[2]-c2[1])
sol = rbind(
	solve.ps(c2[1] - c2[2], R1, n=2),
	solve.ps(c2[2] - c2[1], R2, n=2) )

poly.calc(sol$x)
x = sol$x
25 + 100*x + 200*x^2 + 200*x^3 + 110*x^4 + 20*x^5 + x^8



