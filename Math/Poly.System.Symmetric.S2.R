
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: P2
### Decompositions of Symmetric Systems
### v.0.3h


####################
### Introduction ###
####################

# - various types of basic Symmetric Systems;


### A.1. Simple Symmetric
# x^n + y^n = R1
# x*y = R2

### A.2. Shifted Symmetric
# (x - d)^n + (y - d)^n = R1
# x*y = R2


### B. Non-Linearly Entangled

### B.1. Order 1 Ent
# x^n + y^n = R1
# x*y*(x+y) = R2

### B.2. Order 2 Ent
# x^n + y^n = R1
# x*y*(x+y)^2 = R2

### B.3. Div Ent
# x^n + y^n = R1
# x*y / (x+y) = R2


###############

###############
### History ###
###############

### draft v.0.3g - v.0.3h:
# - multiplicative & div entanglements for order 2;
#   [was skipped in v.0.2f - v.0.3a]
# - moved Section [Derived Polynomials] to the end;
### draft v.0.3e - v.0.3f:
# - 2-Shift systems [Order 3]:
#   same shift:
#   (x + d)^n + (y + d)^n = R1;
#   (x - d)^n + (y - d)^n = R2;
#   OR different shifts:
#   (x + d2)^n + (y + d2)^n = R2;
### draft v.0.3d - v.0.3d-cases:
# - classic Polynomials for P[3] M-type: P[6];
# - more particular cases for P[6]; (v.0.3d-cases)
### draft v.0.3a-v.0.3c:
# - systematic approach to entanglements:
#  -- multiplicative: x*y*(x+y) = R;
#  -- dividing: x*y/(x+y)^j = R;
#  -- basic order 5: multiplicative variant, x*y*(x+y) (in v.0.3c);
# - TODO: all variants (more variants in v.0.3b);
### draft v.0.2f:
# - entanglement: x*y*(x+y) = R;
### draft v.0.2e:
# - solved a symmetric "liniar" etension;
#   x^3 + y^3 + b1*(x+y) = R1;
#   x*y + b2*(x+y) = R2;
### draft v.0.2d:
# - solved: x^6 + b4*x^4 + b3*x^3 + b4*b2*x^2 + b2^3;
#   e.g. x^6 - x^4 - x^3 + x^2 - 1 = 0;
#   Note: this is a generalized symmetric poly;
#   [see file: Polynomials.Derived.P6.Symmetric.R]
# draft v.0.2b & v.0.2c:
# - some examples of derived polynomials of order 8
#   based on the order 2 system & cos(2*pi/5)-entanglement;
# draft v.0.2:
# - added solutions for order 2 & order 4;
# draft v.0.1:
# - solution to the Symmetric order 3 system;
# - used to solve also:
#   x^6 + 3*x^5 + 3*x^4 + b*x^3 + 3*c*x^2 + 3*c^2*x + c^3 = 0;
#   and
#   x^6 + 3*a*x^5 + 3*a^2*x^4 + b*x^3 + 3*a^2*x^2 + 3*a*x + 1 = 0;


######################

### helper functions

# TODO: use exact solver for P3
library(polynom)
library(pracma) # for polynomials with complex coefficients


solve.shift.ps = function(s1, s2, R, n=3) {
	# product is shifted as well
	sol = solve.ps(s1 - s2, R, n=n)
	return(sol - s2)
}
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
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}
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

########################

########################
### Symmetric System ###
########################

### Root Shift ###

### n = 2
# (x + a)^2 + (y + a)^2 = R1
# x*y = R2

### Solution:

x^2 + y^2 + 2*a*(x + y) + 2*a^2 - R1 # = 0
# S = x + y =>
S^2 + 2*a*S + 2*a^2 - R1 - 2*R2 # = 0
# Step 1: Solve for S
# Step 2: x + y = S

### Examples:

a = 1
R = c(1, 1)
### Solution
sol = solve.ps(a, R, n=2)

### Test
(sol$x + a)^2 + (sol$y + a)^2
sol$x * sol$y

### Ex 2:
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


####################
####################

### n = 3
# (x + a)^3 + (y + a)^3 = R1
# x*y = R2

### Solution:
x^3 + y^3 + 3*a*(x^2 + y^2) + 3*a^2*(x + y) + 2*a^3 - R1 # = 0
# X = x + y =>
X^3 - 3*R2*X + 3*a*X^2 - 6*a*R2 + 3*a^2*X + 2*a^3 - R1 # = 0
X^3 + 3*a*X^2 + 3*(a^2-R2)*X + 2*a^3 - R1 - 6*a*R2 # = 0
# Step 1: Solve for X
# Step 2: x + y = X

### Examples:

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

### Case a = 1
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

### Case c = 1
###
a = 1
b = 3
#
sol = solve.ps(a, c(2*a^3 - b, 1), n=3)
x = sol$x
err = x^6 + 3*a*x^5 + 3*a^2*x^4 + b*x^3 + 3*a^2*x^2 + 3*a*x + 1
round0(err)
sol

###
a = 1
b = -3
#
sol = solve.ps(a, c(2*a^3 - b, 1), n=3)
x = sol$x
err = x^6 + 3*a*x^5 + 3*a^2*x^4 + b*x^3 + 3*a^2*x^2 + 3*a*x + 1
round0(err)
sol

###
a = 2
b = 8
#
sol = solve.ps(a, c(2*a^3 - b, 1), n=3)
x = sol$x
err = x^6 + 3*a*x^5 + 3*a^2*x^4 + b*x^3 + 3*a^2*x^2 + 3*a*x + 1
round0(err)
sol
# shift origin
x = sol$x + 1
6 - 18*x + 15*x^2 - 3*x^4 + x^6


##################
##################

### n = 4
# (x + a)^4 + (y + a)^4 = R1
# x*y = R2

### Solution:
(x + a)^4 + (y + a)^4 - R1 # = 0
x^4 + y^4 + 4*a*(x^3 + y^3) + 6*a^2*(x^2 + y^2) + 4*a^3*(x + y) + 2*a^4 - R1 # = 0
# X = x + y =>
X^4 - 4*R2*(x^2 + y^2) - 6*R2^2 +
	+ 4*a*X^3 - 12*a*R2*X + 6*a^2*X^2 - 12*a^2*R2 + 4*a^3*X + 2*a^4 - R1 # = 0

X^4 + 4*a*X^3 + (6*a^2-4*R2)*X^2 + 4*(a^3-3*a*R2)*X + 2*a^4 - R1 + 2*R2^2 - 12*a^2*R2 # = 0
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


### Equivalent Polynomial:
x = sol$x
x^4*(x + a)^4 + (x*a + R[2])^4 - R[1]*x^4 # = 0


############################
############################

### Liniar/Polynomial Shifts

### Order 3
# x^3 + y^3 + b1*(x + y) = R1
# x*y = R2

### Solution
### Step 1:
# s = x + y; c = x*y = R2 =>
# s^3 - (3*c - b1)*s - R1 = 0
### Step 2:
# x + y = s
# x*y = R2

b = c(-1)
R = c(1, -1)
### Solution
s = roots(c(1, 0, - (3*R[2] - b[1]), -R[1]))
sm = sqrt(s^2 - 4*R[2] + 0i)
x = (s + sm)/2
y = (s - sm)/2
sol = cbind(x,y)
sol = rbind(sol, cbind(y,x))
sol

### Test
x^3 + y^3 + b[1]*(x + y)
x*y

### Classic
x = sol[,1]
x^6 + b[1]*x^4 - R[1]*x^3 + b[1]*R[2]*x^2 + R[2]^3

# for b1 = -1; R = c(1, -1):
x^6 - x^4 - x^3 + x^2 - 1


###################

###
# x^3 + y^3 + b1*(x + y) = R1
# x*y + b2*(x + y) = R2

### Solution
### Step 1:
# s = x + y; c = x*y = R2
# => c = R2 - b2*s
# s^3 - (3*c - b1)*s - R1 = 0
# s^3 - (3*R2 - 3*b2*s - b1)*s - R1 = 0
# s^3 + 3*b2*s^2 - (3*R2 - b1)*s - R1 = 0
### Step 2:
# x + y = s
# x*y = R2 - b2*s

solve.p2p3liniar = function(b, R) {
	b = as.vector(unlist(b))
	p.coeff = c(1, 3*b[2], - (3*R[2] - b[1]), -R[1])
	# print(p.coeff)
	s = roots(p.coeff)
	xy = R[2] - b[2]*s
	sm = sqrt(s^2 - 4*xy + 0i)
	x = (s + sm)/2
	y = (s - sm)/2
	sol = cbind(x,y)
	sol = rbind(sol, cbind(y,x))
	# coeffs Classic
	coeffs = c(1, 3*b[2], (b[1]+3*b[2]^2), (2*b[1]*b[2]-R[1]),
		(b[1]*b[2]^2 + b[1]*R[2] + 3*b[2]^2*R[2] - 3*b[2]*R[1]),
		(2*b[1]*b[2]*R[2] - 3*b[2]*R[2]^2 - 3*b[2]^2*R[1]),
		R[2]^3 + b[1]*b[2]^2*R[2] - b[2]^3*R[1] )
	### Test
	eq1 = x^3 + y^3 + b[1]*(x + y)
	eq2 = x*y + b[2]*(x + y)
	#
	return(list(sol=sol, coeffs=coeffs, test=rbind(round0(eq1)), round0(eq2)))
}


############
### Examples

### Ex. 1
b = c(-3, 1)
R = c(1, -1)
### Solution
sol = solve.p2p3liniar(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x + y)
x*y + b[2]*(x + y)


### Classic
x = sol$sol[,1]
coeffs = c(1, 3*b[2], (b[1]+3*b[2]^2), (2*b[1]*b[2]-R[1]),
  (b[1]*b[2]^2 + b[1]*R[2] + 3*b[2]^2*R[2] - 3*b[2]*R[1]),
  (2*b[1]*b[2]*R[2] - 3*b[2]*R[2]^2 - 3*b[2]^2*R[1]),
  R[2]^3 + b[1]*b[2]^2*R[2] - b[2]^3*R[1] )
coeffs
x^6 + 3*b[2]*x^5 + (b[1]+3*b[2]^2)*x^4 + (2*b[1]*b[2]-R[1])*x^3 + (b[1]*b[2]^2 + b[1]*R[2] + 3*b[2]^2*R[2] - 3*b[2]*R[1])*x^2 + (2*b[1]*b[2]*R[2] - 3*b[2]*R[2]^2 - 3*b[2]^2*R[1])*x + R[2]^3 + b[1]*b[2]^2*R[2] - b[2]^3*R[1]
#
x^6 + 3*x^5 - 7*x^3 - 6*x^2 - 3*x - 3

### Derivation
x^3*(x + b[2])^3 + b[1]*(x^2 + R[2])*(x + b[2])^2 - R[1]*(x + b[2])^3 + (R[2]-b[2]*x)^3 # == 0

#########
### Ex. 2
b = c(-3, 1)
R = c(1, -2)
### Solution
sol = solve.p2p3liniar(b, R)
sol


###
bg = expand.grid(-3:3, -3:3)
R = c(1, 1)
#
p.l = sapply(1:nrow(bg), function(id) print(solve.p2p3liniar(bg[id,], R)$coeffs))
# contains various trivial examples
# like x^6 + b1*x^5 + b2*x^4 + b2*x^2 - b1*x + 1; x = {1i, -1i, ...}
# but also many non-trivial examples;


### trivial example
b = c(0,-1)
R = c(0, 1)
sol = solve.p2p3liniar(b, R)
sol
x = sol$sol[,1]
x^6 - 3*x^5 + 3*x^4 + 3*x^2 + 3*x + 1


### TODO: all variants
# - liniar/polynomial shifts;
# - root shifts;



############################
############################

############################
### Non-Liniar Entanglements

###############
### Order 2 ###
###############


### x^2 + y^2 + b*(x+y) = R1
### x*y*(x+y) = R2

### Solution:

### Eq 1:
S^2 - 2*x*y + b*S - R1 # = 0
S^3 + b*S^2 - R1*S - 2*R2 # = 0

### Solver:
solve.Ent.S2P2 = function(R, b, debug=TRUE) {
	coeff = c(1, b[1], - R[1], - 2*R[2])
	S = roots(coeff)
	if(debug) print(S);
	#
	xy = R[2] / S;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	sol = rbind(sol, sol[,2:1])
	return(sol)
}

### Examples:

R = c(-2, 2)
b = 2;
sol = solve.Ent.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + y^2 + b[1]*(x+y) # - R[1]
x*y*(x+y) # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x))

R1 = R[1]; R2 = R[2];
R2^2 +
	- R2*(2*R1 + b^2)*x +
	+ (R1^2 + 3*b*R2)*x^2 +
	- (R1*b - 2*R2)*x^3 +
	- 3*R1*x^4 + 2*b*x^5 + 2*x^6


#########
### Ex 2:
R = c(-2, -4)
b = 2
sol = solve.Ent.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + y^2 + b[1]*(x+y) # - R[1]
x*y*(x+y) # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x))


###############

### Div-Variant
### (x + y)^3 + b2*x*y + b1*(x+y) = R1
### x*y / (x+y)^2 = R2

### Solution:

### Eq 1:
S^3 + b2*R2*S^2 + b1*S - R1 # = 0


### Solver:
solve.Ent2.S2 = function(R, b=c(0, 1), debug=TRUE) {
	coeff = c(1, b[2]*R[2], b[1], - R[1])
	S = roots(coeff)
	if(debug) print(S);
	#
	xy = R[2] * S^2;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	sol = rbind(sol, sol[,2:1])
	return(sol)
}
test.Ent2.S2 = function(sol, R, b=c(0, 1)) {
	x = sol[,1]; y = sol[,2];
	err1 = (x + y)^3 + b[2]*x*y + b[1]*(x+y) # - R[1]
	err2 = x*y / (x+y)^2 # - R[2]
	err = round0(rbind(err1, err2))
	return(err)
}

### Examples:

R = c(-2, 1)
b = c(1, 1);
sol = solve.Ent2.S2(R, b)
x = sol[,1]; y = sol[,2];

### Test
test.Ent2.S2(sol, R, b)

### Classic Polynomial:
round0.p(poly.calc(x))
err = 4 + 2*x - x^2 - 3*x^3 + x^5 + x^6
round0(err)


#########
### Ex 2:
R = c(-4, 1)
b = c(4, 2);
sol = solve.Ent2.S2(R, b)
x = sol[,1]; y = sol[,2];

### Test
test.Ent2.S2(sol, R, b)

### Classic Polynomial:
round0.p(poly.calc(x))
err = 16 + 16*x + 8*x^2 + 2*x^5 + x^6
round0(err)


R1 = R[1]; R2 = R[2]; b1 = b[1]; b2 = b[2];
(R1^2*R2^3) +
- b1*R1*R2^2*x^1 +
+ R2^2*(2*R1*R2*b2 - R1*b2 + b1^2)*x^2 +
+ (3*R1*R2 + b1*b2*R2^2 - R1)*x^3 +
+ (R2^3*b2^2 - 2*R2*b1 + b1)*x^4 +
+ b2*R2*x^5 + x^6


###############

###############
### Order 3 ###
###############

### Type Multiplicative:
# x^3 + y^3 + b1*(x+y) = R1
# x*y*(x+y) = R2

### Type Div:
# x^3 + y^3 + b1*(x+y) = R1
# x*y / (x+y) = R2

### Extensions:
# x^3 + y^3 + b3*(x+y)^2 + b1*(x+y) = R1
# x*y {*,/} (x+y) + b2*(x+y) = R2

### Solution:

### Step 1:
# S = x + y =>
S^3 + b1*S - R1 - 3*R2 # = 0
### Step 2:
# x + y = S
# x*y = R2/S

solve.p2p3ent = function(b, R, type="mult", n=3) {
	# type "mult": x*y*(x+y) + b[2]*(x+y) = R[2]
	# type "div":  x*y/(x+y) + b[2]*(x+y) = R[2]
	# type "div2": x*y/(x+y)^2 + b[2]*(x+y) = R[2]
	# simple type: b[2] = 0
	if(length(b) < 2) { b = c(b, 0) }
	if(length(b) < 3) { b = c(b, 0) }
	# Solve:
	if(type == "mult" && n == 3) {
		s = roots(c(1, b[3], b[1] + 3*b[2], - R[1] - 3*R[2]))
		xy = R[2]/s - b[2]
	} else if(type == "div") {
		s = roots(c(1 + 3*b[2], b[3] - 3*R[2], b[1], - R[1]))
		xy = R[2] * s - b[2] * s^2
	} else if(type == "div2") {
		s = if(b[2] == 0) roots(c(1 - 3*R[2], b[3], b[1], - R[1]))
			else roots(c(3*b[2], 1 - 3*R[2], b[3], b[1], - R[1]))
		xy = R[2] * s^2 - b[2] * s^3
	} else if(type == "mult" && n == 2) {
		# TODO: b[3]
		s = roots(c(1, b[1], 2*b[2] - R[1], -2*R[2]))
		xy = R[2]/s - b[2]
	} else {
		stop("Type NOT yet supported!")
	}
	s.diff = sqrt(s^2 - 4*xy + 0i)
	x = (s + s.diff)/2
	y = (s - s.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	### Test
	S = (x+y);
	t1 = x^n + y^n + b[1]*S + b[3]*S^2;
	t2 = if(type == "mult") { x*y*S + b[2]*S; }
	else if(type == "div") { x*y/S + b[2]*S; }
	else if(type == "div2"){ x*y/S^2 + b[2]*S; }
	#
	return(list(sol=sol, test=rbind(t1, t2)))
}

### Example 1
b = c(-1)
R = c(-2, 1)
#
sol = solve.p2p3ent(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x+y)
x*y*(x+y)

### Classic
poly.calc(x) # trivial: (x^3 - x + 1)^2
1 - 2*x + x^2 + 2*x^3 - 2*x^4 + x^6


### the 1-coeff variants are trivial
b = c(2)
R = c(-2, 1)
#
sol = solve.p2p3ent(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x+y)
x*y*(x+y)

### Classic
poly.calc(x) # (x^3 + 2*x + 1)^2
1 + 4*x + 4*x^2 + 2*x^3 + 4*x^4 + x^6

### Classic Polynomials:

### Simple + Extensions A1:
# b2 = 0;
(R2^3) +
R2^2*(2*b1 - b3^2)*x^1 +
R2*(2*R1*b3 + b1^2 + 5*R2*b3)*x^2 +
(- 3*R1*R2 - R1^2 + R2*b1*b3)*x^3 +
b1*(R1 + 4*R2)*x^4 +
b3*(R1 + 3*R2)*x^5 +
(R1 + 3*R2)*x^6

### b2 != 0
(- R1*b2^3 + R2*b1*b2^2 + R2^2*b2*b3 + R2^3) +
(- 3*R1*R2*b2 + R1*b2^2*b3 - R2*b1*b2*b3 + 2*R2^2*b1 - 3*R2^2*b2 - R2^2*b3^2)*x^1 +
(2*R1*R2*b3 - R1*b1*b2 + R2*b1*b2 + R2*b1^2 + 3*R2*b2^2 + 5*R2^2*b3)*x^2 +
(- 3*R1*R2 - 2*R1*b2*b3 - R1^2 + R2*b1*b3 - 3*R2*b2*b3)*x^3 +
(R1*b1 + 4*R2*b1 + 3*R2*b2)*x^4 +
b3*(R1 + 3*R2)*x^5 +
(R1 + 3*R2)*x^6

### Special Cases:
### R1 = -2*R2
(R2^2 + R2*b2*b3 + b1*b2^2 + 2*b2^3) +
(2*R2*b1 + 3*R2*b2 - R2*b3^2 - b1*b2*b3 - 2*b2^2*b3)*x +
(R2*b3 + b1^2 + 3*b1*b2 + 3*b2^2)*x^2 +
(2*R2 + (b1 + b2)*b3)*x^3 +
(2*b1 + 3*b2)*x^4 +
b3*x^5 + x^6
### b1 = -3/2 * b2
(2*R2^2 + 2*R2*b2*b3 + b2^3) +
- (2*R2*b3^2 + b2^2*b3)*x +
1/2 * (4*R2*b3 + 3*b2^2)*x^2 +
(4*R2 - b2*b3)*x^3 + 2*b3*x^5 + 2*x^6
### b1 = b2 = 0
R2^2 - R2*b3^2*x + R2*b3*x^2 + 2*R2*x^3 + b3*x^5 + x^6
(x^3 + R2)^2 + b3*x^2*(x^3 + R2) - R2*b3^2*x

### R1 = -4*R2
(- R2^2 - R2*b2*b3 - b1*b2^2 - 4*b2^3) +
(- 2*R2*b1 - 9*R2*b2 + R2*b3^2 + b1*b2*b3 + 4*b2^2*b3)*x +
(3*R2*b3 - 5*b1*b2 - b1^2 - 3*b2^2)*x^2 +
(4*R2 - b1*b3 - 5*b2*b3)*x^3 +
- 3*b2*x^4 + b3*x^5 + x^6


#########
### Ex 1:
R = c(4, -1)
b = c(-1, 0, 4)
#
sol = solve.p2p3ent(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x+y) + b[3]*(x+y)^2
x*y*(x+y)

### Classic
round0.p(poly.calc(x))
err = -1 - 18*x - 13*x^2 + 4*x^5 + x^6
round0(err)


#########
### Ex 2:
R = c(-2, 1)
b = c(-1, 0, 2)
#
sol = solve.p2p3ent(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x+y) + b[3]*(x+y)^2
x*y*(x+y)

### Classic
poly.calc(x)
err = 1 - 6*x + 3*x^2 - 2*x^4 + 2*x^5 + x^6
round0(err)


#########
### Ex 3:
b3 = 6;
R = c(1/6 * b3^3, -1/12 * b3^3)
b = c(1/2*b3^2, -1/3 * b3^2, b3)
#
sol = solve.p2p3ent(b, R)
sol

### Test
x = sol$sol[,1]; y = sol$sol[,2]
x^3 + y^3 + b[1]*(x+y) + b[3]*(x+y)^2
x*y*(x+y) + b[2]*(x+y)

### Classic
round0.p(poly.calc(x))
# ...*b3^6 + 1/36*b3^5*x + b3*x^5 + x^6
err = 21*6^2 + 6^3*x + 6*x^5 + x^6
round0(err)


########
### only 1 Coeff: Trivial
b = c(2)
R = c(-2, 1)
p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(b, R)$sol[,1]))))
# (x^3 + b*x + 1)^2
1 - 10*x + 25*x^2 + 2*x^3 - 10*x^4 + x^6 
1 - 8*x + 16*x^2 + 2*x^3 - 8*x^4 + x^6 
1 - 6*x + 9*x^2 + 2*x^3 - 6*x^4 + x^6 
1 - 4*x + 4*x^2 + 2*x^3 - 4*x^4 + x^6 
1 - 2*x + x^2 + 2*x^3 - 2*x^4 + x^6 
1 - 0 + 0 + 2*x^3 - 0 + x^6 
1 + 2*x + x^2 + 2*x^3 + 2*x^4 + x^6 
1 + 4*x + 4*x^2 + 2*x^3 + 4*x^4 + x^6 
1 + 6*x + 9*x^2 + 2*x^3 + 6*x^4 + x^6 
1 + 8*x + 16*x^2 + 2*x^3 + 8*x^4 + x^6 
1 + 10*x + 25*x^2 + 2*x^3 + 10*x^4 + x^6 
# parametric example (but trivial)
b = 3
R = c(-2, 1); # fixed
p = solve.p2p3ent(b, R)
x = p$sol[,1]
err = 1 + 2*b*x + b^2*x^2 + 2*x^3 + 2*b*x^4 + x^6
round0(err)

### (x^3 + b*x + 1)^2 - b*x^4 - b*x
R = c(-2, 1)
p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(c(2*b, -b), R)$sol[,1]))))
1 - 5*x + 25*x^2 + 2*x^3 - 5*x^4 + x^6 
1 - 4*x + 16*x^2 + 2*x^3 - 4*x^4 + x^6 
1 - 3*x + 9*x^2 + 2*x^3 - 3*x^4 + x^6 
1 - 2*x + 4*x^2 + 2*x^3 - 2*x^4 + x^6 
1 - x + x^2 + 2*x^3 - x^4 + x^6 
1 - 0 + 0 + 2*x^3 - 0 + x^6 
1 + x + x^2 + 2*x^3 + x^4 + x^6 
1 + 2*x + 4*x^2 + 2*x^3 + 2*x^4 + x^6 
1 + 3*x + 9*x^2 + 2*x^3 + 3*x^4 + x^6 
1 + 4*x + 16*x^2 + 2*x^3 + 4*x^4 + x^6 
1 + 5*x + 25*x^2 + 2*x^3 + 5*x^4 + x^6 
# parametric example
b = 3
b.p = c(2*b, -b)
R = c(-2, 1); # fixed
p = solve.p2p3ent(b.p, R)
x = p$sol[,1]
err = 1 + b*x + b^2*x^2 + 2*x^3 + b*x^4 + x^6
round0(err)


###
R = c(-4,1)
p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(b, R, type="mult")$sol[,1]))))
# x^6 - R[1]*x^3 - (b*x - 1)^2; b0 = 1/(R[1]+3), more terms for other R = c(..., 1)
-1 + 10*x - 25*x^2 + 4*x^3 + x^6 
-1 + 8*x - 16*x^2 + 4*x^3 + x^6 
-1 + 6*x - 9*x^2 + 4*x^3 + x^6 
-1 + 4*x - 4*x^2 + 4*x^3 + x^6 
-1 + 2*x - x^2 + 4*x^3 + x^6 
-1 + 0 + 4*x^3 + x^6 
-1 - 2*x - x^2 + 4*x^3 + x^6 
-1 - 4*x - 4*x^2 + 4*x^3 + x^6 
-1 - 6*x - 9*x^2 + 4*x^3 + x^6 
-1 - 8*x - 16*x^2 + 4*x^3 + x^6
-1 - 10*x - 25*x^2 + 4*x^3 + x^6

###
R = c(-2,1)
p = sapply(-8:8, function(b) print(round0.p(poly.calc(solve.p2p3ent(b, R, type="div")$sol[,1]))))
#
-2 + 6*x - 12*x^2 + 14*x^3 - 3*x^4 - 3*x^5 + x^6 
-2 + 6*x - 11*x^2 + 12*x^3 - 2*x^4 - 3*x^5 + x^6 
-2 + 6*x - 10*x^2 + 10*x^3 - x^4 - 3*x^5 + x^6 
-2 + 6*x - 9*x^2 + 8*x^3 - 0 - 3*x^5 + x^6 
-2 + 6*x - 8*x^2 + 6*x^3 + x^4 - 3*x^5 + x^6 
-2 + 6*x - 7*x^2 + 4*x^3 + 2*x^4 - 3*x^5 + x^6 
-2 + 6*x - 6*x^2 + 2*x^3 + 3*x^4 - 3*x^5 + x^6 
-2 + 6*x - 5*x^2 + 0 + 4*x^4 - 3*x^5 + x^6 
-2 + 6*x - 4*x^2 - 2*x^3 + 5*x^4 - 3*x^5 + x^6 
-2 + 6*x - 3*x^2 - 4*x^3 + 6*x^4 - 3*x^5 + x^6 
-2 + 6*x - 2*x^2 - 6*x^3 + 7*x^4 - 3*x^5 + x^6 
-2 + 6*x - x^2 - 8*x^3 + 8*x^4 - 3*x^5 + x^6 
-2 + 6*x - 0 - 10*x^3 + 9*x^4 - 3*x^5 + x^6 
-2 + 6*x + x^2 - 12*x^3 + 10*x^4 - 3*x^5 + x^6


###
R = c(-2,1)
p = sapply(-8:8, function(b) print(round0.p(poly.calc(solve.p2p3ent(b, R, type="div2")$sol[,1]))))
#
1 - 4*x + 16*x^2 + 2*x^3 - 4*x^4 + x^6 
1 - 3.5*x + 12.25*x^2 + 2*x^3 - 3.5*x^4 + x^6 
1 - 3*x + 9*x^2 + 2*x^3 - 3*x^4 + x^6 
1 - 2.5*x + 6.25*x^2 + 2*x^3 - 2.5*x^4 + x^6 
1 - 2*x + 4*x^2 + 2*x^3 - 2*x^4 + x^6 
1 - 1.5*x + 2.25*x^2 + 2*x^3 - 1.5*x^4 + x^6 
1 - x + x^2 + 2*x^3 - x^4 + x^6 
1 - 0.5*x + 0.25*x^2 + 2*x^3 - 0.5*x^4 + x^6 
1 - 0 + 0 + 2*x^3 - 0 + x^6
1 + 0.5*x + 0.25*x^2 + 2*x^3 + 0.5*x^4 + x^6 
1 + x + x^2 + 2*x^3 + x^4 + x^6 
1 + 1.5*x + 2.25*x^2 + 2*x^3 + 1.5*x^4 + x^6 
1 + 2*x + 4*x^2 + 2*x^3 + 2*x^4 + x^6 
1 + 2.5*x + 6.25*x^2 + 2*x^3 + 2.5*x^4 + x^6 
1 + 3*x + 9*x^2 + 2*x^3 + 3*x^4 + x^6 
1 + 3.5*x + 12.25*x^2 + 2*x^3 + 3.5*x^4 + x^6 
1 + 4*x + 16*x^2 + 2*x^3 + 4*x^4 + x^6


###
p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(c(b, 2), c(2, 2), type="mult", n=2)$sol[,1]))))

p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(c(b, b), c(2, 2), type="mult", n=2)$sol[,1]))))

p = sapply(-6:6, function(b) print(round0.p(poly.calc(solve.p2p3ent(c(b, b), c(4, 4), type="mult", n=2)$sol[,1]))))

# b0 = b^2 / 2
p = sapply( (-6:6)[-7], function(b) print(round0.p(poly.calc(solve.p2p3ent(c(-b, -b), c(b, b), type="mult", n=2)$sol[,1]))))

sol = solve.p2p3ent(c(-2, -2), c(2, 2), type="mult", n=2)
x = sol$sol[,1]
2 - 3*x^4 - 2*x^5 + x^6

sol = solve.p2p3ent(c(-4, -4), c(4, 4), type="mult", n=2)
x = sol$sol[,1]
8 - 4*x^3 - 6*x^4 - 4*x^5 + x^6


#######

###
sol = solve.p2p3ent(c(0,-1/3), c(9,1/3), type="div2")
x = sol$sol[,1]
sol
poly.calc(x)
10 + 12*x + 30*x^2 - x^4 - 12*x^5 + x^8
###
sol = solve.p2p3ent(c(0, 1/3), c(9,1/3), type="div2")
sol
x = sol$sol[,1]
poly.calc(x)
-8 - 12*x + 30*x^2 + x^4 - 12*x^5 + x^8

###
sol = solve.p2p3ent(c(0, -1/3), c(9, -2/3), type="div2")
x = sol$sol[,1]
1 - 6*x + 21*x^2 - 27*x^3 + 29*x^4 - 12*x^5 + 3*x^6 - 3*x^7 + x^8
###
sol = solve.p2p3ent(c(0, 1/3), c(9, -2/3), type="div2")
x = sol$sol[,1]
-17 + 6*x + 21*x^2 - 27*x^3 - 29*x^4 - 12*x^5 + 3*x^6 + 3*x^7 + x^8



#############


### more complex example

m = complex(re=cos(2*pi/3), im=sin(2*pi/3))
m = m^(0:2)

K = 2
# works, but often with fractions
k = K^(1/3)
k = k*m

b = cbind(4*(k^2 - k), 0)
R = c(2, 2)
sol = as.vector(sapply(1:3, function(id) solve.p2p3ent(b[id,], R, type="mult", n=3)$sol[,1]))
poly.calc(sol)
#
x = sol
1 + 48*x^2 + 26*x^3 + 576*x^4 + 1008*x^5 + 823*x^6 - 1152*x^7 + 438*x^8 +  
724*x^9 + 576*x^10 + 132*x^11 + 145*x^12 + 102*x^14 - 6*x^15 + x^18


b = cbind(4*(k^2 - k), 2*k)
R = c(2, 2)
sol = as.vector(sapply(1:3, function(id) solve.p2p3ent(b[id,], R, type="mult", n=3)$sol[,1]))
poly.calc(sol)
# still one fraction
x = sol
-307 - 1728*x + 3120*x^2 - 9586*x^3 + 7680*x^4 - 7920*x^5 + 4486*x^6 +  
+ 2517*x^8 + 1003*x^9 + 1200*x^10 + 366*x^11 + 300.25*x^12 + 57*x^14 - 6*x^15 + x^18



################################
################################

###############
### Order 5 ###
###############

### Simple/Basic System
# x^5 + y^5 = R1
# x*y*(x+y) = R2

# (x+y)^5 - 5*x*y*(x+y)^3 + 5*(x*y)^2*(x+y) - R1 = 0
# X = x + y =>
# X^5 - 5*R2*X^2 + 5*R2^2/X - R1
# X^6 - 5*R2*X^3 - R1*X + 5*R2^2

R = c(25, 5)
x.sum = roots(c(1,0,0, - 5*R[2], 0, - R[1], 5*R[2]^2))
xy = R[2] / x.sum
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^5 + y^5
x*y*(x+y)

### Classic
round0.p(poly.calc(sol[,1]))
err = 125 - 125*x^4 - 25*x^5 - 25*x^7 + 5*x^9 + x^10 + x^12
round0(err)

### TODO:
# - parametric classic polynomial;


################################
################################

########################
### Symmetric System ###
### Variants         ###
########################

####################
### Double Shift ###
####################

### Order 3: Simple
# - same shift;
# - parameter d is specified;

# (x + d)^3 + (y + d)^3 = R1
# (x - d)^3 + (y - d)^3 = R2

### Solution:

### Diff =>
6*d*(x^2 + y^2) + 4*d^3 - R1 + R2 # = 0
6*d*(S^2 - 2*x*y) + 4*d^3 - R1 + R2 # = 0
6*d*S^2 - 12*d*x*y + 4*d^3 - R1 + R2 # = 0
# 12*d*x*y = 6*d*S^2 + 4*d^3 - R1 + R2

### Sum =>
2*(x^3 + y^3) + 6*d^2*(x+y) - R1 - R2 # = 0
2*S^3 - 6*x*y*S + 6*d^2*S - R1 - R2 # = 0
4*d*S^3 - 12*d*x*y*S + 12*d^3*S - 2*d*R1 - 2*d*R2 # = 0
4*d*S^3 - (6*d*S^2 + 4*d^3 - R1 + R2)*S + 12*d^3*S - 2*d*R1 - 2*d*R2 # = 0
2*d*S^3 - (8*d^3 + R1 - R2)*S + 2*d*R1 + 2*d*R2 # = 0

### Solver
solve.Shift2.S2P3 = function(R, d, debug=TRUE) {
	coeff = c(2*d, 0, - (8*d^3 + R[1] - R[2]), 2*d*R[1] + 2*d*R[2])
	S = roots(coeff)
	if(debug) print(S);
	#
	R1 = R[1]; R2 = R[2];
	xy = (6*d*S^2 + 4*d^3 - R1 + R2) / 12 / d;
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d) / 2;
	y = (S - xy.d) / 2;
	sol = cbind(as.vector(x), as.vector(y))
	return(rbind(sol, sol[,2:1]));
}

# - trivial roots:
#   if R1 = 0 or R2 = 0: S = -2*d or + 2*d;
#   if R1 + R2 = 0: S = 0;

### Examples:

R = c(1,2)
d = 1
sol = solve.Shift2.S2P3(R, d)
x = sol[,1]; y = sol[,2];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x)) * 16*27


#########
### Ex 2:
R = c(-2, 1)
d = 2
sol = solve.Shift2.S2P3(R, d)
x = sol[,1]; y = sol[,2];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x)) * 128 * 27


##########################

### Order 3:
# - different shifts;
# - parameters d are specified;

# (x + d1)^3 + (y + d1)^3 = R1
# (x + d2)^3 + (y + d2)^3 = R2

### Solution:

### Eq 1:
(x^3 + y^3) + 3*d1*(x^2 + y^2) + 3*d1^2*(x + y) + 2*d1^3 - R1 # = 0
S^3 - 3*x*y*S + 3*d1*(S^2 - 2*x*y) + 3*d1^2*S + 2*d1^3 - R1 # = 0
S^3 - 3*x*y*S - 6*d1*x*y + 3*d1*S^2 + 3*d1^2*S + 2*d1^3 - R1

### Eq 2:
S^3 - 3*x*y*S + 3*d2*(S^2 - 2*x*y) + 3*d2^2*S + 2*d2^3 - R2 # = 0

### Eq 1 - Eq 2 =>
3*(d1 - d2)*(S^2 - 2*x*y) + 3*(d1^2 - d2^2)*S + 2*d1^3 - 2*d2^3 - R1 + R2 # = 0
3*(d1 - d2)*S^2 - 6*(d1 - d2)*x*y + 3*(d1^2 - d2^2)*S + 2*d1^3 - 2*d2^3 - R1 + R2
# 6*(d1 - d2)*x*y = 3*(d1 - d2)*S^2 + 3*(d1^2 - d2^2)*S + 2*d1^3 - 2*d2^3 - R1 + R2

### =>
2*(d1 - d2)*S^3 - 6*(d1 - d2)*x*y*S - 2*d1*6*(d1 - d2)*x*y +
	+ 6*d1*(d1 - d2)*S^2 + 6*d1^2*(d1 - d2)*S + 4*(d1 - d2)*d1^3 - 2*(d1 - d2)*R1
2*(d1 - d2)*S^3 - (3*(d1 - d2)*S^2 + 3*(d1^2 - d2^2)*S + 2*d1^3 - 2*d2^3 - R1 + R2)*S +
	- 2*d1*(3*(d1 - d2)*S^2 + 3*(d1^2 - d2^2)*S + 2*d1^3 - 2*d2^3 - R1 + R2) +
	+ 6*d1*(d1 - d2)*S^2 + 6*d1^2*(d1 - d2)*S + 4*(d1 - d2)*d1^3 - 2*(d1 - d2)*R1
(d1 - d2)*S^3 + 3*(d1 - d2)*(d1 + d2)*S^2 + (2*d1^3 - 2*d2^3 + 6*d1^2*d2 - 6*d1*d2^2 - R1 + R2)*S +
	+ 4*d1*d2*(d1 - d2)*(d1 + d2) - 2*d2*R1 + 2*d1*R2

### Solver
solve.Shift2.S2P3 = function(R, d, debug=TRUE) {
	isSimple = (length(d) == 1) || (d[1] == -d[2]);
	if(isSimple) {
		coeff = c(2*d, 0, - (8*d^3 + R[1] - R[2]), 2*d*R[1] + 2*d*R[2])
	} else {
		if(d[1] == d[2]) stop("System is NOT solvable: d1 == d2!");
		d.diff = d[1] - d[2];
		d.sum = d[1] + d[2];
		coeff = c(d.diff, 3*d.diff*d.sum,
			(2*d.diff*(d.sum^2 + 2*d[1]*d[2]) - R[1] + R[2]),
			4*d[1]*d[2]*d.diff*d.sum - 2*d[2]*R[1] + 2*d[1]*R[2])
	}
	S = roots(coeff)
	if(debug) print(S);
	#
	R1 = R[1]; R2 = R[2];
	if(isSimple) {
		xy = (6*d*S^2 + 4*d^3 - R1 + R2) / 12 / d;
	} else {
		xy = (3*d.diff*S^2 + 3*d.diff*d.sum*S + 2*(d[1]^3 - d[2]^3) - R1 + R2)
		xy = xy / 6 / d.diff;
	}
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d) / 2;
	y = (S - xy.d) / 2;
	sol = cbind(as.vector(x), as.vector(y))
	return(rbind(sol, sol[,2:1]));
}

### Examples:

R = c(1, 2)
d = c(-1, 2)
sol = solve.Shift2.S2P3(R, d)
x = sol[,1]; y = sol[,2];

### Test
(x + d[1])^3 + (y + d[1])^3 # - R[1]
(x + d[2])^3 + (y + d[2])^3 # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x)) * 3^6 * 2


#########
### Ex 2:

R = c(1, 1)
d = c(-2, 1)
sol = solve.Shift2.S2P3(R, d)
x = sol[,1]; y = sol[,2];

### Test
(x + d[1])^3 + (y + d[1])^3 # - R[1]
(x + d[2])^3 + (y + d[2])^3 # - R[2]

### Classic Polynomial:
round0.p(poly.calc(x)) * 2


##########################
### Asymmetric Systems ###
##########################

### TODO:
# - move to separate file;

### Basic P2 System

# x^3 + b1*y^3 = R1
# x*y = R2

m3 = unity(3, all=T)

b = c(2)
R = c(1, 1)
#
det = sqrt(R[1]^2 - 4*b[1]*R[2]^3 + 0i)
x1 = m3 * ((R[1] + det)/2)^(1/3)
x2 = m3 * ((R[1] - det)/2)^(1/3)
x = c(x1, x2)
y = R[2]/x
sol = cbind(x, y)

### Test
x^3 + b[1]*y^3
x*y


###########################
###########################

###########################
### Derived Polynomials ###

### Note:
# - all strictly symmetric Polynomials P[2*n]
#   can be decomposed into P[2] o P[n];
# - these were some early experiments;


c2 = 2*cos(2*pi/5 * 1:2)
m = complex(re=cos(2*pi/3), im=sin(2*pi/3))
r3 = (1 + c(1i, -1i)*sqrt(3))/2

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


### non-symmetric
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



###
R1 = c(c2[1]-c2[2], 1)
R2 = c(c2[2]-c2[1], 1)
sol = rbind(
	solve.shift.ps(0, c2[2]-c2[1], R1, n=2),
	solve.shift.ps(0, c2[1]-c2[2], R2, n=2) )

poly.calc(sol$x)
x = sol$x
131 - 80*x + 50*x^2 + 20*x^3 - 33*x^4 + x^8


###
R1 = c(-1, 1)
R2 = c(-1, 1)
sol = rbind(
	solve.shift.ps(c2[2], c2[1]+1, R1, n=2),
	solve.shift.ps(c2[1], c2[2]+1, R2, n=2) )

poly.calc(sol$x)
x = sol$x
16 + 32*x^2 + 50*x^3 + 39*x^4 + 30*x^5 + 8*x^6 + x^8


### Experimental
R1 = c(c2[1], 1-c2[1])
R2 = c(c2[2], 1-c2[2])
sol = rbind(
	solve.shift.ps(1, -1+c2[1], R1, n=2),
	solve.shift.ps(1, -1+c2[2], R2, n=2) )

poly.calc(sol$x)
x = sol$x


################
###
R1 = c( 1+1i, 1)
R2 = c( 1-1i, 1)
sol = rbind(
	solve.shift.ps(-1-1i, -1i, R1, n=2),
	solve.shift.ps(-1+1i, +1i, R2, n=2) )

poly.calc(sol$x)
x = sol$x
1 + 8*x + 10*x^2 - 20*x^3 + 31*x^4 - 20*x^5 + 10*x^6 - 4*x^7 + x^8


###
R1 = c( 0, 1)
R2 = c( 0, 1)
sol = rbind(
	solve.shift.ps(-m, m^2, R1, n=2),
	solve.shift.ps(-m^2, m, R2, n=2) )

poly.calc(sol$x)
x = sol$x
1 - 2*x - 2*x^3 + 5*x^4 + 2*x^5 + 4*x^6 + x^8


###
R1 = c( 2*m-1, m^2)
R2 = c( 2*m^2-1, m)
sol = rbind(
	solve.shift.ps(m, 0, R1, n=2),
	solve.shift.ps(m^2, 0, R2, n=2) )

poly.calc(sol$x)
x = sol$x
1 + x - 3*x^2 - x^3 + 2*x^4 - x^5 + 3*x^6 - 2*x^7 + x^8


### (1 +/- i*sqrt(3))/2

###
R1 = c( 1, 2*r3[2])
R2 = c( 1, 2*r3[1])
sol = rbind(
	solve.shift.ps(-2, 1, R1, n=2),
	solve.shift.ps(-2, 1, R2, n=2) )

poly.calc(sol$x)
x = sol$x


