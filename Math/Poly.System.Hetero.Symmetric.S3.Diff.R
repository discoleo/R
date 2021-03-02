
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Hetero-Symmetric Differences
###
### draft v.0.3a-class


### Hetero-Symmetric Differences
### Polynomial Systems: 3 Variables

### Example:
x^n - y^n + b*x*y = R
y^n - z^n + b*y*z = R
z^n - x^n + b*z*x = R

# - for full derivation of solutions, see:
#   Poly.System.Hetero.Symmetric.S3.Diff.Derivation.R;

#####################

###############
### History ###

### draft v.0.3a - v.0.3a-class:
# - started work on + classic Polynomial (P[24]):
#   x^3 + y^3 - z^3 + b*x = R;
### draft v.0.2d:
# - variant system: Multiplicative
#   x^2 - y^2 + b*x*y*(x*y*z) = R;
### draft v.0.2c - v.0.2c-ext:
# - variant system: "Additive"
#   x^2 - y^2 + b*x*y*(x+y+z) = R;
# - including A1 type extensions (powers 1 & 2); [v.0.2c-ext]
### draft v.0.2b:
# - implemented various power extensions of type A1:
#   powers 1, 2 & 3: system P3 + b[k+1]*(x+y+z)^k;
### draft v.0.2a:
# - moved derivation of formulas to:
#   Poly.System.Hetero.Symmetric.S3.Diff.Derivation.R;
### draft v.0.1b:
# - solved: x^3 - y^3 + b1*x*y = R;
### draft v.0.1a-ext2:
# - extension A1 power 2: + b[3]*(x+y+z)^2;
### draft v.0.1a:
# - moved [Difference] section to this new file;
### old file: Poly.System.Hetero.Symmetric.S3.R
### draft v.0.2b - v.0.2b-ext:
# - first concepts / solved [v.0.2b-sol]:
#   x^2 - y^2 + b*x*y = R;
# - solved extension A1 (Pow 1): + b[2]*(x+y+z); [v.0.2b-ext]


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

##########################


##########################
### Difference Types   ###
##########################

# quasi-[Negative Correlations]


###############
### Order 2 ###
###############

# x^2 - y^2 + b*x*y = R
# y^2 - z^2 + b*y*z = R
# z^2 - x^2 + b*x*z = R

### Solution:

### Case 1:
det = sqrt(R[1]/b[1])
x = y = z = c(det, -det)

### Case 2:
# (x, y, z) NOT equal;

### Sum =>
b*E2 - 3*R # = 0

### Sum(z^2*...) =>
# b*x*y*z*S = R*(x^2+y^2+z^2)
b*E3*S - R*(S^2 - 2*E2) # = 0

### Sum(x*y*...) =>
b^3*S^6 + R*(b^4 - 15*b^2)*S^4 + R^2*(27*b - 18*b^3)*S^2 + 81*R^3*(b^2 + 3)
(b*S^2 - 9*R)^2 * (b*S^2 + R*(b^2 + 3))

### Extension A1: pow 1:
(b[1]*S^2 + 9*b[2]*S - 9*R)^2 * (b[1]*S^2 - b[2]*(b[1]^2 + 3)*S + R*(b[1]^2 + 3))
### Extension A1: pow 2:
((b[1] + 9*b[3])*S^2 + 9*b[2]*S - 9*R)^2 * ((b[1] - b[3]*(b[1]^2 + 3))*S^2 - b[2]*(b[1]^2 + 3)*S + R*(b[1]^2 + 3))


### Solution

solve.Ht3Diff = function(R, b) {
	if(length(b) == 1) {
		coeff = c(b[1], 0, R[1]*(b[1]^2 + 3))
	} else if(length(b) == 2) {
		coeff = c(b[1], - b[2]*(b[1]^2 + 3), R[1]*(b[1]^2 + 3))
	} else if(length(b) == 3) {
		coeff = c((b[1] - b[3]*(b[1]^2 + 3)), - b[2]*(b[1]^2 + 3), R[1]*(b[1]^2 + 3))
	}
	det = sqrt(coeff[2]^2 - 4*coeff[1]*coeff[3] + 0i)
	S = c( - coeff[2] + det, - coeff[2] - det) / 2 / coeff[1]
	print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3: pow 2;
	R1 = R[1] - b2*S - b3*S^2
	E2 = 3*R1 / b[1]
	E3 = R1*(S^2 - 2*E2) / b[1] / S
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y2 = (yz.s^2 + R1 - (2 + b[1])*yz) / 2
	y = (R1 - x^2 + y2) / b[1] / x
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	### x = 0
	if(any(isZero)) {
		print("Solution: x == 0")
		# cleanup
		sol = sol[ ! isZero , ];
		R1 = R1[isZero]; yz.s = yz.s[isZero];
		# y2 = - R1; z2 = R1;
		yz = 3*R1 / b[1];
		yz.d = (R1 - b[1]*yz) / yz.s;
		y = (yz.s + yz.d) / 2; z = yz.s - y;
		sol = rbind(0, unique(cbind(x,y,z)))
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = c(1,-1) * sqrt(R[1] / b[1] + 0i)
	} else if(length(b) == 2) {
		x = y = z = roots(c(b[1], 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(b[1] + 9*b[3], 3*b[2], -R[1]))
	}
	sol = rbind(sol, cbind(x,y,z))
	return(sol)
}

### Examples:

R = 1
b = 1
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y # - R
y^2 - z^2 + b[1]*y*z # - R
z^2 - x^2 + b[1]*x*z # - R

round0.p(poly.calc(x[1:6]))

#########
### Ex 2:
R = 1
b = 5
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y # - R
y^2 - z^2 + b[1]*y*z # - R
z^2 - x^2 + b[1]*x*z # - R

round0.p(poly.calc(x[1:6])) * 5

###############
### Extensions:

R = 1
b = c(1, -2)
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y + b[2]*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z + b[2]*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z + b[2]*(x+y+z) # - R

round0.p(poly.calc(x[1:6]))

err = 25 + 60*x - 131*x^2 - 284*x^3 - 38*x^4 + 8*x^5 + x^6
round0(err)

##########
### Ext 2:

R = 1
b = c(1, -2, 1)
#
sol = solve.Ht3Diff(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z + b[3]*(x+y+z)^2 + b[2]*(x+y+z) # - R

round0.p(poly.calc(x[1:6]))

err = 25 + 60*x - 131*x^2 - 284*x^3 - 38*x^4 + 8*x^5 + x^6
round0(err)


###########################
###########################

###############
### Order 3 ###
###############

# x^3 - y^3 + b*x*y = R
# y^3 - z^3 + b*y*z = R
# z^3 - x^3 + b*x*z = R

### Solution:

# - for full derivation, see:
#   Poly.System.Hetero.Symmetric.S3.Diff.Derivation.R;

### Case 1:
det = sqrt(R[1]/b[1])
x = y = z = c(det, -det)


### Case 2:
# (x, y, z) NOT equal;

### Sum =>
b*E2 - 3*R # = 0

### Sum(z^3*...) =>
# b*x*y*z*(x^2 + y^2 + z^2) = R*(x^3 + y^3 + z^3)
b*(b*S^2 - 9*R)*E3 - R*(b*S^3 - 9*R*S) # = 0

### Eq 3:
### sq: x^3 - y^3 = R - b[1]*x*y =>
x^6 + y^6 - 2*x^3*y^3 - (R^2 + b[1]^2*x^2*y^2 - 2*b[1]*R*x*y) # = 0
### Sum(...) =>
(27*R^3*b1^4 + 729*R^4*b1) +
	+ (- 12*R^2*b1^5 - 567*R^3*b1^2)*S^2 +
	+ (R*b1^6 + 162*R^2*b1^3)*S^4 +
	+ (- 21*R*b1^4)*S^6 + b1^5*S^8
(b[1]*S^2 - 9*R) * (b[1]^4*S^6 - 12*R*b[1]^3*S^4 + (R*b[1]^5 + 54*R^2*b[1]^2)*S^2 - (3*R^2*b[1]^4 + 81*R^3*b[1]))
(b[1]*S^2 - 9*R) * (b[1]*S^2 - 3*R) * (b[1]^3*S^4 - 9*R*b[1]^2*S^2 + (R*b[1]^4 + 27*R^2*b[1]))


coeff = c(b[1]^3, 0, - 9*R[1]*b[1]^2, 0, (R[1]*b[1]^4 + 27*R[1]^2*b[1]))


### Extension A1: power 1:
(b[1]*S^2 + 9*b[2]*S - 9*R) *
	(b[1]^3*S^4 + 9*b[1]^2*b[2]*S^3 +
	+ (27*b[1]*b[2]^2 - 9*R*b[1]^2)*S^2  - (b[1]^4*b[2] + 54*b[1]*b[2]*R)*S + (R*b[1]^4 + 27*R^2*b[1]))
### Extension A1: power 2:
((b[1] + 9*b[3])*S^2 + 9*b[2]*S - 9*R) *
	((b[1]^3 + 9*b[1]^2*b[3] + 27*b[1]*b[3]^2)*S^4 + 9*(b[1]^2*b[2] + 6*b[1]*b[2]*b[3])*S^3 +
	+ (27*b[1]*b[2]^2 - 9*R*b[1]^2 + 54*b[1]*b[2]*b[3] - b[1]^4*b[3] + 54*b[1]*b[3]*R)*S^2  +
	- (b[1]^4*b[2] + 54*b[1]*b[2]*R)*S + (R*b[1]^4 + 27*b[1]*R^2))
### Extension A1: power 3:
(9*b[4]*S^3 + (b[1] + 9*b[3])*S^2 + 9*b[2]*S - 9*R) *
	(27*b[1]*b[4]^2*S^6 + 9*(b[1]^2*b[4] + 6*b[1]*b[3]*b[4])*S^5 +
	+ (b[1]^3 + 9*b[1]^2*b[3] + 27*b[1]*b[3]^2 + 54*b[1]*b[2]*b[4])*S^4 +
	+ (9*b[1]^2*b[2] + 54*b[1]*b[2]*b[3] - b[1]^4*b[4] - 54*b[1]*b[4]*R)*S^3 +
	+ (27*b[1]*b[2]^2 - 9*R*b[1]^2 - b[1]^4*b[3] - 54*b[1]*b[3]*R)*S^2  +
	- (b[1]^4*b[2] + 54*b[1]*b[2]*R)*S + (R*b[1]^4 + 27*b[1]*R^2))


### Solution
solve.Ht3DiffP3 = function(R, b) {
	if(b[1] == 0) stop("System undefined when b[1] == 0!")
	if(length(b) == 1) {
		coeff = c(b[1]^3, 0, - 9*R[1]*b[1]^2, 0, (R[1]*b[1]^4 + 27*R[1]^2*b[1]))
	} else if(length(b) == 2) {
		coeff = c(b[1]^3, 9*b[1]^2*b[2], (27*b[1]*b[2]^2 - 9*R[1]*b[1]^2),
			- (b[1]^4*b[2] + 54*b[1]*b[2]*R[1]), (R[1]*b[1]^4 + 27*R[1]^2*b[1]))
	} else if(length(b) == 3) {
		coeff = c((b[1]^3 + 9*b[1]^2*b[3] + 27*b[1]*b[3]^2), 9*(b[1]^2*b[2] + 6*b[1]*b[2]*b[3]),
			(27*b[1]*b[2]^2 - 9*R[1]*b[1]^2 + 54*b[1]*b[2]*b[3] - b[1]^4*b[3] + 54*b[1]*b[3]*R[1]),
			- (b[1]^4*b[2] + 54*b[1]*b[2]*R[1]), (R*b[1]^4 + 27*R[1]^2*b[1]))
	} else if(length(b) == 4) {
		coeff = c(27*b[1]*b[4]^2, 9*(b[1]^2*b[4] + 6*b[1]*b[3]*b[4]),
			(b[1]^3 + 9*b[1]^2*b[3] + 27*b[1]*b[3]^2 + 54*b[1]*b[2]*b[4]), # S^4
			9*(b[1]^2*b[2] + 6*b[1]*b[2]*b[3] - 6*b[1]*b[4]*R[1]) - b[1]^4*b[4], # S^3
			(27*b[1]*b[2]^2 - 9*R[1]*b[1]^2 - b[1]^4*b[3] - 54*b[1]*b[3]*R[1]), # S^2
			- (b[1]^4*b[2] + 54*b[1]*b[2]*R[1]), (R[1]*b[1]^4 + 27*b[1]*R[1]^2))
	}
	S = roots(coeff)
	print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A1: pow 2;
	b4 = if(length(b) > 3) b[4] else 0; # Ext A1: pow 3;
	R1 = R[1] - b2*S - b3*S^2 - b4*S^3
	E2 = 3*R1 / b[1]
	E3 = R1*(b[1]*S^3 - 9*R1*S) / b[1] / (b[1]*S^2 - 9*R1)
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	# x = x[ ! isZero] # TODO: changes length of x => S;
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y3 = (yz.s^3 + R1 - (3*yz.s + b[1])*yz) / 2
	y = (R1 - x^3 + y3) / b[1] / x
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	### x = 0
	if(any(isZero)) {
		print("Solution: x == 0")
		# cleanup
		sol = sol[ ! isZero , ];
		R1 = R1[isZero]; yz.s = yz.s[isZero];
		# y3 = - R1; z3 = R1;
		yz = 3*R1 / b[1];
		yz.d = - 2*R1 / (yz.s^2 - yz)
		y = (yz.s + yz.d)/2;
		z = yz.s - y;
		sol2 = cbind(0, y, z)
		sol2 = unique(sol2) # remove 3x duplicates;
		sol = rbind(sol, sol2)
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = c(1,-1) * sqrt(R[1] / b[1] + 0i)
	} else if(length(b) == 2) {
		x = y = z = roots(c(b[1], 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(b[1] + 9*b[3], 3*b[2], -R[1]))
	} else if(length(b) == 4) {
		x = y = z = roots(c(27*b[4], (b[1] + 9*b[3]), 3*b[2], - R[1]))
	}
	sol = rbind(sol, cbind(x,y,z))
	return(sol)
}


### Examples:

R = 2
b = 3
#
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y # - R[1]
y^3 - z^3 + b[1]*y*z # - R[1]
z^3 - x^3 + b[1]*x*z # - R[1]


#########
### Ex 2: x = 0
R = -1; b = 3;
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y # - R[1]
y^3 - z^3 + b[1]*y*z # - R[1]
z^3 - x^3 + b[1]*x*z # - R[1]


###############
### Extensions:

### Ext. power 1:
R = 1
b = c(1, -2)
#
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y + b[2]*(x+y+z) # - R
y^3 - z^3 + b[1]*y*z + b[2]*(x+y+z) # - R
z^3 - x^3 + b[1]*x*z + b[2]*(x+y+z) # - R

round0.p(poly.calc(x[1:12]))

err = 28 + 330*x + 1563*x^2 + 4526*x^3 + 20325*x^4 + 47980*x^5 + 34722*x^6 - 20524*x^7 + 7654*x^8 +
	- 1510*x^9 + 219*x^10 - 18*x^11 + x^12
round0(err)


### Ext. power 2:
R = 1
b = c(1, -2, -2)
#
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R
y^3 - z^3 + b[1]*y*z + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R
z^3 - x^3 + b[1]*x*z + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R

round0.p(poly.calc(x[1:12])) ) * 7^3 * 13^3

err = 28 + 330*x + 2217*x^2 + 10206*x^3 + 34283*x^4 + 83236*x^5 + 147578*x^6 + 217756*x^7 + 410766*x^8 +
	+ 920194*x^9 + 1572389*x^10 + 1639638*x^11 + 753571*x^12
round0(err)


### Ext. power 3:
R = 1
b = c(1, -2, -2, 1)
#
sol = solve.Ht3DiffP3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 - y^3 + b[1]*x*y + b[2]*(x+y+z) + b[3]*(x+y+z)^2 + b[4]*(x+y+z)^3 # - R
y^3 - z^3 + b[1]*y*z + b[2]*(x+y+z) + b[3]*(x+y+z)^2 + b[4]*(x+y+z)^3 # - R
z^3 - x^3 + b[1]*x*z + b[2]*(x+y+z) + b[3]*(x+y+z)^2 + b[4]*(x+y+z)^3 # - R


##########################

##########################
### Difference Types   ###
### Variants           ###
##########################

### x^2 - y^2 + b*x*y*(x+y+z)


###############
### Order 2 ###
###############

# x^2 - y^2 + b*x*y*(x+y+z) = R
# y^2 - z^2 + b*y*z*(x+y+z) = R
# z^2 - x^2 + b*x*z*(x+y+z) = R


### Solution:

### E2:
# b*E2*S = 3*R
### E3:
# b*E3*S^2 = R*(S^2 - 2*E2)

### Eq:
(b[1]*S^3 - 9*R[1]) * (b[1]*S^3 + R[1]*b[1]^2*S^2 + 3*R[1])

### Extension A1: power 1:
(b[1]*S^3 + 9*b[2]*S - 9*R[1]) *
	((b[1] - b[1]^2*b[2])*S^3 + R[1]*b[1]^2*S^2 - 3*b[2]*S + 3*R[1])
### Extension A1: power 2:
(b[1]*S^3 + 9*b[3]*S^2 + 9*b[2]*S - 9*R[1]) *
	(- b[1]^2*b[3]*S^4 + (b[1] - b[1]^2*b[2])*S^3 + (R[1]*b[1]^2 - 3*b[3])*S^2 - 3*b[2]*S + 3*R[1])


solve.Ht3DiffV1 = function(R, b, debug=TRUE) {
	if(R[1] == 0) stop("Currently NOT implemented: R[1] == 0!")
	if(length(b) == 1) {
		coeff = c(b[1], R[1]*b[1]^2, 0, 3*R[1])
	} else if(length(b) == 2) {
		coeff = c((b[1] - b[1]^2*b[2]), R[1]*b[1]^2, - 3*b[2], 3*R[1])
	} else if(length(b) == 3) {
		coeff = c(- b[1]^2*b[3], (b[1] - b[1]^2*b[2]), (R[1]*b[1]^2 - 3*b[3]), - 3*b[2], 3*R[1])
	}
	S = roots(coeff)
	if(debug) print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3: pow 2;
	R1 = R[1] - b2*S - b3*S^2
	E2 = 3*R1 / b[1] / S;
	E3 = R1*(S^2 - 2*E2) / b[1] / S^2;
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y2 = (yz.s^2 + R1 - 2*yz - b[1]*yz*S) / 2
	y = (R1 - x^2 + y2) / b[1] / x / S
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	### x = 0
	if(any(isZero)) {
		print("Solution: x == 0")
		# cleanup
		sol = sol[ ! isZero , ];
		R1 = R1[isZero]; yz.s = yz.s[isZero];
		# y2 = - R1; z2 = R1;
		# TODO: check!
		yz = 3*R1 / b[1] / yz.s;
		yz.d = (R1 - b[1]*yz^2) / yz.s;
		y = (yz.s + yz.d) / 2; z = yz.s - y;
		sol = rbind(0, unique(cbind(x,y,z)))
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = roots(c(3*b[1], 0, 0, -R[1]))
	} else if(length(b) == 2) {
		x = y = z = roots(c(3*b[1], 0, 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(3*b[1], 9*b[3], 3*b[2], -R[1]))
	}
	sol = rbind(sol, cbind(x,y,z))
	return(sol)
}

### Examples:

R = 1
b = 1
#
sol = solve.Ht3DiffV1(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z*(x+y+z) # - R


### Ext A1: power 1

R = 1
b = c(-1, -1)
#
sol = solve.Ht3DiffV1(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y*(x+y+z) + b[2]*(x+y+z) # - R
y^2 - z^2 + b[1]*y*z*(x+y+z) + b[2]*(x+y+z) # - R
z^2 - x^2 + b[1]*x*z*(x+y+z) + b[2]*(x+y+z) # - R


### Ext A1: power 2

R = 1
b = c(-1, -1, 2)
#
sol = solve.Ht3DiffV1(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y*(x+y+z) + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R
y^2 - z^2 + b[1]*y*z*(x+y+z) + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R
z^2 - x^2 + b[1]*x*z*(x+y+z) + b[2]*(x+y+z) + b[3]*(x+y+z)^2 # - R


#############################
#############################

### x^2 - y^2 + b*x*y*(x*y*z)


###############
### Order 2 ###
###############

# x^2 - y^2 + b*x*y*(x*y*z) = R
# y^2 - z^2 + b*y*z*(x*y*z) = R
# z^2 - x^2 + b*x*z*(x*y*z) = R


### Solution:

### E2:
# b*E2*E3 = 3*R
### E3:
# E3 =
#    (6*R^2*S^8*b^4 + 162*R^3*S^3*b^3 - 99*R^4*S^4*b^4 - 6*R^5*S^5*b^5 + 486*R^6*b^4) /
#    (R*S^10*b^5 - 63*R^2*S^5*b^4 - 12*R^3*S^6*b^5 + 243*R^4*S*b^4 - R^4*S^7*b^6 + 135*R^5*S^2*b^5);

### Eq:
(b[1]*S^5 - 243*R[1]) *
	(b[1]*S^7 - b[1]^2*R[1]^3*S^4 - 8*b[1]*R[1]^2*S^3 - 3*R[1]*S^2 - 4*b[1]^2*R[1]^5)

coeff = c(R*b^3, - 12*b^2, 0, - 2*R^4*b^4, 4*R^3*b^3, 255*R^2*b^2, (252*R*b + R^7*b^5), # S^14
	4*R^6*b^4, - 30*R^5*b^3, - 3456*R^4*b^2, (- 2214*R^3*b + 4*R^9*b^5), # S^10
	(- 648*R^2 - 81*R^8*b^4), 648*R^7*b^3, 18819*R^6*b^2, 5832*R^5*b, - 324*R^10*b^4,
	3402*R^9*b^3, - 34992*R^8*b^2, - 13122*R^7*b, 0, - 17496*R^11*b^3)

### TODO:
### Extension A1: power 1:
(b[1]*S^3 + 9*b[2]*S - 9*R[1]) *
	()
### Extension A1: power 2:
(b[1]*S^3 + 9*b[3]*S^2 + 9*b[2]*S - 9*R[1]) *
	()


solve.Ht3DiffV2 = function(R, b, debug=TRUE) {
	if(R[1] == 0) stop("Currently NOT implemented: R[1] == 0!")
	if(length(b) == 1) {
		coeff = c(b[1], 0, 0, - b[1]^2*R[1]^3, - 8*b[1]*R[1]^2, - 3*R[1], 0, - 4*b[1]^2*R[1]^5)
	} else if(length(b) == 2) {
		#
	} else if(length(b) == 3) {
		#
	}
	S = roots(coeff)
	if(debug) print(S)
	#
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: pow 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3: pow 2;
	R1 = R[1] - b2*S - b3*S^2
	E3 = (6*R^2*S^8*b^4 + 162*R^3*S^3*b^3 - 99*R^4*S^4*b^4 - 6*R^5*S^5*b^5 + 486*R^6*b^4) /
		(R*S^10*b^5 - 63*R^2*S^5*b^4 - 12*R^3*S^6*b^5 + 243*R^4*S*b^4 - R^4*S^7*b^6 + 135*R^5*S^2*b^5);
	E2 = 3*R1 / b[1] / E3;
	### x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], - E3[id])))
	len = length(S)
	S  = rep(S,  each=3)
	E3 = rep(E3, each=3)
	R1 = rep(R1, each=3)
	isZero = round0(x) == 0
	yz = E3/x
	yz.s = S - x
	### robust: includes Ext A1: powers 1 & 2;
	y2 = (yz.s^2 + R1 - 2*yz - b[1]*yz*E3) / 2
	y = (R1 - x^2 + y2) / b[1] / x / E3
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z), S=S)
	### x = 0
	if(any(isZero)) {
		print("Solution: x == 0")
		# cleanup
		sol = sol[ ! isZero , ];
		R1 = R1[isZero]; yz.s = yz.s[isZero];
		# y2 = - R1; z2 = R1;
		# TODO: check!
		yz = yz.s^2 / 2;
		yz.d = -2*R1 / yz.s;
		y = (yz.s + yz.d) / 2; z = yz.s - y;
		sol = rbind(0, unique(cbind(x,y,z, S[isZero])))
	}
	### x == y == z
	if(length(b) < 2) {
		x = y = z = roots(c(b[1], 0, 0, 0, 0, -R[1]))
	} else if(length(b) == 2) {
		x = y = z = roots(c(b[1], 0, 0, 0, 3*b[2], -R[1]))
	} else if(length(b) == 3) {
		x = y = z = roots(c(b[1], 0, 0, 9*b[3], 3*b[2], -R[1]))
	}
	sol = rbind(sol, cbind(as.vector(x), as.vector(y), as.vector(z), 3*as.vector(x)))
	return(sol)
}

### Examples:

R = 1
b = 1
#
sol = solve.Ht3DiffV2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y*(x*y*z) # - R
y^2 - z^2 + b[1]*y*z*(x*y*z) # - R
z^2 - x^2 + b[1]*x*z*(x*y*z) # - R

#########
### Ex 2:
R = 1
b = -2
#
sol = solve.Ht3DiffV2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 - y^2 + b[1]*x*y*(x*y*z) # - R
y^2 - z^2 + b[1]*y*z*(x*y*z) # - R
z^2 - x^2 + b[1]*x*z*(x*y*z) # - R


###########################
###########################

#######################
### Leading Terms:  ###
### Complete        ###
#######################

###############
### Order 3 ###
###############

#   x^3 + y^3 - z^3 + b*x = R
# - x^3 + y^3 + z^3 + b*y = R
#   x^3 - y^3 + z^3 + b*z = R

### Equivalent system:
2*x^3 + b*(x+z) - 2*R # = 0
2*y^3 + b*(x+y) - 2*R # = 0
2*z^3 + b*(y+z) - 2*R # = 0

### Solution:

### Sum =>
x^3 + y^3 + z^3 + b*(x+y+z) - 3*R # = 0
S^3 - 3*E2*S + 3*E3 + b*S - 3*R

### Sum(x^3*...) =>
(x^6 + y^6 + z^6) + b*(x^4 + y^4 + z^4) - R*(x^3 + y^3 + z^3) # = 0
(- 2*E2^3 + 3*E3^2 - 12*E2*E3*S + 9*E2^2*S^2 + 6*E3*S^3 - 6*E2*S^4 + S^6) +
	+ b*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) - R*(S^3 - 3*E2*S + 3*E3) # = 0

### Sum(x*...) & Rotation =>
# [see Derivation]
(4*E2^2*b^2 - 8*E2^3*b + 4*E2^4) +
(- 16*E2*E3*b + 4*E2*R*b + 16*E2^2*E3 - 4*E2^2*R)*S^1 +
(- 4*E2*b^2 + 20*E2^2*b - 12*E2^3 - 8*E3*R + 43*E3^2 + R^2)*S^2 +
(- 50*E2*E3 + 8*E2*R + 8*E3*b - 2*R*b)*S^3 +
(- 12*E2*b + 19*E2^2 + b^2)*S^4 +
(12*E3 - 2*R)*S^5 +
(- 8*E2 + 2*b)*S^6 + S^8

### E2:
E2Subst = 360*R*S^6*b^3 + 1098*R*S^8*b^2 - 306*R*S^10*b - 1152*R*S^12 + 243*R^2*S*b^4 +
	+ 1458*R^2*S^3*b^3 + 4185*R^2*S^5*b^2 - 2862*R^2*S^7*b - 3024*R^2*S^9 - 1458*R^3*S^2*b^2 +
	- 6318*R^3*S^4*b + 7776*R^3*S^6 - 27*S^3*b^6 - 138*S^5*b^5 - 432*S^7*b^4 - 468*S^9*b^3 +
	+ 219*S^11*b^2 + 606*S^13*b + 240*S^15;
E2Div = - 162*R*S^2*b^4 - 882*R*S^4*b^3 + 396*R*S^6*b^2 + 72*R*S^8*b + 576*R*S^10 +
	- 243*R^2*S*b^3 - 729*R^2*S^3*b^2 - 3348*R^2*S^5*b + 4320*R^2*S^7 + 153*S^3*b^5 +
	+ 575*S^5*b^4 + 571*S^7*b^3 + 33*S^9*b^2 - 724*S^11*b - 608*S^13;
- E2Subst / E2Div;


### TODO: solve!
# - massive overflow!
# - S-Polynomial: actually P[8];


### Test
  x^3 + y^3 - z^3 + b*x # - R
- x^3 + y^3 + z^3 + b*y # - R
  x^3 - y^3 + z^3 + b*z # - R

### Classic Polynomial: P[3] * P[24]
# see file:
# Poly.System.Hetero.Symmetric.S3.Diff.Derivation.R;


### Debug
b = 3; R = 2;
x =  0.2917949139 + 0.4133125124i;
y = -0.3108981004 + 1.4856874177i;
z =  1.1246683294 - 0.4366248387i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;

