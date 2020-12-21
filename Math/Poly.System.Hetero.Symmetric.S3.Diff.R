
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Hetero-Symmetric Differences
###
### draft v.0.2c-ext


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

### draft v.0.2c - v.0.2c-ext:
# - variant system:
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


solve.Ht3DiffV1 = function(R, b) {
	if(R[1] == 0) stop("Currently NOT implemented: R[1] == 0!")
	if(length(b) == 1) {
		coeff = c(b[1], R[1]*b[1]^2, 0, 3*R[1])
	} else if(length(b) == 2) {
		coeff = c((b[1] - b[1]^2*b[2]), R[1]*b[1]^2, - 3*b[2], 3*R[1])
	} else if(length(b) == 3) {
		coeff = c(- b[1]^2*b[3], (b[1] - b[1]^2*b[2]), (R[1]*b[1]^2 - 3*b[3]), - 3*b[2], 3*R[1])
	}
	S = roots(coeff)
	print(S)
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

