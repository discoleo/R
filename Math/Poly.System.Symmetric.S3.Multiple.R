########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Symmetric S3:
### Multiple/Composite Monoms
###
### draft v.0.1a


### Heterogenous Symmetric
### Polynomial Systems: 3 Variables
### Symmetric with Composite Monoms

### Example:
# (x*y)^n + (x*z)^n + (y*z)^n = R1

# Note:
# Heterosymmetric & Mixt Systems:
# - are discussed separately, see:
#   Poly.Heterosymmetric.S3.R &
#   Poly.Heterosymmetric.S3.Mixt.R


###############

###############
### History ###
###############


### draft v.0.1a:
# - solutions to the System order 2:
#   (x*y)^2 + (x*z)^2 + (y*z)^2 = R1;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

#############################

#############################
### Simple Multiple Terms ###

### (x*y)^n + (x*z)^n + (y*z)^n

### Solution:
# - decomposition into elementary polynomials;


##################
### Extensions ###

### Type A:
# A1: (x*y)^n + (x*z)^n + (y*z)^n + b1*(x+y+z) = R1
# A2: x*y + y*z + z*x + b2*(x+y+z) = R2
# A3: x*y*z + b3*(x+y+z) = R3

### Types Other:
# M1: (x*y + y*z + z*x)*(x+y+z) = R2
# M2: x*y*z*(x+y+z) = R3
# D1: (x*y + y*z + z*x)/(x+y+z) = R2
# D2: x*y*z/(x+y+z) = R3


################################

##############
### Simple ###
### n = 2  ###

### (x*y)^2 + (x*z)^2 + (y*z)^2 = R1
# [trivial: P3]

### Order 2: n = 2
(x*y)^2 + (x*z)^2 + (y*z)^2 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:
# [not run]
R2^2 - 2*R3*S - R1 # = 0
S = 1/2 * (R2^2 - R1) / R3

### Solver
solve.sym3 = function(R, b=0) {
	S = 1/2 * (R[2]^2 - R[1]) / R[3]
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=1, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	# TODO: robust
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 2, -1);
sol = solve.sym3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^2 + y^2*z^2 + z^2*x^2 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 1 + 2*x + 1.5*x^2 + x^3
round0(err)

##################
### Extensions ###

###########
### Type A:
# TODO


############
### Type M1:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0

### Solution:
# E2 = R2 / S
R2^2 / S^2 - 2*R3*S - R1 # = 0
2*R3*S^3 + R1*S^2 - R2^2 # = 0


#################
### Type M1 & M2:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z*(x+y+z) - R3 # = 0

### Solution:
# E2 = R2 / S
# E3 = R3 / S
R2^2 / S^2 - 2*R3 - R1 # = 0
(2*R3 + R1)*S^2 - R2^2 # = 0

#################
### Type M1 & D2:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z / (x+y+z) - R3 # = 0

### Solution:
# E2 = R2 / S
# E3 = R3 * S
R2^2 / S^2 - 2*R3*S^2 - R1 # = 0
2*R3*S^4 + R1*S^2 - R2^2 # = 0

#################
### Type D1 & D2:
(x*y + y*z + z*x) / (x+y+z) - R2 # = 0
x*y*z / (x+y+z) - R3 # = 0

### Solution:
# E2 = R2 * S
# E3 = R3 * S
(R2^2 - 2*R3)*S^2 - R1 # = 0


##########
### Solver
solve.sym3 = function(R, b=0) {
	# TODO: all types!
	coeff = c(2*R[3], R[1], 0, - R[2]^2)
	len = 3
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] / x - b2*x, -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	# TODO: robust
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 2, -1);
sol = solve.sym3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^2 + y^2*z^2 + z^2*x^2 # - R[1] # = 0
(x*y + y*z + z*x)*(x+y+z) # - R[2] # = 0
x*y*z # - R[3] # = 0

poly.calc(x)

err = 1 - 1.5*x^2 + 5*x^3 + 0.5*x^4 - 4*x^5 + 11*x^6 - 0.5*x^8 + x^9
round0(err)

