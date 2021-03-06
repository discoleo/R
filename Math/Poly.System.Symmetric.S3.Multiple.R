########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Symmetric S3:
### Multiple/Composite Monoms
###
### draft v.0.1d


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


### draft v.0.1d:
# - Order 2 variant:
#   Eq 2: x^2 + y^2 + z^2 = R2;
### draft v.0.1c:
# - basic examples of M1 & D1 extensions;
### draft v.0.1b-name:
# - renamed extensions:
#   M1 & D1: [new] modified Eq 1;
#   M2 / D2, M3 / D3: the previous M1/D1 & M2/D2;
### draft v.0.1b:
# - some solutions to the system of order 3:
#   [with the non-trivial M2D3 extension (extension renamed)]
#   (x*y)^3 + (x*z)^3 + (y*z)^3 = R1;
### draft v.0.1a:
# - solutions to the system of order 2:
#   (x*y)^2 + (x*z)^2 + (y*z)^2 = R1;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

solve.S = function(S, R, b=0) {
	# generic solver (based on existing S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

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
# M1: ((x*y)^n + (x*z)^n + (y*z)^n)*(x+y+z) = R1
# M2: (x*y + y*z + z*x)*(x+y+z) = R2
# M3: x*y*z*(x+y+z) = R3
# D1: ((x*y)^n + (x*z)^n + (y*z)^n)/(x+y+z) = R1
# D2: (x*y + y*z + z*x)/(x+y+z) = R2
# D3: x*y*z/(x+y+z) = R3


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
S = 1/2 * (R[2]^2 - R[1]) / R[3]

### Solver
solve.sym3 = function(R, b=0) {
	S = 1/2 * (R[2]^2 - R[1]) / R[3]
	print(S)
	solve.S(S, R, b)
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
((x*y)^2 + (x*z)^2 + (y*z)^2) * (x+y+z) - R1 # = 0

### Solution:
(E2^2 - 2*E3*S)*S - R1 # = 0
2*E3*S^2 - E2^2*S + R1 # = 0
2*R[3]*S^2 - R[2]^2*S + R[1] # = 0

### Example:
R = c(1, 2, 1)
S = roots(c(2*R[3], - R[2]^2, R[1]))
sol = solve.S(S, R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
(x^2*y^2 + y^2*z^2 + z^2*x^2)*(x+y+z) # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


############
### Type D1:
((x*y)^2 + (x*z)^2 + (y*z)^2) / (x+y+z) - R1 # = 0

### Solution:
(E2^2 - 2*E3*S)/S - R1 # = 0
2*E3*S + R1*S - E2^2 # = 0
(2*R3+R1)*S - R2^2 # = 0
S = R[2]^2 / (2*R[3] + R[1])

### Example:
R = c(1, 2, 1)
S = R[2]^2 / (2*R[3] + R[1])
sol = solve.S(S, R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
(x^2*y^2 + y^2*z^2 + z^2*x^2)/(x+y+z) # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


############
### Type M2:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0

### Solution:
# E2 = R2 / S
R2^2 / S^2 - 2*R3*S - R1 # = 0
2*R[3]*S^3 + R[1]*S^2 - R[2]^2 # = 0


#################
### Type M2 & M3:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z*(x+y+z) - R3 # = 0

### Solution:
# E2 = R2 / S
# E3 = R3 / S
R2^2 / S^2 - 2*R3 - R1 # = 0
(2*R3 + R1)*S^2 - R2^2 # = 0


######################
### Type M1 & M2 & M3:
((x*y)^2 + (x*z)^2 + (y*z)^2) / (x+y+z) - R1 # = 0
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z*(x+y+z) - R3 # = 0

### Solution:
(2*R3 + R1/S)*S^2 - R2^2 # = 0
(2*R3*S + R1)*S - R2^2 # = 0
2*R[3]*S^2 + R[1]*S - R[2]^2 # = 0


#################
### Type M2 & D3:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z / (x+y+z) - R3 # = 0

### Solution:
# E2 = R2 / S
# E3 = R3 * S
R2^2 / S^2 - 2*R3*S^2 - R1 # = 0
2*R3*S^4 + R1*S^2 - R2^2 # = 0


#################
### Type D2 & D3:
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
	len = length(coeff) - 1
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] / x - b2, -R[3])))
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

################
### Variants ###

# (x*y)^2 + (x*z)^2 + (y*z)^2 = R1
# x^2 + y^2 + z^2 = R2
# x*y*z = R3

### Solution:
# [not run]

### Eq 2 =>
x^2 + y^2 + z^2 - R2 # = 0
S^2 - 2*E2 - R2 # = 0
# 2*E2 = S^2 - R2
# =>
4*E2^2 - 8*E3*S - 4*R1 # = 0
(S^2 - R2)^2 - 8*R3*S - 4*R1 # = 0
S^4 - 2*R2*S^2 - 8*R3*S + R2^2 - 4*R1 # = 0

### Solution:
solve.SymM2.S2MP2 = function(R, b=0) {
	coeff = c(1, 0, - 2*R[2], - 8*R[3], R[2]^2 - 4*R[1])
	# Extensions:
	b1 = if(length(b) > 0) b[1] else 0;
	b2 = if(length(b) > 1) b[2] else 0;
	b3 = if(length(b) > 2) b[3] else 0;
	if(b1 != 0) coeff = coeff + c(0, 0, 0, 4*b1, 0)
	if(b2 != 0) coeff = coeff + c(0, 2*b2, b2^2, -2*b2*R[2], 0)
	if(b3 != 0) coeff = coeff + c(0, 0, 8*b3, 0, 0)
	S = roots(coeff)
	#
	E2 = (S^2 - R[2] + b2*S) / 2
	E3 = R[3] - b3*S;
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	S  = matrix(S, ncol=len, nrow=3, byrow=T)
	E2 = matrix(E2, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	yz.s = S - x;
	yz = E2 - x*yz.s;
	yz.d = sqrt(yz.s^2 - 4*yz) # TODO: add also negative values;
	y = (yz.s + yz.d) / 2;
	z = (yz.s - yz.d) / 2;
	return(cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z)))
}

### Examples:

R = c(1,3,2)
sol = solve.SymM2.S2MP2(R)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test:
(x*y)^2 + (x*z)^2 + (y*z)^2 # - R1
x^2 + y^2 + z^2 # - R2
x*y*z # - R3

round0.p(poly.calc(x[1:3]^2))


### Ex 2:
R = c(0,2,1)
sol = solve.SymM2.S2MP2(R)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test:
(x*y)^2 + (x*z)^2 + (y*z)^2 # - R1
x^2 + y^2 + z^2 # - R2
x*y*z # - R3


### Extensions:
# Classical Polynomial: P12

### Ext A1:
R = c(1,3,2)
b = 1
sol = solve.SymM2.S2MP2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test:
(x*y)^2 + (x*z)^2 + (y*z)^2 + b[1]*(x+y+z) # - R1
x^2 + y^2 + z^2 # - R2
x*y*z # - R3

round0.p(poly.calc(x))


### Ext A2:
R = c(1,3,2)
b = c(1, -1)
sol = solve.SymM2.S2MP2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test:
(x*y)^2 + (x*z)^2 + (y*z)^2 + b[1]*(x+y+z) # - R1
x^2 + y^2 + z^2 + b[2]*(x+y+z) # - R2
x*y*z # - R3

round0.p(poly.calc(x))


### Ext A2:
R = c(1,3,2)
b = c(1, -1, -2)
sol = solve.SymM2.S2MP2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test:
(x*y)^2 + (x*z)^2 + (y*z)^2 + b[1]*(x+y+z) # - R1
x^2 + y^2 + z^2 + b[2]*(x+y+z) # - R2
x*y*z + b[3]*(x+y+z) # - R3

round0.p(poly.calc(x))


#######################
#######################

##############
### Simple ###
### n = 3  ###

### (x*y)^3 + (x*z)^3 + (y*z)^3 = R1
# [trivial: P3]

### Order 2: n = 2
(x*y)^3 + (x*z)^3 + (y*z)^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:
# [not run]
R2^3 - 3*R3*(R2*S - R3) - R1 # = 0
S = (R2^3 + 3*R3^2 - R1) / (3*R2*R3)


#################
### Type M2 & D3:
(x*y + y*z + z*x)*(x+y+z) - R2 # = 0
x*y*z / (x+y+z) - R3 # = 0

### Solution:
# E2 = R2 / S
# E3 = R3 * S
R2^3 / S^3 - 3*R3*S*(R2 - R3*S) - R1 # = 0
R2^3 - 3*R3*S^4*(R2 - R3*S) - R1*S^3 # = 0
3*R3^2*S^5 - 3*R2*R3*S^4 - R1*S^3 + R2^3 # = 0

##########
### Solver
solve.sym3 = function(R, b=0) {
	# TODO: all types!
	coeff = c(3*R[3]^2, - 3*R[2]*R[3], - R[1], 0, 0, R[2]^3)
	len = length(coeff) - 1
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] / x - b2, -R[3] * x)))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3] * S/x
	yz.s = S - x
	# TODO: robust
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}
test.sym3 = function(x, y, z, b=0, n=3, R, type="M2D3") {
	err1 = x^n*y^n + y^n*z^n + z^n*x^n # - R[1] # = 0
	err2 = (x*y + y*z + z*x) # - R[2] # = 0
	err3 = x*y*z # - R[3] # = 0
	S = x+y+z;
	if( ! is.na(pmatch("M2", type))) {
		err2 = err2 * S
	} else if( ! is.na(pmatch("D2", type))) {
		err2 = err2 / S
	}
	if( ! is.na(pmatch("M3", type))) {
		err3 = err3 * S
	} else if( ! is.na(pmatch("D3", type))) {
		err3 = err3 / S
	}
	# TODO: b, R;
	round0(rbind(err1, err2, err2))
}

### Examples:

### Ex 1:
R = c(3, 3, -1);
sol = solve.sym3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.sym3(x, y, z, n=3)

poly.calc(x)

err = -9 - 9*x^2 - 3*x^3 - 18*x^4 + 24*x^5 - 9*x^6 - 30*x^7 + 45*x^8 - 25*x^9 +
	- 15*x^10 + 38*x^11 + 12*x^12 - x^13 + 3*x^14 + x^15
round0(err)


### Ex 2:
R = c(0, 3, -1);
sol = solve.sym3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.sym3(x, y, z, n=3)

poly.calc(x)

err = -9 - 9*x^2 - 9*x^4 + 18*x^5 + 9*x^6 - 27*x^7 + 9*x^8 - 27*x^9 +
	- 6*x^10 + 36*x^11 + 12*x^12 + 3*x^14 + x^15
round0(err)


