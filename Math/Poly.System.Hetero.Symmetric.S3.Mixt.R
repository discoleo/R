########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric S3:
### Mixt Type
###
### draft v.0.2g


### Heterogenous Symmetric
### Polynomial Systems: 3 Variables
### Mixt: Hetero + Symmetric

### Example:
# x^p*y^n + y^p*z^n + z^p*x^n = R1
# x*y + x*z + y*z = R2
# x*y*z = R3


###############

###############
### History ###
###############

### draft v.0.2g:
# - classic Polynomial for the simple Dual system:
#   degenerate P18: pseudo-P6;
### draft v.0.2f:
# - linear extension (type A3) to the Order 2 system;
### draft v.0.2e:
# - combined variant:
#   (x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2) = R1;
### draft v.0.2d:
# - Extensions to the Order 2 system:
#   M3 extension; [v.0.2d]
### draft v.0.2c - v.0.2c-robust:
# - added system with 2 hetero-symmetric equations:
#   Order 2: trivial polynomials;
#   Order 3: intermediary polynomial of order 13; [v.0.2c-ord3]
# - [DONE]: fixed robust roots! [v.0.2c-robust]
### draft v.0.2a - v.0.2b-sp: [08-12-2020]
# - Generalization:
#   x^p*y^n + y^p*z^n + z^p*x^n = R1;
# - solved: n = 3, p = 2; [full in v.0.2b]
# - [TODO][DONE]: robust special case; [v.0.2b-sp]
### draft v.0.1c: [08-12-2020]
# - full robust solution: using dS / dR1;
### draft v.0.1b: [07-12-2020]
# - solved + robust:
#   x*y^3 + y*z^3 + z*x^3 = R1;
### draft v.0.1a-ext1 - v.0.1a-ext2:
# - Extensions:
#   x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1;
#   x*y + x*z + y*z + b2*(x+y+z) = R2;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Hetero.Symmetric.S3.R;
#
#########################
### former branch v.0.1c:
### in Poly.System.Hetero.Symmetric.S3.R
### draft v.0.1c-pre-alpha - v.0.1c-exact:
# - first look & solved + exact/robust solution: [v.0.1c-exact]
#   x*y^2 + y*z^2 + z*x^2 = R1;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions

test.ht3 = function(x, y, z, R, n=2, p=1, b=0) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	### Test
	x.sum = x + y + z
	err1 = x^p*y^n + y^p*z^n + z^p*x^n + b[1]*x.sum
	err2 = x*y + y*z + z*x + if(length(b) < 2) 0 else b[2]*x.sum
	err3 = x*y*z + if(length(b) < 3) 0 else b[3]*x.sum
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		err = err - rep(R, each=length(x))
	}
	err = round0(err)
	return(err)
}
test.ht3Dual = function(x, y, z, R, n=2, p=1, b=0, type) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	x.sum = x + y + z
	e3f = if(missing(type)) 1 else if(match("M3", type) > 0) x.sum else if(match("D3", type) > 0) 1/x.sum;
	if( ! missing(type)) print(paste0("Type: ", type))
	### Test
	err1 = x^p*y^n + y^p*z^n + z^p*x^n + b[1]*x.sum
	err2 = x^p*z^n + y^p*x^n + z^p*y^n + if(length(b) < 2) 0 else b[2]*x.sum
	err3 = x*y*z * e3f + if(length(b) < 3) 0 else b[3]*x.sum
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		err = err - rep(R, each=length(x))
	}
	err = round0(err)
	return(err)
}

############################

### Mixt Systems Type 1 ###

### Generalization:
### x^p*y^n + y^p*z^n + z^p*x^n
### x*y + x*z + y*z
### x*y*z

################################

##############
### Simple ###
### p = 1  ###

### x*y^n + y*z^n + z*x^n = R1

### Order 2: n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Solution:

### Eq 1: * x
x^2*y^2 + x*y*z^2 + z*x^3 - R1*x # = 0
R3^2/z^2 + R3*z + z*x^3 - R1*x # = 0 # *z^2
R3^2 + R3*z^3 + x^3*z^3 - R1*x*z^2 # = 0
# similar:
R3^2 + R3*x^3 + x^3*y^3 - R1*y*x^2 # = 0
R3^2 + R3*y^3 + y^3*z^3 - R1*z*y^2 # = 0

### Sum =>
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) - R1*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*(x^2*y + x^2*z + x*y^2 + y^2*z + x*z^2 + y*z^2) + R1^2 + 3*R3^2 # = 0
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) + R1^2 + 3*R3^2 # = 0
(R1+R3)*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- R1*((S^2 - 2*R2)*S) + R1^2 + 3*R3^2 # = 0
(R1+R3)*(S^3 - 3*R2*S + 3*R3) + (R2^3 - 3*R3*(R2*S - R3)) +
	- R1*((S^2 - 2*R2)*S) + R1^2 + 3*R3^2 # = 0
R3*S^3 - (R1+6*R3)*R2*S + R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0

### Solution
solve.ht3 = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3], 0, - (R[1]+6*R[3])*R[2], R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	} else if(length(b) < 3) {
		coeff = c(R[3], (b[1]^2 + b[1]*R[2]), - (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
		if(length(b) > 1) {
			# Ext 2:
			coeff = coeff - c(b[2]^3 + b[1]*b[2], -b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0)
		}
	} else {
		# Ext 3:
		coeff = c(-b[3], R[3], (b[1]^2 + 3*b[1]*b[3] + 9*b[3]^2 + b[1]*R[2] + 6*b[3]*R[2]),
			- (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1] + 3*b[3]*R[1] + 18*b[3]*R[3]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
		coeff = coeff - c(0, b[2]^3 + b[1]*b[2] + 6*b[2]*b[3],
			- b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0) # Ext 2
	}
	S = roots(coeff)
	len = length(S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext 3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = (R[3] - b3*S) / x
	yz.s = S - x
	### robust:
	# x*y^2 - (x^2+yz)*y + (x^2+yz)*yz.s + b1*S - R1
	# x*y^2 + x*yz - x*y*yz.s = 0 # x*y*(y+z - yz.sum) = 0
	# (x^2+yz - x*yz.s)*y + x*yz - (x^2+yz)*yz.s - b1*S + R1 = 0
	y = - (x*yz - (x^2+yz)*yz.s - b[1]*S + R[1]) / (x^2+yz - x*yz.s)
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}
test.ht3 = function(x, y, z, R, b=0) {
	if(missing(y)) {
		y = x[,2]; z = x[,3]; x = x[,1];
	}
	### Test
	err1 = x*y^2 + y*z^2 + z*x^2 + b[1]*(x+y+z) # - R[1] # = 0
	err2 = x*y + y*z + z*x + if(length(b) < 2) 0 else b[2]*(x+y+z) # - R[2] # = 0
	err3 = x*y*z + if(length(b) < 3) 0 else b[3]*(x+y+z) # - R[3] # = 0
	err = rbind(err1, err2, err3)
	if( ! missing(R)) {
		err = err - rep(R, each=length(x))
	}
	err = round0(err)
	return(err)
}

### Examples

###
R = c(1, 1, 1)
b = 0
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -1 + 3*x - 3*x^2 + 4*x^3 + x^4 - 4*x^5 + 11*x^6 - 4*x^7 + x^9
round0(err)

###
R = c(0, 1, 1)
sol = solve.ht3(R)

### Test
test.ht3(sol)

poly.calc(sol[,1])
x = sol[,1]
err = -1 + 3*x - 3*x^2 + 4*x^3 - 3*x^5 + 7*x^6 - 3*x^7 + x^9
round0(err)


###############
### Extensions:

### Extension A1:
### x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1
### Extension A2:
### x*y + x*z + y*z + b2*(x+y+z) = R2;
### Extension A3:
### x*y*z + b3*(x+y+z) = R3;

### * x*z^2 =>
x^2*y^2*z^2 + x*y*z^4 + z^3*x^3 + b1*(x+y+z)*x*z^2 - R1*x*z^2 # = 0
R3^2 + R3*z^3 + z^3*x^3 + b1*S*x*z^2 - R1*x*z^2 # = 0
### Sum() =>
3*R3^2 + R3*(x^3+y^3+z^3) + (x^3*z^3+x^3*y^3+y^3*z^3) +
	+ (b1*S - R1)*(x*z^2 + x^2*y + y^2*z) # = 0
### Sum - R1*initial Eq
R3*(x^3+y^3+z^3) + (x^3*z^3+x^3*y^3+y^3*z^3) +
	(b1*S - R1)*(x*z^2 + x^2*y + y^2*z + x*y^2 + y*z^2 + z*x^2) + # = (E2*S - 3*E3)
	+ b1*S*(b1*S - R1) - R1*(b1*S - R1) + 3*R3^2 # = 0
R3*(S^3 - 3*R2*S + 3*R3) + (R2^3 - 3*R3*(R2*S - R3)) +
	(b1*S - R1)*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) +
	+ b1^2*S^2 - 2*R1*b1*S + R1^2 + 3*R3^2 # = 0
R3*S^3 + (b1*S - R1)*((x^2+y^2+z^2)*(x+y+z) - (x^3+y^3+z^3)) +
	+ b1^2*S^2 - 6*R2*R3*S - 2*R1*b1*S + R1^2 + R2^3 + 9*R3^2 # = 0
R3*S^3 + (b1*S - R1)*(R2*S - 3*R3) +
	+ b1^2*S^2 - 6*R2*R3*S - 2*R1*b1*S + R1^2 + R2^3 + 9*R3^2 # = 0
R3*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0

### Extension A3:
# - includes A1, but A2 was missed;
#   [see function solve.ht3() for complete variant]
# E3 = R3 - b3*S
E3*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*E3 + 6*R2*E3 + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*E3^2 + 3*R1*E3 # = 0
(R3 - b3*S)*S^3 + (b1^2 + b1*R2)*S^2 - (R1*R2 + 3*b1*(R3 - b3*S) + 6*R2*(R3 - b3*S) + 2*b1*R1)*S +
	+ R1^2 + R2^3 + 9*(R3 - b3*S)^2 + 3*R1*(R3 - b3*S) # = 0
-b3*S^4 + R3*S^3 + (b1^2 + 3*b1*b3 + 9*b3^2 + b1*R2 + 6*b3*R2)*S^2 +
	- (R1*R2 + 3*b1*R3 + 6*R2*R3 + 2*b1*R1 + 3*b3*R1 + 18*b3*R3)*S +
	+ R1^2 + R2^3 + 9*R3^2 + 3*R1*R3 # = 0


### Examples

###
R = c(1, 1, 1); b = 1;
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

err = -1 + 3*x - x^2 + 8*x^4 - 13*x^5 + 15*x^6 - 9*x^7 + 2*x^8 + x^9
round0(err)


###
R = c(0, 1, -1); b = 1;
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

err = 1 + 3*x + x^2 - 5*x^4 - 10*x^5 - 11*x^6 - 6*x^7 - 2*x^8 + x^9
round0(err)


###
# R3 = (b[2]^2 + b[1]*b[2]) +/- 1;
R = c(0, 1, 9); b = c(1, 2);
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

# numeric instability with 1st root!
err = -729 - 19521*x - 8649*x^2 - 6428*x^3 - 4979*x^4 - 2234*x^5 - 653*x^6 - 154*x^7 - 122*x^8 + x^9
round0(err)


###
# R3 = (b[2]^2 + b[1]*b[2]) +/- 1;
R = c(0, 1, -2); b = c(2, -1);
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

poly.calc(x)

err = 8 - 72*x + 60*x^2 - 42*x^3 + 60*x^4 - 12*x^5 + 31*x^6 + 9*x^7 + 21*x^8 + x^9
round0(err)


### Ext 3:
R = c(1, 1, 0); b = c(-7, 1, 1);
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -2 + 15*x - 65*x^2 + 71*x^3 - 15*x^4 + 96*x^5 + 14*x^6 - 95*x^7 - 57*x^8 - 87*x^9 - 36*x^10 + x^12
round0(err)


### Ext 3, ex 2:
R = c(0, 1, 0); b = c(-7, 1, 1);
sol = solve.ht3(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3(x, y, z, b=b)

round0.p(poly.calc(x))

err = -1 + x - 32*x^2 + 67*x^3 - 39*x^4 + 86*x^5 - 4*x^6 - 79*x^7 - 25*x^8 - 75*x^9 - 35*x^10 + x^12
round0(err)


#########################
#########################

### x*y^n + y*z^n + z*x^n = R1

############
### Order 3: n = 3
x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Eq 1: * x*z^3
x^2*y^3*z^3 + x*y*z^6 + z^4*x^4 - R1*x*z^3 # = 0
R3^2*y*z + R3*z^5 + x^4*z^4 - R1*x*z^3 # = 0
# similar:
R3^2*x*z + R3*x^5 + x^4*y^4 - R1*y*x^3 # = 0
R3^2*x*y + R3*y^5 + y^4*z^4 - R1*z*y^3 # = 0

### Sum =>
R3^2*(x*y+x*z+y*z) + R3*(x^5+y^5+z^5) +
	+ (x^4*y^4+x^4*z^4+y^4*z^4) - R1*(y*x^3+z*y^3+x*z^3) # = 0
R3^2*R2 + R3*(S^5 - 5*R2*S^3 + 5*R3*S^2 + 5*R2^2*S  - 5*R2*R3) +
	+ (x^4*y^4+x^4*z^4+y^4*z^4) - R1*(y*x^3+z*y^3+x*z^3) # = 0
R3*S^5 - 5*R2*R3*S^3 + 5*R3^2*S^2 + 5*R2^2*R3*S +
	+ (4*R2*R3^2 + R2^4 - 4*R2^2*R3*S + 2*R3^2*S^2) - R1*(y*x^3+z*y^3+x*z^3) - 4*R2*R3^2 # = 0
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(y*x^3+z*y^3+x*z^3) + R2^4 # = 0
### Sum - R1*Initial_Eq =>
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(y*x^3+z*x^3+x*y^3+z*y^3+x*z^3+y*z^3) + R1^2 + R2^4 # = 0
R3*S^5 - 5*R2*R3*S^3 + 7*R3^2*S^2 + R2^2*R3*S +
	- R1*(R2*(S^2 - 2*R2) - R3*S) + R1^2 + R2^4 # = 0
R3*S^5 - 5*R2*R3*S^3 + (7*R3^2 - R1*R2)*S^2 + (R2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*R2^2 + R2^4 # = 0


### Solution
solve.ht3 = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3], 0, - 5*R[2]*R[3], (7*R[3]^2 - R[1]*R[2]),
			(R[2]^2*R[3] + R[1]*R[3]), R[1]^2 + 2*R[1]*R[2]^2 + R[2]^4)
	} else {
		# TODO
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=5, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	if(R[1] == 0) {
		# with chain rule!
		# x3 = (5*R[3]*S^4 - 15*R[2]*R[3]*S^2 + 14*R[3]^2*S + R[2]^2*R[3]) # * dS/dR
		dS = - R[2]*S^2 + R[3]*S + 2*R[2]^2
		x3 = - dS
	} else {
		x3 = (R[3]*S^5 - 5*R[2]*R[3]*S^3 + 7*R[3]^2*S^2 + R[2]^2*R[3]*S + R[2]^4) / R[1]
	}
	yz.d = (x3 - R[1]) / (x^3 + yz*yz.s - x*(yz.s^2 - yz))
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 1, -2);
sol = solve.ht3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 32 + 80*x + 80*x^2 + 120*x^3 + 130*x^4 + 61*x^5 + 36*x^6 + 6*x^7 - 9.5*x^8 - 17*x^9 +
	- 19*x^10 - 3*x^11 - 3.5*x^12 + x^15
round0(err)


### Ex 2:
R = c(0, 2, -2);
sol = solve.ht3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x*y^3 + y*z^3 + z*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0

poly.calc(x)

err = 32 + 160*x + 320*x^2 + 400*x^3 + 400*x^4 + 272*x^5 + 104*x^6 + 8*x^7 - 48*x^8 - 48*x^9 +
	- 44*x^10 - 16*x^11 - 4*x^12 + x^15
round0(err)


#############################
#############################

### Generalization:
### x^p*y^n + y^p*z^n + z^p*x^n

#####################
### Higher Powers ###
### p > 1         ###

### p = 2
### x^2*y^n + y^2*z^n + z^2*x^n = R1

### Order 3: n = 3, p = 2
x^2*y^3 + y^2*z^3 + z^2*x^3 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

### Eq 1: * x^3*y^2
x^5*y^5 + x^3*y^5*z^3 + x^6*y^2*z^2 - R1*x^3*y^2 # = 0
x^5*y^5 + R3^3*y + R3^2*x^4 - R1*x^3*y^2 # = 0
# similar:
x^5*z^5 + R3^3*x + R3^2*z^4 - R1*x^2*z^3 # = 0
y^5*z^5 + R3^3*z + R3^2*y^4 - R1*y^3*z^2 # = 0

### Sum =>
R3^3*(x+y+z) + R3^2*(x^4+y^4+z^4) +
	+ (x^5*y^5+x^5*z^5+y^5*z^5) - R1*(y^2*x^3+z^2*y^3+x^2*z^3) # = 0
R3^3*S + R3^2*(S^4 - 4*R2*S^2 + 4*R3*S + 2*R2^2) +
	+ (5*R2^2*R3^2 + R2^5 - 5*R2^3*R3*S - 5*R3^3*S + 5*R2*R3^2*S^2) +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3) # = 0
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3) + R2^5 + 7*R2^2*R3^2 # = 0

# Sum() - R1*(initial Eq) + R1^2 =>
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*(y^2*x^3+z^2*y^3+x^2*z^3 + x^2*y^3 + y^2*z^3 + z^2*x^3) +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 # = 0
R3^2*S^4 + R2*R3^2*S^2 - 5*R2^3*R3*S +
	- R1*((R2^2 - 2*R3*S)*S - R3*R2) +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 # = 0
R3^2*S^4 + (2*R1*R3 + R2*R3^2)*S^2 - (R1*R2^2 + 5*R2^3*R3)*S +
	+ R1^2 + R2^5 + 7*R2^2*R3^2 + R1*R2*R3 # = 0

### Solution
solve.ht3 = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3]^2, 0, (2*R[1]*R[3] + R[2]*R[3]^2), - (R[1]*R[2]^2 + 5*R[2]^3*R[3]),
			R[1]^2 + R[2]^5 + 7*R[2]^2*R[3]^2 + R[1]*R[2]*R[3])
	} else {
		# TODO
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=4, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	if(R[1] == 0) {
		# with chain rule!
		dS = 2*R[3]*S^2 - R[2]^2*S + R[2]*R[3]
		x3 = - dS
	} else {
		x3 = (R[3]^2*S^4 + R[2]*R[3]^2*S^2 - 5*R[2]^3*R[3]*S + R[2]^5 + 7*R[2]^2*R[3]^2) / R[1]
	}
	yz.d = (x3 - R[1]) / (x^3*yz.s + yz^2 - x^2*(yz.s^2 - yz))
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 1, -1);
sol = solve.ht3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^3 + y^2*z^3 + z^2*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 1 + 4*x + 6*x^2 + 8*x^3 + 12*x^4 + 10*x^5 + 13*x^6 + 14*x^7 + 12*x^8 + 8*x^9 + 3*x^10 + x^12
round0(err)


### Ex 2:

R = c(0, 2, -1);
sol = solve.ht3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2*y^3 + y^2*z^3 + z^2*x^3 # - R[1] # = 0
x*y + y*z + z*x # - R[2] # = 0
x*y*z # - R[3] # = 0


poly.calc(x)

err = 1 + 8*x + 24*x^2 + 36*x^3 + 42*x^4 + 56*x^5 + 86*x^6 + 108*x^7 + 92*x^8 + 44*x^9 + 10*x^10 + x^12
round0(err)


###############
### Extensions:

### TODO: A1, A2, A3;


#############################
#############################

###############
### Dual Eq ###
### p = 1   ###

### x*y^n + y*z^n + z*x^n = R1
### x*z^n + y*x^n + z*y^n = R2

############
### Order 2: n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*z^2 + y*x^2 + z*y^2 - R2 # = 0
x*y*z - R3 # = 0

### Extensions:
### Extension M3:
# x*y*z*(x+y+z) = R3
### Extension D3:
# x*y*z / (x+y+z) = R3

### Solution:
### Eq 1 + Eq 2 =.
x*y^2 + y*z^2 + z*x^2 + x*z^2 + y*x^2 + z*y^2 - R1 - R2 # = 0
E2*S - 3*E3 - R1 - R2 # = 0
E2 = (R1 + R2 + 3*E3) / S

### see Section "Simple System":
R3*S^3 - (R1+6*R3)*E2*S + R1^2 + E2^3 + 9*R3^2 + 3*R1*R3 # = 0
R3*S^3 - (R1+6*R3)*(R1 + R2 + 3*R3) + R1^2 + (3*R3 + R1 + R2)^3 / S^3 + 9*R3^2 + 3*R1*R3 # = 0
R3*S^6 - ((R1+6*R3)*(R1 + R2 + 3*R3) - R1^2 - 9*R3^2 - 3*R1*R3)*S^3 + (R1 + R2 + 3*R3)^3 # = 0
R3*S^6 - (R1*R2 + 6*R2*R3 + 9*R3^2 + 6*R1*R3)*S^3 + (R1 + R2 + 3*R3)^3 # = 0

### Extension M3:
R3*S^5 - (R1*R2 + 6*R2*R3/S + 9*R3^2/S^2 + 6*R1*R3/S)*S^3 + (R1 + R2 + 3*R3/S)^3 # = 0
R3*S^8 - (R1*R2 + 6*R2*R3/S + 9*R3^2/S^2 + 6*R1*R3/S)*S^6 + ((R1 + R2)*S + 3*R3)^3 # = 0
R3*S^8 - R1*R2*S^6 - 6*(R1*R3+R2*R3)*S^5 - 9*R3^2*S^4 + (R1 + R2)^3*S^3 +
	+ 9*R3*(R1 + R2)^2*S^2 + 27*R3^2*(R1 + R2)*S + 27*R3^3 # = 0

### Solution
solve.ht3Dual = function(R, b=0, type) {
	# TODO: R1 == R2
	if(missing(type)) {
		type = 1;
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(R[3], 0, 0, - (R[1]*R[2] + 6*R[2]*R[3] + 9*R[3]^2 + 6*R[1]*R[3]), 0, 0, (R[1] + R[2] + 3*R[3])^3)
	} else {
		# TODO
	}
	} else if(match("M3", type) > 0) {
		type = 3;
		R12 = R[1] + R[2]
		coeff = c(R[3], 0, - R[1]*R[2], - 6*R[3]*(R12), - 9*R[3]^2, (R12)^3,
			9*R[3]*(R12)^2, 27*R[3]^2*(R12), 27*R[3]^3)
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	len = length(S)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	E3 = if(type == 1) R[3] else if(type == 3) R[3] / S; # TODO: Ext 2;
	E2 = (R[1] + R[2] + 3*E3) # TODO: Ext 2;
	if(any(S == 0)) {
		print("Warning: Div by 0!")
		E2.zero = - rootn(R[1]^2 + 9*R[3]^2 + 3*R[1]*R[3], 3);
		E2 = ifelse(S == 0, E2.zero, E2 / S);
		# E3 cannot be 0!
		E3 = R[3]
		x = sapply(seq(length(S)), function(id) roots(c(1, -S[id], E2[id], -E3)))
		len = length(S)
	} else {
		E2 = E2 / S;
		E3 = if(type == 1) rep(E3, length(S)) else E3;
		x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	}
	S  = matrix(S,  ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	yz = E3/x
	yz.s = S - x
	### robust:
	yz.d = (R[2] - R[1]) / (x^2 + yz - x*yz.s)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 3, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z)

round0.p(poly.calc(x))

err = 1 - 37*x^6 + 76*x^9 - 37*x^12 + x^18 # trivial;
round0(err)

### Classic Polynomial:
R1 = R[1]; R2 = R[2]; R3 = R[3];
R3*x^18 + (3*R3^2 - R1*R2)*x^15 + (R1^3 + R2^3 - 5*R1*R2*R3 + 6*R3^3)*x^12 +
	(-R1^2*R2^2 + 2*R1^3*R3 + 2*R2^3*R3 - 6*R1*R2*R3^2 + 7*R3^4)*x^9 +
	(R1^3*R3^2 + R2^3*R3^2 - 5*R1*R2*R3^3 + 6*R3^5)*x^6 +
	(-R1*R2*R3^4 + 3*R3^6)*x^3 + R3^7


### Ex 2:
R = c(1, 2, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=2)

# there are some issues!
poly.calc(x[1:9]^3)
poly.calc(x[10:18])

err = 1 + 3*x^3 - 4*x^6 + x^9 # trivial;
round0(err)

#################
### Extension M3:

### Ex 1:
R = c(1, 3, -1);
sol = solve.ht3Dual(R, type="M3")
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=2, type="M3")

round0.p(poly.calc(x)) * 27

err = 1 + 20*x + 166*x^2 + 736*x^3 + 1862*x^4 + 2596*x^5 + 1305*x^6 - 2488*x^7 - 7403*x^8 +
	- 8428*x^9 - 2887*x^10 + 3848*x^11 + 9086*x^12 + 7032*x^13 - 81*x^14 - 1160*x^15 - 3402*x^16 +
	- 108*x^17 - 180*x^18 - 432*x^19 - 27*x^20 + 81*x^22 + 27*x^24;
round0(err)


### Ex 2:
R = c(0, 3, -3);
sol = solve.ht3Dual(R, type="M3")
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=2, type="M3")

round0.p(poly.calc(x))

err = 27 + 135*x + 270*x^2 + 270*x^3 + 108*x^4 - 135*x^5 - 405*x^6 - 522*x^7 - 315*x^8 - 18*x^9 +
	+ 135*x^10 + 231*x^11 + 222*x^12 + 72*x^13 - 8*x^15 - 27*x^16 - 9*x^17 + 3*x^18 - 9*x^19 +
	- 3*x^20 + x^24
round0(err)


### TODO: R1 == R2;


### Classic Polynomial:
# TODO


############
### Order 2: n = 2
### Combined Variant
(x*y^2 + y*z^2 + z*x^2) - b*(x*z^2 + y*x^2 + z*y^2) - R1 # = 0
x*y + x*z + y*z - R2 # = 0
x*y*z - R3 # = 0

### Solution:
(b+1)*(x*y^2 + y*z^2 + z*x^2) - b*(E2*S - 3*E3) - R1 # = 0
(b+1)*(x*z^2 + y*x^2 + z*y^2) - (E2*S - 3*E3) + R1 # = 0
# =>
R3*(x^3+y^3+z^3) + (x^3*y^3+x^3*z^3+y^3*z^3) +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*(S^3 - 3*E2*S + 3*E3) + (R2^3 - 3*R3*(R2*S - R3)) +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) + 3*R3^2 # = 0
R3*S^3 - 6*R2*R3*S + R2^3 + 9*R3^2 +
	- 1/(b+1) * (b*(E2*S - 3*E3) + R1)*(x^2*y + y^2*z + x*z^2) # = 0
R3*S^3 - 6*R2*R3*S + R2^3 + 9*R3^2 +
	- 1/(b+1)^2 * (b*(E2*S - 3*E3) + R1)*(E2*S - 3*E3 - R1) # = 0
(b+1)^2*R3*S^3 - (b+1)^2*6*R2*R3*S + (b+1)^2*R2^3 + 9*(b+1)^2*R3^2 +
	- (b*(E2*S - 3*E3) + R1)*(E2*S - 3*E3 - R1) # = 0
(b+1)^2*R3*S^3 - (b+1)^2*6*R2*R3*S + (b+1)^2*R2^3 + 9*(b+1)^2*R3^2 +
	- b*(R2*S - 3*R3)^2 + (b-1)*R1*(R2*S - 3*R3) + R1^2 # = 0
(b+1)^2*R3*S^3 - b*R2^2*S^2 - 6*(b+1)^2*R2*R3*S + 6*b*R2*R3*S + (b-1)*R1*R2*S +
	+ R1^2 + (b+1)^2*R2^3 + 9*((b+1)^2-b)*R3^2 - 3*(b-1)*R1*R3 # = 0


### Solution
solve.ht3Combi = function(R, a, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c((a+1)^2*R[3], - a*R[2]^2, - 6*(a+1)^2*R[2]*R[3] + 6*a*R[2]*R[3] + (a-1)*R[1]*R[2],
			R[1]^2 + (a+1)^2*R[2]^3 + 9*((a+1)^2-a)*R[3]^2 - 3*(a-1)*R[1]*R[3])
	} else {
		# TODO
		print("Not yet implemented!")
	}
	S = roots(coeff)
	len = length(S)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x
	yz.s = S - x
	### robust:
	R1A = (a*(R[2]*S - 3*R[3]) + R[1]) / (a+1) # TODO: a = -1
	R1B = ((R[2]*S - 3*R[3]) - R[1]) / (a+1) # TODO: a = -1
	yz.d = - (R1A - R1B) / (x^2 + yz - x*yz.s)
	y = (yz.s + yz.d)/2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Example:
a = 1
R = c(1,1,1)
sol = solve.ht3Combi(R, a=a)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
(x*y^2 + y*z^2 + z*x^2) - a*(x*z^2 + y*x^2 + z*y^2) # - R[1] # = 0
x*y + x*z + y*z # - R[2] # = 0
x*y*z # - R[3] # = 0

round0.p(poly.calc(x))

err = -1 + 3*x - 3.25*x^2 + 4.5*x^3 - 1.75*x^4 - x^5 + 4.5*x^6 - 1.5*x^7 - 0.25*x^8 + x^9
round0(err)


##################
##################

############
### Order 3: n = 3
x*y^3 + y*z^3 + z*x^3 - R1 # = 0
x*z^3 + y*x^3 + z*y^3 - R2 # = 0
x*y*z - R3 # = 0

### Solution:
### Eq 1 + Eq 2 =>
x*y^3 + y*z^3 + z*x^3 + x*z^3 + y*x^3 + z*y^3 - R1 - R2 # = 0
E2*(S^2 - 2*E2) - R3*S - R1 - R2 # = 0
2*E2^2 - E2*S^2 + R3*S + R1 + R2 # = 0
# E2 = (S^2 +/- sqrt(Det))/4

### see Section "Simple System":
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
- 128*R1*R2 + 64*R1^2 + 64*R2^2 - 32*Dsq*R1*S^2 - 32*Dsq*R2*S^2 + 1728*R3^2*S^2 - 320*Dsq*R3*S^3 +
	- 64*R1*S^4 - 64*R2*S^4 - 96*R3*S^5 + 8*Dsq*S^6 + 8*S^8 # = 0
4*Dsq*R1*S^2 + 4*Dsq*R2*S^2 + 40*Dsq*R3*S^3 - Dsq*S^6 =
	S^8 - 12*R3*S^5 - 8*R1*S^4 - 8*R2*S^4 + 216*R3^2*S^2 - 16*R1*R2 + 8*R1^2 + 8*R2^2
R3*S^13 - 26*R3^2*S^10 - (13*R1*R3 + 13*R2*R3)*S^9 - R1*R2*S^8 + 119*R3^3*S^7 +
	+ (186*R1*R3^2 + 186*R2*R3^2)*S^6 + (90*R1*R2*R3 + 39*R1^2*R3 + 39*R2^2*R3)*S^5 +
	+ (8*R1*R2^2 + 8*R1^2*R2 + 729*R3^4)*S^4 +
	(- 108*R1*R2*R3^2 + 54*R1^2*R3^2 + 54*R2^2*R3^2)*S^2 +
	(- 4*R1*R2^3 + 6*R1^2*R2^2 - 4*R1^3*R2 + R1^4 + R2^4)
### 2*E2^2 = E2*S^2 - R3*S - R1 - R2
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (E2^2*R3 + R1*R3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
R3*S^5 - 5*E2*R3*S^3 + (7*R3^2 - R1*E2)*S^2 + (1/2*(E2*S^2 - R3*S - R1 - R2)*R3 + R1*R3)*S +
	+ R1^2 + R1*(E2*S^2 - R3*S - R1 - R2) + 1/4 * (E2*S^2 - R3*S - R1 - R2)^2 # = 0
E2*(S^6 - 40*R3*S^3 - 4*S^2*(R1 + R2)) +
	7*R3*S^5 - (R1 + R2)*S^4 + 54*R3^2*S^2 + 2*(R1 - R2)^2 # = 0


### Solution
solve.ht3Dual = function(R, b=0) {
	if(length(b) == 1 && b[1] == 0) {
		coeff = c(
			R[3], 0, 0, - 26*R[3]^2, - (13*R[1]*R[3] + 13*R[2]*R[3]), - R[1]*R[2], 119*R[3]^3,
			(186*R[1]*R[3]^2 + 186*R[2]*R[3]^2), (90*R[1]*R[2]*R[3] + 39*R[1]^2*R[3] + 39*R[2]^2*R[3]),
			(8*R[1]*R[2]^2 + 8*R[1]^2*R[2] + 729*R[3]^4), 0,
			(- 108*R[1]*R[2]*R[3]^2 + 54*R[1]^2*R[3]^2 + 54*R[2]^2*R[3]^2), 0,
			(- 4*R[1]*R[2]^3 + 6*R[1]^2*R[2]^2 - 4*R[1]^3*R[2] + R[1]^4 + R[2]^4))
	} else {
		# TODO
	}
	if(length(b) > 1) {
		# Ext 2:
		# TODO
	}
	S = roots(coeff)
	len = length(S)
	print(S)
	b2 = if(length(b) > 1) b[2] else 0; # TODO: Ext 2;
	# TODO: Ext 2;
	# E2 = (S^2 + sign*sqrt(S^4 - 8*R[3]*S - 8*R[1] - 8*R[2]))/4
	# robust
	div = - (S^6 - 40*R[3]*S^3 - 4*S^2*(R[1] + R[2]));
	E2 = (7*R[3]*S^5 - (R[1] + R[2])*S^4 + 54*R[3]^2*S^2 + 2*(R[1] - R[2])^2) / div
	#
	x = sapply(seq(length(S)), function(id) roots(c(1, -S[id], E2[id], -R[3])))
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	### TODO: Div by 0
	yz = R[3]/x
	yz.s = S - x
	### robust:
	yz.d = (R[2] - R[1]) / (x^3 + yz*yz.s - x*(yz.s^2 - yz))
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Examples:

R = c(1, 3, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=3)

# order 39
round0.p(poly.calc(x))

# order 39
err = 1 + 5*x^3 - 6*x^4 + 10*x^6 - 9*x^7 + 9*x^8 + 11*x^9 - 41*x^10 + 4*x^11 - 52*x^12 + 168*x^13 +
	+ 100*x^14 + 20*x^15 - 65*x^16 - 41*x^17 + 28*x^18 + 250*x^19 + 85*x^20 - 27*x^21 - 25*x^22 +
	+ 3*x^23 + 33*x^24 + 65*x^25 - 54*x^26 - 30*x^27 + 3*x^28 - x^30 - 16*x^32 - 4*x^33 + 3*x^34 + x^39
round0(err)


### Ex 2:
R = c(0, -1, -1);
sol = solve.ht3Dual(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.ht3Dual(x, y, z, n=3)

# order 39
round0.p(poly.calc(x))


