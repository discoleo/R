########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric S3:
### Mixt Type
###
### draft v.0.1b


### Heterogenous Symmetric
### Polynomial Systems: 3 Variables
### Mixt: Hetero + Symmetric

#####################

### History

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

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

######################

### Mixt Systems Type 1

### x*y^n + y*z^n + z*x^n = R1

### Order 2: n = 2
x*y^2 + y*z^2 + z*x^2 - R1 # = 0
x*y + y*z + z*x - R2 # = 0
x*y*z - R3 # = 0

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
	} else {
		coeff = c(R[3], (b[1]^2 + b[1]*R[2]), - (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	}
	if(length(b) > 1) {
		# Ext 2:
		coeff = coeff - c(b[2]^3 + b[1]*b[2], -b[2]*(R[1] + 3*b[2]*R[2] + 6*R[3]), 3*b[2]*R[2]^2, 0)
	}
	S = roots(coeff)
	b2 = if(length(b) > 1) b[2] else 0; # Ext 2;
	x = sapply(S, function(x) roots(c(1,-x, R[2] - b2*x, -R[3])))
	S = matrix(S, ncol=3, nrow=3, byrow=T)
	yz = R[3]/x
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
	err3 = x*y*z # - R[3] # = 0
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
sol = solve.ht3(R)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
test.ht3(x, y, z)

poly.calc(x)

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

### Extension 1:
### x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1
### Extension 2:
### x*y + x*z + y*z + b2*(x+y+z) = R2;

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


#########################
#########################

### x*y^n + y*z^n + z*x^n = R1

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
		# TODO
	} else {
		x3 = (R[3]*S^5 - 5*R[2]*R[3]*S^3 + 7*R[3]^2*S^2 + R[2]^2*R[3]*S + R[2]^4) / R[1]
		yz.d = (x3 - R[1]) / (x^3 + yz*yz.s - x*(yz.s^2 - yz))
	}
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


