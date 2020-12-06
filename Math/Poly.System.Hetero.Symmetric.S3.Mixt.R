########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogenous Symmetric S3:
### Mixt Type
###
### draft v.0.1a-ext1


### Heterogenous Symmetric
### Polynomial Systems: 3 Variables
### Mixt: Hetero + Symmetric

#####################

### History

### draft v.0.1a-ext1:
# - Extension:
#   x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Hetero.Symmetric.S3.Mixt.R;
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

### n = 2
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
	if(b == 0) {
		coeff = c(R[3], 0, - (R[1]+6*R[3])*R[2], R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	} else {
		coeff = c(R[3], (b[1]^2 + b[1]*R[2]), - (R[1]*R[2] + 3*b[1]*R[3] + 6*R[2]*R[3] + 2*b[1]*R[1]),
			R[1]^2 + R[2]^3 + 9*R[3]^2 + 3*R[1]*R[3])
	}
	S = roots(coeff)
	x = sapply(S, function(x) roots(c(1,-x,R[2],-R[3])))
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
	err2 = x*y + y*z + z*x # - R[2] # = 0
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

### x*y^2 + y*z^2 + z*x^2 + b1*(x+y+z) = R1

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


