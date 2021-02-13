########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogenous Symmetric
### with Composite Leading Term
###
### draft v.0.1a


### Hetero-Symmetric
### Polynomial Systems: 3 Variables
### Composite Leading Term

### Example:
x^n*y^m + P(x, y, z) = R
y^n*z^m + P(y, z, x) = R
z^n*x^m + P(z, x, y) = R


###############
### History ###
###############

### draft v.0.1a:
# - moved to this file from:
#   Poly.System.Hetero.Symmetric.S3.R;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R


################################
################################

######################
### X*Y High Power ###
######################

### x*y + b*y = R

# x*y + b1*y = R
# y*z + b1*z = R
# z*x + b1*x = R

### Extensions:

### Simple Extension: type A1
# x*y + b1*y + be2*(x+y+z)^2 + be1*(x+y+z) = R
### Structural Extension:
# a*x*y*z + x*y + b1*y = R


### Solution:

### Note:
# - simple system: has NO solutions x != y != z;
# - type A extensions: have such solutions;

### Sum =>
E2 + b1*S - 3*R # = 0

### Sum(z*...) =>
3*E3 + b1*E2 - R*S # = 0
# b1*E2 = R*S - 3*E3

### Diff =>
# y*(x-z) = -b1*(y-z)
# z*(x-y) = -b1*(x-z)
# x*(y-z) =  b1*(x-y)
### Prod =>
# E3 = b1^3

### Eq:
b1*E2 + b1^2*S - 3*b1*R # = 0
R*S - 3*E3 + b1^2*S - 3*b1*R # = 0
(R + b1^2)*S - 3*b1*(b1^2 + R) # = 0
S = 3*b1; # is a FALSE solution;


### Solver:

solve.CHP.S3P1 = function(R, b, b.ext=0, a=0, debug=TRUE) {
	be1 = b.ext[1];
	be2 = if(length(b.ext) > 1) b.ext[2] else 0;
	a1 = a[1];
	a2 = if(length(a) > 1) a[2] else 0;
	if(be1 == 0 && be2 == 0) {
		stop("NO solutions: x != y != z")
		S = 3*b[1];
	} else {
		# S = roots(c(-be2, 3*b[1]*be2 - be1, R + b[1]^2 + 3*b[1]*be1, -3*b[1]^3 - 3*b[1]*R))
		# isWrong = round0(S - 3*b[1]) == 0
		# S = S[ ! isWrong]
		S = roots(c(be2, be1, -R[1] - b[1]^2 + a1*b[1]^3 + a2*b[1]^6))
	}
	if(debug) print(S);
	E3 = b[1]^3 - 0*S;
	R1 = R[1] - be1*S - be2*S^2 - a1*E3 - a2*E3^2;
	E2 = 3*R1 - b[1]*S;
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	### robust
	if(len == 1) {
		S = rep(S, 3); R1 = rep(R1, 3);
	} else {
		S  = matrix(S, ncol=len, nrow=3, byrow=T)
		R1 = matrix(R1, ncol=len, nrow=3, byrow=T)
		E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	}
	#
	yz = E3 / x;
	z = (R1 - yz) / b[1];
	y = S - x - z;
	# Note: 1 set of solutions is incorrect!
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol)
}
test.CHP.S3P1 = function(sol, R, b, b.ext=0, a=0) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0);
	if(length(a) < 2) a = c(a, 0)
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	xyz = x*y*z; a.ext = a[1]*xyz + a[2]*xyz^2;
	s = (x+y+z); s.ext = b.ext[1]*s + b.ext[2]*s^2;
	ext = a.ext + s.ext;
	err1 = x*y + b[1]*y + ext # - R
	err2 = y*z + b[1]*z + ext # - R
	err3 = z*x + b[1]*x + ext # - R
	round0(rbind(err1, err2, err3))
}

### Examples:
R = -1
b = 3
b.ext = c(1, 0)
sol = solve.CHP.S3P1(R, b, b.ext=b.ext)

### Test
test.CHP.S3P1(sol, R, b, b.ext)


### Ex 2:
R = -1
b = 3
b.ext = c(1, 3)
sol = solve.CHP.S3P1(R, b, b.ext=b.ext)

### Test
test.CHP.S3P1(sol, R, b, b.ext)
round0.p(poly.calc(sol[1:6, 1]))


############
### Str Ext: Ex 1
R = -1
b = 3
b.ext = c(1, 3); a = 2
sol = solve.CHP.S3P1(R, b, b.ext=b.ext, a=a)

### Test
test.CHP.S3P1(sol, R, b, b.ext, a=a)
round0.p(poly.calc(sol[1:6, 1]))


############
### Str Ext: Ex 2
R = -1
b = 3
b.ext = c(1, 3); a = c(-1, 1)
sol = solve.CHP.S3P1(R, b, b.ext=b.ext, a=a)

### Test
test.CHP.S3P1(sol, R, b, b.ext, a=a)
round0.p(poly.calc(sol[1:6, 1]))


############
### Str Ext: Ex 3
R = -1
b = 3
b.ext = c(0, -1); a = c(2, -1)
sol = solve.CHP.S3P1(R, b, b.ext=b.ext, a=a)

### Test
test.CHP.S3P1(sol, R, b, b.ext, a=a)
round0.p(poly.calc(sol[1:6, 1]))


#######################
#######################


#######################
### Mixt-Order: 2+1 ###

### x[i]^2*x[j] + b*x[j]

### TODO:


#######################

### Trivial
### x[i]^2*x[j] + b*Sum

# x^2*y + b1*(x+y+z) = R
# y^2*z + b1*(x+y+z) = R
# z^2*x + b1*(x+y+z) = R

### Solution:

# - trivial P9 polynomial;

### Diff =>
# x^2 = y*z
# y^2 = x*z
# z^2 = x*y

# y = x^2/z
# => x^4/z^2 = x*z
# => x^3 = z^3
# => y^3 = z^3

# Case x != y != z
# y = x*m
# z = x*m^2
# => 1 + m + m^2 = 0!
# x^3*m = R

m3 = unity(3, all=FALSE)

### Example 1:

b = 3
R = 1
#
x = (R/m3)^(1/3) * c(1, m3, m3^2)
y = x * m3
z = x * m3^2
sol = cbind(x,y,z)
sol

### Test
x^2*y + b[1]*(x+y+z)
y^2*z + b[1]*(x+y+z)
z^2*x + b[1]*(x+y+z)

### Classical Polynomial
x^9 - R^3


##########################
##########################

### Shifted
### x[i]^2 * (x[j] - shift) + b*Sum

# x^2*(y - s) + b1*(x+y+z) = R
# y^2*(z - s) + b1*(x+y+z) = R
# z^2*(x - s) + b1*(x+y+z) = R

### Solution

# Diff =>
# x^2*(y - s) = y^2*(z - s)
# y^2*(z - s) = z^2*(x - s)
# z^2*(x - s) = x^2*(y - s)
# =>
# x^4*(y-s)^2 = y^2*z^2*(x-s)*(z-s)

### TODO


#######################
#######################

#######################
### X*Y*Z-type Term ###
#######################

### x^2*y*z + b*y

# x^2*y*z + b1*y = R
# x*y^2*z + b1*z = R
# x*y*z^2 + b1*x = R

### Solution

### Diff =>
# x*y*z*(x-y) = - b1*(y-z)
# x*y*z*(y-z) = - b1*(z-x)
# x*y*z*(z-x) = - b1*(x-y)
### Prod =>
# (x*y*z)^3 = -b1^3
# E3 = -b1 * m3;

### Sum =>
# x*y*z*(x+y+z) + b1*(x+y+z) = 3*R
# (E3 + b1)*S = 3*R
# S = 3*R / (E3 + b1)

### Sum(x[i] * ...) =>
# x*y*z*(x^2 + y^2 + z^2) + b1*E2 = R*S
# E3*(S^2 - 2*E2) + b1*E2 = R*S
# (2*E3 - b1)*E2 = E3*S^2 - R*S
# E2 = (E3*S^2 - R*S) / (2*E3 - b1)

m3.all = unity(3, all=T)

### Example:
b = 1/3
R = 1
#
E3 = - b[1] * m3.all [-1]; # real root has to be removed
S = 3*R / (E3 + b[1]); # E3 != -b1
E2 = (E3*S^2 - R*S) / (2*E3 - b[1])
x = as.vector(sapply(1:2, function(id) roots(c(1, -S[id], E2[id], -E3[id]))))
E3 = rep(E3, each=3)
y = (R - x*E3) / b[1]
z = (R - y*E3) / b[1]
sol = cbind(x,y,z)
sol

### Test
x^2*y*z + b[1]*y
x*y^2*z + b[1]*z
x*y*z^2 + b[1]*x
