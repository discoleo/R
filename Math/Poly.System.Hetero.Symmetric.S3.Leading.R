########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogenous Symmetric
### with Composite Leading Term
###
### draft v.0.1d


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


### draft v.0.1c - v.0.1d:
# - solved:
#   x^2*y + b*x = R; [v.0.1c]
#   x^2*y + b*z = R; [v.0.1d]
# - both solutions are trivial, only for: x == y == z;
### draft v.0.1b - v.0.1b-sp-case:
# - solved: x^2*y + b*y = R;
# - [TODO] factorize P[6]; [DONE]
# - classic Polynomial: P[6]; [v.0.1b-poly]
# - solved special case: R = 0; [v.0.1b-sp-case]
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
(R + b1^2)*S - 3*b1*(R + b1^2) # = 0
(R + b1^2)*(S - 3*b1) # = 0
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

# x^2*y + b*y = R
# y^2*z + b*z = R
# z^2*x + b*x = R

### Solution:

### Sum =>
(x^2*y + y^2*z + z^2*x) + b*S - 3*R # = 0 # Eq 1-bis
(x^2*y + y^2*z + z^2*x)*(x*y^2 + y*z^2 + z*x^2) +
	+ (b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2) # = 0
E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 +
	+ (b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2) # = 0 # Eq 2-bis
### Sum Eq 1-bis + Eq 2-bis =>
(b*S - 3*R)*(x*y^2 + y*z^2 + z*x^2 + x^2*y + y^2*z + z^2*x) +
	+ (b*S - 3*R)^2 + E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 # = 0
(b*S - 3*R)*(E2*S - 3*E3) +
	+ b^2*S^2 - 6*b*R*S + 9*R^2 + E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2 # = 0
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
x*y*z*(x+y+z) + b*E2 - R*S # = 0
E3*S + b*E2 - R*S # = 0

### Sum(y*...) =>
(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*(x^2 + y^2 + z^2) - R*S # = 0
E2^2 - 2*E3*S + b*(S^2 - 2*E2) - R*S # = 0
E2^2 - 2*E3*S + b*S^2 - 2*b*E2 - R*S # = 0

### Auxiliary:
E3Subst = - 27*R*S*b^5 - 12*R*S^3*b^4 + 6*R^2*S^4*b^2 + 9*S^2*b^6 + 5*S^4*b^5;
E3Div = - 27*R*S^2*b^3 - 6*R*S^4*b^2 + 3*S^3*b^4 - S^5*b^3;
# E3 = - E3Subst / E3Div;


### Eq:
((R^2 + b^3)*S^2 - b^2*R*S + b^4) * S^4 * (S^3 + 9*b*S - 27*R) * P[9]
### P[9]: false solution;
(- 6561*R^2*b^3) +
(2916*R*b^4)*S^1 +
(- 6561*R^2*b^2 - 243*b^5)*S^2 +
(2916*R*b^3 - 729*R^3)*S^3 +
(- 972*R^2*b - 189*b^4)*S^4 +
(621*R*b^2)*S^5 +
(- 9*b^3)*S^6 +
(54*R*b)*S^7 + b^2*S^8 + R*S^9


### Solver:
solve.CompositeL.S3P21 = function(R, b, debug=TRUE) {
	if(R[1] == 0) {
		x = sqrt(-b[1] + 0i)
		# the (0, 0, 0) solution is not included;
		return(solve.En(c(x, -x), n=3, duplicates=TRUE))
	}
	coeff = c((R^2 + b[1]^3), - b[1]^2*R, b[1]^4)
	S = roots(coeff)
	if(debug) print(S);
	#
	E3Subst = - 27*R*S*b^3 - 12*R*S^3*b^2 + 6*R^2*S^4 + 9*S^2*b^4 + 5*S^4*b^3;
	E3Div = - 27*R*S^2*b - 6*R*S^4 + 3*S^3*b^2 - S^5*b;
	E3 = - E3Subst / E3Div;
	E2 = (R - E3)*S / b
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	rep.m = function(x) matrix(x, ncol=len, nrow=3, byrow=TRUE)
	S = rep.m(S); E3 = rep.m(E3);
	yz.s = S - x;
	y = R / (x^2 + b);
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z), S=as.vector(S));
	return(sol);
}

### Examples:
R = 2
b = -1
sol = solve.CompositeL.S3P21(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
x^2*y + b*y # - R
y^2*z + b*z # - R
z^2*x + b*x # - R

### Classic Polynomial:
round0.p(poly.calc(x)) * (R^2 + b[1]^3)


#########
### Ex 2:
R = 4
b = 2
sol = solve.CompositeL.S3P21(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Classic Polynomial:
round0.p(poly.calc(x)) * 9 # * (R^2 + b[1]^3)


(2*R^2*b^3 + R^4 + b^6) +
	- R*b^4*x +
	+ b^2*(4*R^2 + 3*b^3)*x^2 +
	- R*(2*b^3 - R^2)*x^3 +
	+ 3*b*(R^2 + b^3)*x^4 +
	- R*b^2*x^5 +
	+ (R^2 + b^3)*x^6


### Debug:
R = 2; b = -1;
x =  1.5907409211 - 0.9236008907i;
y =  0.1489943605 + 0.6462891180i;
z = -1.4064019508 - 0.1940927468i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


#######################

### Variant:
### Mixt-Order: 2+1

### x[i]^2*x[j] + b*x[i]

# x^2*y + b*x = R
# y^2*z + b*y = R
# z^2*x + b*z = R

### Solution:
# - only trivial solution: x == y == z;

### Sum + Rotation =>
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
x*y*z*(x+y+z) + b*E2 - R*S # = 0
E3*S + b*E2 - R*S # = 0

### Sum(y*...) =>
(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*E2 - R*S # = 0
E2^2 - 2*E3*S + b*E2 - R*S # = 0

### Auxiliary:
E3Subst = 6*R*S^3*b^2 - R*S^5*b + 6*R^2*S^4 - S^4*b^3;
E3Div   = - 27*R*S^2*b - 6*R*S^4 - 27*S*b^3 - 6*S^3*b^2;
# E3 = - E3Subst / E3Div;

### Eq: only x == y == z;
(S^3 + 9*b*S - 27*R)^2 * S^3 * P[9]
### P[9]: is false solution
(243*b^9) +
(729*R*b^7)*S^1 +
(729*R^2*b^5 + 135*b^8)*S^2 +
(432*R*b^6 + 243*R^3*b^3)*S^3 +
(486*R^2*b^4 + 18*b^7)*S^4 +
(54*R*b^5 + 216*R^3*b^2)*S^5 +
(54*R^2*b^3 + 27*R^4 + b^6)*S^6 +
(3*R*b^4 + 18*R^3*b)*S^7 +
(3*R^2*b^2)*S^8 +
(R^3)*S^9

### Solver:
solve.CompositeLv.S3P21 = function(R, b, b.ext=0, debug=TRUE) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0)
	if(R[1] == 0) {
		if(b.ext[2] == 0) {
			x = sqrt(-b[1] - 3*b.ext[1] + 0i); x = c(x, -x);
			# the (0, 0, 0) solution is not included;
			# return(solve.En(c(x, -x), n=3, duplicates=TRUE))
		} else {
			x = roots(c(1, 9*b.ext[2], 3*b.ext[1] + b[1]));
		}
		return(cbind(x=x, y=x, z=x, S=3*x));
	}
	coeff = c(1, 0, 9*b[1], - 27*R[1])
	if(any(b.ext != 0)) {
		coeff = coeff + c(0, 27*b.ext[2], 27*b.ext[1], 0);
	}
	S = roots(coeff)
	if(debug) print(S);
	### Note: numerically unstable
	# - result is inaccurate as x == y == z;
	x = S / 3;
	return(cbind(x=x, y=x, z=x, S=S));
	# [unstable]
	R1 = R - b.ext[1]*S - b.ext[2]*S^2;
	E3Subst = 6*R1*S^3*b^2 - R1*S^5*b + 6*R1^2*S^4 - S^4*b^3;
	E3Div = - 27*R1*S^2*b - 6*R1*S^4 - 27*S*b^3 - 6*S^3*b^2;
	E3 = - E3Subst / E3Div;
	E2 = (R1 - E3)*S / b
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	rep.m = function(x) matrix(x, ncol=len, nrow=3, byrow=TRUE)
	S = rep.m(S); E3 = rep.m(E3); R1 = rep.m(R1);
	yz.s = S - x;
	y = (R1 - b*x) / x^2;
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z), S=as.vector(S));
	return(sol);
}

### Examples:
R = 2
b = -1
b.ext = c(1, -1)
sol = solve.CompositeLv.S3P21(R, b, b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
S = x+y+z; ext1 = b.ext[1]*S; ext2 = b.ext[2]*S^2;
x^2*y + b*x + ext1 + ext2 # - R
y^2*z + b*y + ext1 + ext2 # - R
z^2*x + b*z + ext1 + ext2 # - R

### Classic Polynomial:
round0.p(poly.calc(x))


#######################

### Variant 2:
### Mixt-Order: 2+1

### x[i]^2*x[j] + b*x[i]

# x^2*y + b*z = R
# y^2*z + b*x = R
# z^2*x + b*y = R

### Solution:
# - only trivial solution: x == y == z;

### Sum + Rotation =>
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
x*y*z*(x+y+z) + b*(x^2 + y^2 + z^2) - R*S # = 0
E3*S + b*S^2 - 2*b*E2 - R*S # = 0

### Sum(y*...) =>
(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*E2 - R*S # = 0
E2^2 - 2*E3*S + b*E2 - R*S # = 0

### Auxiliary:
E3Subst = 108*R*b^3 - 3*R*S^2*b^2 - 11*R*S^4*b + 6*R^2*S^3 - 36*S*b^4 - 7*S^3*b^3 + 5*S^5*b^2;
E3Div   = 7*S^4*b - 6*R*S^3 - 69*S^2*b^2 + 54*R*S*b + 108*b^3;
# E3 = - E3Subst / E3Div;

### Eq:
(S^3 + 9*b*S - 27*R)^2 * S^3 * P[12]
# P[12]: is a false solution
(15552*b^9) +
(15552*R*b^7)*S^1 +
(7776*R^2*b^5 - 25488*b^8)*S^2 +
(- 23760*R*b^6 + 1944*R^3*b^3)*S^3 +
(- 9072*R^2*b^4 + 18036*b^7)*S^4 +
(13644*R*b^5 - 864*R^3*b^2)*S^5 +
(2880*R^2*b^3 - 27*R^4 - 6667*b^6)*S^6 +
(- 3243*R*b^4 + 171*R^3*b)*S^7 +
(- 357*R^2*b^2 + 1231*b^5)*S^8 +
(310*R*b^3 - R^3)*S^9 +
(3*R^2*b - 97*b^4)*S^10 +
(- 3*R*b^2)*S^11 +
(b^3)*S^12


### Solver:
solve.CompositeLvz.S3P21 = function(R, b, b.ext=0, debug=TRUE) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0)
	if(R[1] == 0) {
		if(b.ext[2] == 0) {
			x = sqrt(-b[1] - 3*b.ext[1] + 0i); x = c(x, -x);
			# the (0, 0, 0) solution is not included;
			# return(solve.En(c(x, -x), n=3, duplicates=TRUE))
		} else {
			x = roots(c(1, 9*b.ext[2], 3*b.ext[1] + b[1]));
		}
		return(cbind(x=x, y=x, z=x, S=3*x));
	}
	coeff = c(1, 0, 9*b[1], - 27*R[1])
	if(any(b.ext != 0)) {
		coeff = coeff + c(0, 27*b.ext[2], 27*b.ext[1], 0);
	}
	S = roots(coeff)
	if(debug) print(S);
	### Note: numerically unstable
	# - result is inaccurate as x == y == z;
	x = S / 3;
	return(cbind(x=x, y=x, z=x, S=S));
	# [unstable]
	R1 = R - b.ext[1]*S - b.ext[2]*S^2;
	E3Subst = 108*R*b^3 - 3*R*S^2*b^2 - 11*R*S^4*b + 6*R^2*S^3 - 36*S*b^4 - 7*S^3*b^3 + 5*S^5*b^2;
	E3Div   = 7*S^4*b - 6*R*S^3 - 69*S^2*b^2 + 54*R*S*b + 108*b^3;
	E3 = - E3Subst / E3Div;
	E2 = (R1 - E3)*S / b
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	rep.m = function(x) matrix(x, ncol=len, nrow=3, byrow=TRUE)
	S = rep.m(S); E3 = rep.m(E3); R1 = rep.m(R1);
	yz.s = S - x;
	y = (R1 - b*x)*x / E3;
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z), S=as.vector(S));
	return(sol);
}

### Examples:
R = 2
b = -1
b.ext = c(0, 0)
sol = solve.CompositeLvz.S3P21(R, b, b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
S = x+y+z; ext1 = b.ext[1]*S; ext2 = b.ext[2]*S^2;
x^2*y + b*z + ext1 + ext2 # - R
y^2*z + b*x + ext1 + ext2 # - R
z^2*x + b*y + ext1 + ext2 # - R

### Classic Polynomial:
round0.p(poly.calc(x))


#######################
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
