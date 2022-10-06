########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogeneous Symmetric
### with Composite Leading Term
###
### draft v.0.2L


### Hetero-Symmetric
### Polynomial Systems: 3 Variables (S3)
### Composite Leading Term (L[n,m])

### Example:
x^n*y^m + P(x, y, z) = R
y^n*z^m + P(y, z, x) = R
z^n*x^m + P(z, x, y) = R


###############
### History ###
###############


### draft v.0.2k:
# - cleaned & updated L11 & L22;
### draft v.0.2h - v.0.2h-clPoly:
# - solved structural extension to S3L33M:
#   x^3*y^3 + b2*x*y*z + b1*x*y = R;
# - Classic Polynomial for the basic S3L33M
#   with A1-Extension; [v.0.2h-clPoly]
### draft v.0.2g - v.0.2g-ext:
# - solved S3L33M system & with simple extensions:
#   x^3*y^3 + b*x*y = R;
# - Eq S & classic Poly for S3L44 Simple:
#   x^4*y^4 + b*z = R;
#   [see more details in the Derivation]
### draft v.0.2e-clPoly3 - v.0.2e-sol-robust:
# - classic Polynomial for:
#   x^3*y*3 + b^z = R;
# - solved: x^3*y*3 + b^z = R;
# - [DONE] robust solution;
### draft v.0.2e - v.0.2e-clPoly:
# - solved: x^2*y^2 + b*z = R;
# - classic Polynomial (for case x == y);
### draft v.0.2d-clean:
# - [cleanup] started moving derivations to file:
#   Poly.System.Hetero.Symmetric.S3.Leading.Derivations.R;
### draft v.0.2c:
# - [started work] Mixed Order 2+2:
#   x^2*y^2 + b*z = R;
### draft v.0.2b - v.0.2b-S11:
# - [started work] Mixed Order 3+1:
#   x^3*y + b*z = R;
# - [DONE] S11;
### draft v.0.2a:
# - solved (part):
#   x^2*y + a*x*y^2 = R;
# - TODO: solve Case 2 & clean solution;
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

### Helper Functions


source("Polynomials.Helper.R")


# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions

test.CHP.S3.Symmetric = function(sol, b, b.ext=0, R=NULL, a=0, n=1) {
	test.S3Ht.LnnChY(sol, b=b, b.ext=b.ext, R=R, a=a, n=n);
}
test.S3Ht.LnnChY = function(sol, b, b.ext=0, R=NULL, a=0, n=1) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0);
	if(length(a) < 2) a = c(a, 0)
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	xyz = x*y*z; a.ext = a[1]*xyz + a[2]*xyz^2;
	s = (x+y+z); s.ext = b.ext[1]*s + b.ext[2]*s^2;
	ext = a.ext + s.ext;
	err1 = (x*y)^n + b[1]*y + ext # - R
	err2 = (y*z)^n + b[1]*z + ext # - R
	err3 = (z*x)^n + b[1]*x + ext # - R
	round0(rbind(err1, err2, err3))
}
test.S3Ht.LSymmetricChz = function(sol, b, b.ext=0, R=NULL, n=1) {
	test.S3Ht.LnnChZ(sol, b=b, b.ext=b.ext, R=R, n=n)
}
test.S3Ht.LnnChZ = function(sol, b, b.ext=0, R=NULL, n=1) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0);
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	S = (x+y+z); s.ext = b.ext[1]*S + b.ext[2]*S^2;
	ext = s.ext;
	err1 = (x*y)^n + b[1]*z + ext # - R
	err2 = (y*z)^n + b[1]*x + ext # - R
	err3 = (z*x)^n + b[1]*y + ext # - R
	round0(rbind(err1, err2, err3));
}

### Classic Polynomial

classicPoly.Ln1z = function(n, stats=TRUE) {
	p1 = toPoly.pm("x^n*y + b*z - R");
	p2 = toPoly.pm("y^n*z + b*x - R");
	p3 = toPoly.pm("z^n*x + b*y - R");
	
	# inefficient:
	# pR = solve.lpm(p1, p2, p3, xn=c("z", "y"))
	pR = solve.lpm(p1, p2, p3, xn=c("y", "z"));
	if(stats) {
		print(str(pR[[2]]));
		print(table(pR[[2]]$Rez$x))
	}
	#
	n = n + 1;
	pR2 = div.pm(pR[[2]]$Rez, "x^n + b*x - R", "x");
	invisible(pR2$Rez);
}

################################
################################

########################
### Leading: (X*Y)^n ###
########################

###############
### Order 1 ###
###############

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

### Diff & Prod =>
# E3 = b1^3

### Auxiliary Eqs:
# E3 = b1^3
# b1*E2 = R*S - 3*E3

### Eq:
(R + b1^2)*(S - 3*b1) # = 0
# Note:
# S = 3*b1 is a FALSE solution;
# R + b1^2: can have solutions
# (in the extensions, e.g.: be*S = R + b1^2);


### Solver:
solve.S3Ht.L11 = function(R, b, b.ext=0, a=0, debug=TRUE) {
	be1 = b.ext[1];
	be2 = if(length(b.ext) > 1) b.ext[2] else 0;
	a1 = a[1];
	a2 = if(length(a) > 1) a[2] else 0;
	if(be1 == 0 && be2 == 0) {
		stop("NO solutions: x != y != z")
		S = 3*b[1];
	} else {
		S = roots(c(be2, be1, -R[1] - b[1]^2 + a1*b[1]^3 + a2*b[1]^6))
	}
	if(debug) print(S);
	len = length(S)
	E3 = b[1]^3 - 0*S;
	R1 = R[1] - be1*S - be2*S^2 - a1*E3 - a2*E3^2;
	E2 = 3*R1 - b[1]*S;
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	#
	R1 = rep(R1, each=3);
	S  = rep(S, each=3);
	E3 = rep(E3, each=3);
	### robust
	yz = E3 / x;
	z = (R1 - yz) / b[1];
	y = S - x - z;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol)
}
test.S3Ht.L11 = function(sol, b, b.ext=0, R=NULL, a=0) {
	test.S3Ht.LnnChY(sol, b=b, b.ext=b.ext, R=R, a=a, n=1);
}

### Examples:
R = -1
b = 3
b.ext = c(1, 0)
sol = solve.S3Ht.L11(R, b, b.ext=b.ext)

### Test
test.S3Ht.L11(sol, b, b.ext=b.ext)


### Ex 2:
R = -1
b = 3
b.ext = c(1, 3)
sol = solve.S3Ht.L11(R, b, b.ext=b.ext)

### Test
test.S3Ht.L11(sol, b, b.ext=b.ext)
round0.p(poly.calc(sol[ , 1]))


############
### Str Ext: Ex 1
R = -1
b = 3
b.ext = c(1, 3); a = 2;
sol = solve.S3Ht.L11(R, b, b.ext=b.ext, a=a)

### Test
test.S3Ht.L11(sol, b, b.ext=b.ext, a=a)
round0.p(poly.calc(sol[ , 1]))


############
### Str Ext: Ex 2
R = -1
b = 3
b.ext = c(1, 3); a = c(-1, 1)
sol = solve.S3Ht.L11(R, b, b.ext=b.ext, a=a)

### Test
test.S3Ht.L11(sol, b, b.ext=b.ext, a=a)
round0.p(poly.calc(sol[ , 1]))


############
### Str Ext: Ex 3
R = -1
b = 3
b.ext = c(0, -1); a = c(2, -1)
sol = solve.S3Ht.L11(R, b, b.ext=b.ext, a=a)

### Test
test.S3Ht.L11(sol, b, b.ext=b.ext, a=a)
round0.p(poly.calc(sol[ , 1]))


#################

### Classic Poly:
# - without extensions;
(x^2 + b*x - R)

### FALSE roots:
# - but distinct solutions possible under very special conditions;
(R + b^2)

### Classic Solver:
# - roots possible only under very special conditions;

### Special Cases:
solve.S3L11y.Special = function(R, b) {
	# R = 2*b^2; # Case: x=y=z = b, all equal!
	spec = round0(c(R + b^2, R - 2*b^2));
	if( ! any(spec == 0)) {
		return(NA);
	}
	x = b; y = R / (2*b); z = (R - b^2) / b;
	sol = cbind(x, y, z);
	return(sol);
}

### Ex 1: distinct;
b = -1; R = -b^2;
sol = solve.S3L11y.Special(R, b=b);

print(R);
test.S3Ht.L11(sol, b)

### Ex 2: distinct;
b = 3; R = -b^2;
sol = solve.S3L11y.Special(R, b=b);

print(R);
test.S3Ht.L11(sol, b)


##############
##############

### Variant z:

x*y + b*z # = R
y*z + b*x # = R
z*x + b*y # = R

### Solution:

### Eq 1:
# Sum =>
E2 + b*S - 3*R # = 0

### Eq 2:
# Sum(z*...) =>
3*E3 + b*(S^2 - 2*E2) - R*S # = 0

### Eq 3:
# Sum(z^2*...) =>
E3*S + b*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0

### Eq S:
(S^2 + 3*b*S - 9*R) * (b*S - R - b^2)

### Solver:

solve.S3Ht.L11z = function(R, b, debug=TRUE) {
	S = (R + b^2) / b;
	if(debug) print(S);
	#
	E2 = 3*R - b*S; # = 2*R - b^2;
	E3 = (2*b*E2 - b*S^2 + R*S) / 3; # = b*(R - b^2);
	# Robust:
	# x = roots(c(1, -S, E2, -E3)); # non-robust;
	# x = roots(c(b, -R, E3));
	# - the Quadratic has 2 roots,
	#   but x = b is unstable!
	x = (R - b^2) / b;
	y = (x - b)*R / (x^2 - b^2);
	z = (R - x*y) / b;
	sol = cbind(x, y, z);
	### x = b
	# - are only rotations of the previous root;
	x = b;
	s = R / b; e2 = R - b^2;
	y = roots(c(1, -s, e2));
	z = s - y;
	#
	sol = rbind(sol, cbind(x, y, z));
	return(sol);
}
test.S3Ht.L11z = function(sol, b, R=NULL) {
	err = test.S3Ht.LnnChZ(sol, b=b, R=R, n=1);
	return(err)
}

### Examples:

### Ex 1:
R = -1; b = 3;
sol = solve.S3Ht.L11z(R, b=b)

test.S3Ht.L11z(sol, b)


### Ex 2:
R = -3; b = -4;
sol = solve.S3Ht.L11z(R, b=b)

test.S3Ht.L11z(sol, b)


### Classic Poly:
# - without extensions;
# - Note: (b*x - R) is a FALSE solution;
(x^2 + b*x - R) *
(x - b)^2 * (b*x - R + b^2) * (b*x - R)

### Classic Solver:
solve.S3.L11z.cl = function(R, b) {
	coeff = c(b, - R + b^2);
	x = roots(coeff);
	### x != b; avoid Error during division;
	y = (x - b)*R / (x^2 - b^2);
	z = (R - x*y) / b;
	sol = cbind(x, y, z);
	### x = b
	x = b;
	s = R / b; e2 = R - b^2;
	y = roots(c(1, -s, e2));
	z = s - y;
	sol = rbind(sol, cbind(x, y, z));
	return(sol);
}

###
R = -1; b = 3;
sol = solve.S3.L11z.cl(R, b=b)

test.S3Ht.L11z(sol, b)


########################
########################

########################
### Mixed-Order: 2+2 ###
########################

### x[i]^2*x[j]^2 + b*x[k]

# x^2*y^2 + b*z = R
# y^2*z^2 + b*x = R
# z^2*x^2 + b*y = R

### Solution:

### Case: (x, y, z) distinct;
# - NO solutions;

### Case: x == y, z distinct;

### Sum =>
E2^2 - 2*E3*S + b*S - 3*R # = 0

### Sum(z*...) =>
E2*E3 + b*(S^2 - 2*E2) - R*S # = 0

### Sum(z^2*...) =>
3*E3^2 + 3*b*E3 - 3*b*E2*S + 2*R*E2 + b*S^3 - R*S^2 # = 0
# Reduction =>
3*E3^2 - E2*E3*S + 3*b*E3 - b*E2*S + 2*R*E2 # = 0


### Note:
# Invalid Solution: Diff & Sum =>
E2*S - 3*E3 - 3*b # = 0


### Eq S:
b^2*S^6 - 2*b*R*S^5 + R^2*S^4 - 9*b^3*S^3 + 9*b^2*R*S^2 - 3*b*R^2*S + 27*b^4 - R^3
# alternatively: solve directly for x
x^6 - b*x^3 - R*x^2 + b^2;
# S = 2*x + (R - x^4)/b;

### Note: Case S = 0
# - seems to have NO solutions (except in special cases);


### Solver:

# Robust:
solve.S3L22.Simple = function(R, b, be=0, debug=TRUE) {
	coeff = coeff.S3L22.Simple(R, b, be=be);
	S = roots(coeff);
	if(debug) print(S);
	#
	R  = R - be[1]*S; # Extension
	#
	E2x0 = 2*b*R*S^5 - 2*R^2*S^4 - 9*b^3*S^3 + 27*b^2*R*S^2 + 36*b*R^2*S - 243*b^4 - 18*R^3;
	E2Div = 30*b^2*S^4 - 48*b*R*S^3 + 20*R^2*S^2 - 189*b^3*S + 81*b^2*R;
	E2 = E2x0 / E2Div;
	E3 = - (b*S^2 - 2*b*E2 - R*S) / E2;
	# [alternative]
	# E3 = ...;
	# but: S = 0 => E3 = 2*b;
	# E2 = (b*S^2 - R*S) / (2*b - E3);
	### Case: y != z
	x = (E3 + b) / E2;
	s = S - x; e2 = E2 - x*s;
	len = length(s);
	yz  = sapply(seq(len), function(id) {
		roots(c(1, -s[id], e2[id]));
	});
	yz = t(yz);
	y = yz[,1]; z = yz[,2];
	sol = cbind(x, y, z);
	return(sol);
}
coeff.S3L22.Simple = function(R, b, be=0) {
	coeff = c(b^2, - 2*b*R, R^2, - 9*b^3, 9*b^2*R, - 3*b*R^2, 27*b^4 - R^3);
	if(any(be != 0)) {
		coeff = coeff +
			c(2*b*be[1] + be[1]^2, -2*R*be[1], 0, be[1]^3 - 3*b*be[1]^2 - 9*b^2*be[1],
				-3*R*be[1]^2 + 6*b*R*be[1], 3*R^2*be[1], 0);
	}
	return(coeff);
}
test.S3Ht.L22z = function(sol, b, be=0, R=NULL) {
	err = test.S3Ht.LnnChZ(sol, b=b, b.ext=be, R=R, n=2);
	return(err)
}

### Examples;
R = -1;
b = 3;
be = 0;
sol = solve.S3L22.Simple(R, b, be=be);
x = sol[,1]; y = sol[,2]; z = sol[,3];
S = (x+y+z);

### Test:
test.S3Ht.L22z(sol, b=b)


### Extensions:
### Ex 2:
R = -1;
b = 3;
be = 2;
sol = solve.S3L22.Simple(R, b, be=be);

test.S3Ht.L22z(sol, b=b, be=be)


### Ex 3: S = 0
b = 8;
R = 3*b^(4/3);
be = 0;
sol = solve.S3L22.Simple(R, b, be=be);

print(R)
test.S3Ht.L22z(sol, b=b, be=be)


### Classic Polynomial:
# - for case x == y:
x^10 - 2*R*x^6 + R^2*x^2 + b^3*x - b^2*R
# - for case y == z:
(b^2*x^10 - 2*b*R*x^9 + R^2*x^8 + 2*b*R^2*x^5 - 2*R^3*x^4 + 4*b^3*R*x^3 +
	- 4*b^2*R^2*x^2 + b^5*x - b^4*R + R^4)
# (but one has to know this)
(x^4 + b*x - R) * (x^6 - b*x^3 - R*x^2 + b^2) *
	(b^2*x^6 - 2*b*R*x^5 + R^2*x^4 - b^3*x^3 + 3*b^2*R*x^2 - b*R^2*x + b^4 - R^3)



########################
########################

########################
### Mixed-Order: 3+3 ###
########################

### x[i]^3*x[j]^3 + b*x[k]

# x^3*y^3 + b*z = R
# y^3*z^3 + b*x = R
# z^3*x^3 + b*y = R

### Solution:

### Sum =>
E2^3 - 3*E3*E2*S + 3*E3^2 + b*S - 3*R

### Sum(z*...) =>
2*E3^2*S - E2^2*E3 + 2*b*E2 - b*S^2 + R*S # = 0

### Diff =>
# - if distinct solutions exist?
E2^2 - E3*S # = 0

### Eq S:
b^3*S^5 - 3*b^2*R*S^4 + 3*b*R^2*S^3 - R^3*S^2 + b^4

### Solver:
solve.S3Ht.L33ChZ = function(R, b, be=0, debug=TRUE) {
	coeff = coeff.S3Ht.L33ChZ(R, b, be=be);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2x0 = (32*b^2*S^7 - 56*b*R*S^6 + 24*R^2*S^5 + 144*b^3*S^2 - 216*b^2*R*S);
	E2Div = (104*b^2*S^5 - 168*b*R*S^4 + 72*R^2*S^3 + 216*b^3);
	E2 = E2x0 / E2Div;
	E3 = E2^2 / S;
	# Step 2:
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	# robust:
	# x^3*yz*y^2 + b*z^2 - R*z = 0
	# x^3*yz*y^2 + b*(yz.s - y)^2 - R*(yz.s - y) = 0
	x = as.vector(x); yz.s = as.vector(yz.s);
	y = sapply(seq_along(x), function(id) {
		roots(c(x[id]^3*yz[id] + b, -2*b*yz.s[id] + R, b*yz.s[id]^2 - R*yz.s[id]));
	});
	x = rep(x, each=2); yz.s = rep(yz.s, each=2);
	y = as.vector(y);
	z = yz.s - y;
	sol = cbind(x=x, y=y, z=z);
	# TODO:
	# add also the roots with pairwise equal variables;
	return(sol);
}
coeff.S3Ht.L33ChZ = function(R, b, be=0) {
	coeff = c(b^3, - 3*b^2*R, 3*b*R^2, - R^3, 0, b^4);
	return(coeff);
}
test.S3Ht.L33 = function(sol, b, be=0, R=NULL) {
	test.S3Ht.LnnChZ(sol, b=b, b.ext=be, R=R, n=3)
}

### Examples:
R = -1;
b = 3;
sol = solve.S3Ht.L33ChZ(R, b);

test.S3Ht.L33(sol, b)


### Ex 2:
R = 2;
b = -3;
sol = solve.S3Ht.L33ChZ(R, b);

test.S3Ht.L33(sol, b)


### Test
x = sol[,1]; y = sol[,2]; z = sol[,3];
x^3*y^3 + b*z # - R
y^3*z^3 + b*x # - R
z^3*x^3 + b*y # - R


### Classic Polynomial:
# P[51] = P[6] * P[15] * P[15] * P[15];
### Case x == y:
x^21 - 3*R*x^15 + 3*R^2*x^9 - R^3*x^3 - b^4*x + b^3*R
### Case y == z:
b^3*x^21 - 3*b^2*R*x^20 + 3*b*R^2*x^19 - R^3*x^18 + 3*b^2*R^2*x^14 - 6*b*R^3*x^13 + 3*R^4*x^12 +
	- 2*b^5*x^11 + 4*b^4*R*x^10 - 2*b^3*R^2*x^9 + 3*b*R^4*x^7 - 3*R^5*x^6 + 6*b^4*R^2*x^4 +
	- 6*b^3*R^3*x^3 + b^7*x - b^6*R + R^6
# (but one has to know this)
### Case: all distinct:
b^3*x^15 - 3*b^2*R*x^14 + 3*b*R^2*x^13 - R^3*x^12 + b^4*x^10 - b^3*R*x^9 +
	- b*R^3*x^7 + R^4*x^6 + b^5*x^5 + b^4*R*x^4 - 2*b^3*R^2*x^3 + b^6
### All Cases:
(x^6 + b*x - R) * (x^15 - b*x^10 - 2*R*x^9 + b^2*x^5 + b*R*x^4 + R^2*x^3 - b^3) *
(b^3*x^15 - 3*b^2*R*x^14 + 3*b*R^2*x^13 - R^3*x^12 - b^4*x^10 + 4*b^3*R*x^9 - 3*b^2*R^2*x^8 +
	- 2*b*R^3*x^7 + 2*R^4*x^6 - b^5*x^5 - b^4*R*x^4 + 5*b^3*R^2*x^3 - b^2*R^3*x^2 - b*R^4*x + b^6 - R^5)
# * P[15][all distinct]


########################
########################

########################
### Mixed-Order: 4+4 ###
########################

### x^4*y^4 + b*z

# x^4*y^4 + b*z = R
# y^4*z^4 + b*x = R
# z^4*x^4 + b*y = R

### Solution:

### Case: (x,y,z) all distinct

### Sum =>
4*E2*E3^2 + 2*E3^2*S^2 - 4*E2^2*E3*S + E2^4 + b*S - 3*R # = 0
# Reduction (Diff-based!) =>
3*E2*E3^2 + 2*E3^2*S^2 - 2*E2^2*E3*S + b*S - 3*R # = 0

### Sum(z*...) =>
3*E3^3 - 3*E2*E3^2*S + E2^3*E3 - 2*b*E2 + b*S^2 - R*S # = 0
# Reduction (Diff-based!) =>
3*E2*E3^2*S - 2*E2^3*E3 - 2*b*E2 + b*S^2 - R*S # = 0
2*E3^2*S^3 - 2*E2^2*E3*S^2 + 2*E2^3*E3 + 2*b*E2 - 2*R*S # = 0

### Diff =>
E3^2 - 2*E2*E3*S + E2^3 # = 0


### Eq S:
# P[14] (true roots):
-R^7 - b^8 + 9*b*R^6*S - 39*b^2*R^5*S^2 + 103*b^3*R^4*S^3 - 175*b^4*R^3*S^4 + 187*b^5*R^2*S^5 +
	- 113*b^6*R*S^6 + 29*b^7*S^7 + R^6*S^8 - 6*b*R^5*S^9 + 15*b^2*R^4*S^10 - 20*b^3*R^3*S^11 +
	+ 15*b^4*R^2*S^12 - 6*b^5*R*S^13 + b^6*S^14


### Auxiliary Eqs:
# TODO


### Solver:

solve.S3Ht.L44ChZ = function(R, b, be=0, debug=TRUE) {
	coeff = coeff.S3Ht.L44ChZ(R, b, be=be);
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	E2 = sapply(seq(len), function(id) {
		E2 = Ex.S3Ht.L44ChZ(S[id], R, b);
		E2 = - E2[1];
		return(E2);
	});
	E3 = (3*E2^4 + 2*E2^3*S^2 - b*S + 3*R) / (4*E2*S*(S^2 + E2));
	if(debug) { cat("E2:\n"); print(E2); cat("E3:\n"); print(E3); }
	# Step 2:
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	x = as.vector(x);
	S = rep(S, each=3); E3 = rep(E3, each=3);
	s = S - x; yz = E3 / x;
	# robust:
	len = length(x);
	y = sapply(seq(len), function(id) {
		roots(c(1, -s[id], yz[id]));
	})
	y = as.vector(y);
	x = rep(x, each=2); s = rep(s, each=2);
	z = s - y;
	sol = cbind(x=x, y=y, z=z);
	# TODO:
	# add also the roots with pairwise equal variables;
	return(sol);
}
coeff.S3Ht.L44ChZ = function(R, b, be=0) {
	coeff = c(b^6, - 6*b^5*R, 15*b^4*R^2, - 20*b^3*R^3, 15*b^2*R^4,
		- 6*b*R^5, R^6, 29*b^7, - 113*b^6*R, 187*b^5*R^2, - 175*b^4*R^3,
		103*b^3*R^4, - 39*b^2*R^5,	9*b*R^6, - b^8 - R^7);
	return(coeff);
}
Ex.S3Ht.L44ChZ = function(S, R, b) {
	cc = c((18*R^2*S^2 - 12*R*S^3*b + 2*S^4*b^2) / 27,
		(27*R^2 - 18*R*S*b + 3*S^2*b^2) / 27,
		(- 48*R + 16*b*S)*S^6 / 27,
		(- 96*R + 32*b*S)*S^4 / 27, 0,
		(54*R - 18*S*b) / 27,
		- 8*S^6 / 27, - 4*S^4 / 27, 30*S^2 / 27, 1);
	# [excluded] Case: S = 0
	dd = c((27*R^2 - 18*b*R*S + 3*b^2*S^2) / 3,
		16*(S*b - R)*S^6 / 3,
		- 32*R*S^4 / 3, - (4*R*S^2 + 52*b*S^3) / 3,
		(30*R - 42*b*S) / 3, 0, - 4*S^4 / 3, - 4*S^2 / 3, 1);
	cc = cc - c(0, dd);
	cc = round0(cc); dd = round0(dd);
	len = length(cc);
	strip = function(x) {
		pos = 0;
		for(i in seq(length(x), 1, by=-1)) {
			if(x[i] == 0) { pos = i; }
			else break;
		}
		if(pos > 0) x = x[seq(pos - 1)];
		x = x / x[length(x)];
		return(x);
	}
	for(i in seq(len - 1)) {
		# if(debug) print(paste0("Step: ", i)); # Debug
		cc = strip(cc);
		dd = strip(dd);
		len1 = length(cc); len2 = length(dd);
		if(len2 <= 2) return(dd);
		while(len1 >= len2) {
			cc = cc - c(rep(0, len1 - len2), dd);
			cc = cc[ - len1];
			cc = strip(cc); len1 = length(cc);
		}
		if(len1 <= 2) return(cc);
		while(len2 >= len1) {
			dd = dd - c(rep(0, len2 - len1), cc);
			dd = dd[ - len2];
			dd = strip(dd); len2 = length(dd);
		}
		if(len2 <= 2) return(dd);
	}
	warning("Some error!");
	return(NA);
}
test.S3Ht.L44 = function(sol, b, be=0, R=NULL) {
	test.S3Ht.LnnChZ(sol, b=b, b.ext=be, R=R, n=4)
}

### Examples:
R = -1;
b = 3;
sol = solve.S3Ht.L44ChZ(R, b);

test.S3Ht.L44(sol, b)


### Ex:
R = 2;
b = -3;
sol = solve.S3Ht.L44ChZ(R, b);

test.S3Ht.L44(sol, b)


### Classic Polynomial:

### Case: (x,y,z) all distinct
b^12 - 3*b^8*R^3*x^4 + b^9*R^2*x^5 + b^10*R*x^6 + b^11*x^7 + 3*b^4*R^6*x^8 - 2*b^5*R^5*x^9 - b^6*R^4*x^10 +
	- (4*b^8*R^2 + R^9)*x^12 + b*R*(4*b^8 + R^7)*x^13 - b^4*R^5*x^16 + 5*b^5*R^4*x^17 - 2*b^6*R^3*x^18 +
	- 6*b^7*R^2*x^19 + (4*b^8*R+3*R^8)*x^20 - 8*b*R^7*x^21 + 6*b^2*R^6*x^22 - 5*b^4*R^4*x^24 +
	+ 12*b^5*R^3*x^25 - 12*b^6*R^2*x^26 + 4*b^7*R*x^27 - 3*R^7*x^28 + 13*b*R^6*x^29 - 21*b^2*R^5*x^30 +
	+ 15*b^3*R^4*x^31 - 5*b^4*R^3*x^32 + 3*b^5*R^2*x^33 - 3*b^6*R*x^34 + b^7*x^35 + R^6*x^36 - 6*b*R^5*x^37 +
	+ 15*b^2*R^4*x^38 - 20*b^3*R^3*x^39 + 15*b^4*R^2*x^40 - 6*b^5*R*x^41 + b^6*x^42

### Case: x == y
# - see file:
#   Poly.System.Hetero.Symmetric.S3.Leading.Derivation.R;


#########################
#########################

#########################
### Variants:         ###
### (x*y)^n + b*(x*y) ###
#########################

###############
### Order 3 ###
###############

# x^3*y^3 + b*x*y = R
# y^3*z^3 + b*y*z = R
# z^3*x^3 + b*z*x = R

### Solution:

### Case: (x,y,z) all distinct

### Diff =>
E2 # = 0

### Sum =>
E3^2 - R # = 0

### Sum(z^3*...) =>
3*E3^3 - 2*b*E2*E3 + b*E3*S^2 - 3*R*E3 + 3*R*E2*S - R*S^3 # = 0

### Eq S:
R*S^2 - b^2 # = 0
# S = 0 is NOT a solution for the base-system;

### Auxiliary Eqs:
# E2 = 0
# b*E3 = R*S


### Solver:
solve.S3L33M = function(R, b, be=0, all.sol="Eq2", debug=TRUE) {
	sol.type = pmatch(all.sol, c("Distinct", "Eq2", "All"));
	coeff = c(R, 0, - b^2);
	if(any(be != 0)) {
		coeff = c(- rev(be), coeff);
	}
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	# Extensions:
	R1 = R[1] - sapply(S, function(S) sum(be*S^seq(length(be))));
	E2 = rep(0, len);
	E3 = R1*S / b;
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	### robust
	yz.d = sqrt(yz.s^2 - 4*yz + 0i);
	y = (yz.s + yz.d)/2;
	z = (yz.s - yz.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = rbind(sol, sol[,c(1,3,2)]);
	if(sol.type > 1) {
		sol = rbind(sol, solve.S3L33MEq(R, b, be=be));
	}
	return(sol);
}
solve.S3L33MEq = function(R, b, be=0) {
	if(all(be == 0)) {
		coeff = c(1, 0,0,0, b, 0, -R[1]);
		x = roots(coeff);
		z = sapply(x, function(x) roots(c(x^3, 0, b*x, -R[1])));
		x = rep(x, each=3);
		return(cbind(x=x, y=x, z=as.vector(z)));
	}
	# TODO
}

### Examples:
R = -3
b = -1
#
sol = solve.S3L33M(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
x^3*y^3 + b*x*y # - R
y^3*z^3 + b*y*z # - R
z^3*x^3 + b*z*x # - R


### Extensions:
R = -2
b = -2
be = 1/3
#
sol = solve.S3L33M(R, b, be=be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
S = (x+y+z); ext = be[1]*S;
x^3*y^3 + b*x*y + ext # - R
y^3*z^3 + b*y*z + ext # - R
z^3*x^3 + b*z*x + ext # - R

### Classic Poly:
# with Extension: + be1*(x+y+z)
be1 = be[1]; b = b[1];
be1*x^9 - R*x^8 + b^2*x^6 - 3*b*be1*x^5 + 2*b*R*x^4 - be1*R*x^3 + R^2*x^2 + b*be1^2

round0.p(poly.calc(x[1:9]) * 3)
# * be[1] for integer be[1];

# Special Case: be1 = 1/b =>
x^9 - b*R*x^8 + b^3*x^6 - 3*b*x^5 + 2*b^2*R*x^4 - R*x^3 + b*R^2*x^2 + 1


##########
### Ext 2:
R = -2
b = -2
be = c(-1/3, 1/3)
#
sol = solve.S3L33M(R, b, be=be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
S = (x+y+z); ext = be[1]*S + be[2]*S^2;
x^3*y^3 + b*x*y + ext # - R
y^3*z^3 + b*y*z + ext # - R
z^3*x^3 + b*z*x + ext # - R


### Classic Poly:
# TODO

round0.p(poly.calc(x[1:12]) * 3)


##########################

### Structural Extensions:
### + b*x*y*z

# x^3*y^3 + b2*x*y*z + b1*x*y = R
# y^3*z^3 + b2*x*y*z + b1*y*z = R
# z^3*x^3 + b2*x*y*z + b1*z*x = R

### Solution:

### Case: (x,y,z) all distinct

### Sum =>
E3^2 + b2*E3 - R # = 0

### Sum(z^3*...) =>
b2*E3*S^3 + b1*E3*S^2 - R*S^3 # = 0

### Auxiliary Eqs:
# (b2*S + b1)*E3 = R*S;

### Eq S:
R*S^2 - b1*b2*S - b1^2 # = 0


### Solver:
solve.S3L33M_Ext = function(R, b, be=0, all.sol="Distinct", debug=TRUE) {
	sol.type = pmatch(all.sol, c("Distinct", "Eq2", "All"));
	coeff = c(R, -b[1]*b[2], - b[1]^2);
	if(any(be != 0)) {
		coeff = c(- rev(be), coeff);
	}
	S = roots(coeff);
	if(debug) print(S);
	len = length(S);
	# Extensions:
	R1 = R[1] - sapply(S, function(S) sum(be*S^seq(length(be))));
	E2 = rep(0, len);
	E3 = R1*S / (b[2]*S + b[1]);
	R1 = R1 - b[2]*E3; # [not needed]
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	### robust
	yz.d = sqrt(yz.s^2 - 4*yz + 0i);
	y = (yz.s + yz.d)/2;
	z = (yz.s - yz.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = rbind(sol, sol[,c(1,3,2)]);
	if(sol.type > 1) {
		sol = rbind(sol, solve.S3L33MEq(R, b, be=be));
	}
	return(sol);
}

### Examples:
R = -3
b = c(-1, -2)
#
sol = solve.S3L33M_Ext(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
x^3*y^3 + b[2]*x*y*z + b[1]*x*y # - R
y^3*z^3 + b[2]*x*y*z + b[1]*y*z # - R
z^3*x^3 + b[2]*x*y*z + b[1]*z*x # - R


### Extensions:
R = -3
b = c(-1, -2)
be = 1/3
#
sol = solve.S3L33M_Ext(R, b, be=be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
S = (x+y+z); ext = be[1]*S;
x^3*y^3 + b[2]*x*y*z + b[1]*x*y + ext # - R
y^3*z^3 + b[2]*x*y*z + b[1]*y*z + ext # - R
z^3*x^3 + b[2]*x*y*z + b[1]*z*x + ext # - R


### Ext 2:
R = -3
b = c(-1, -2)
be = c(-2, 1/3)
#
sol = solve.S3L33M_Ext(R, b, be=be)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test:
S = (x+y+z); ext = be[1]*S + be[2]*S^2;
x^3*y^3 + b[2]*x*y*z + b[1]*x*y + ext # - R
y^3*z^3 + b[2]*x*y*z + b[1]*y*z + ext # - R
z^3*x^3 + b[2]*x*y*z + b[1]*z*x + ext # - R


########################
########################
########################

########################
### Mixed-Order: 2+1 ###
########################

### x[i]^2*x[j] + b*x[j]

# x^2*y + b*y = R
# y^2*z + b*z = R
# z^2*x + b*x = R

### Solution:

# - shortcut: using eq. from the HtSymmetric Mixed system;
#   (see x^2*y + y^2*z + z^2*x = R1)

### Sum =>
# shortcut: Eq1[Mixed](R1 = 3*R - b*S) =>
E3*S^3 - ((3*R - b*S) + 6*E3)*E2*S + (3*R - b*S)^2 + E2^3 + 9*E3^2 + 3*(3*R - b*S)*E3 # = 0
# long (see Derivation):
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
E3*S + b*E2 - R*S # = 0

### Sum(y*...) =>
E2^2 - 2*E3*S + b*S^2 - 2*b*E2 - R*S # = 0

### Auxiliary:
E3Subst = - 27*R*S*b^5 - 12*R*S^3*b^4 + 6*R^2*S^4*b^2 + 9*S^2*b^6 + 5*S^4*b^5;
E3Div = - 27*R*S^2*b^3 - 6*R*S^4*b^2 + 3*S^3*b^4 - S^5*b^3;
# E3 = - E3Subst / E3Div;


### Eq S:
((R^2 + b^3)*S^2 - b^2*R*S + b^4) * S^4 * (S^3 + 9*b*S - 27*R) * P[9]
### P[9]: false solution;


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
round0.p(poly.calc(x)) * 3 # * (R^2 + b[1]^3)


(2*R^2*b^3 + R^4 + b^6) +
	- R*b^4*x +
	+ b^2*(4*R^2 + 3*b^3)*x^2 +
	- R*(2*b^3 - R^2)*x^3 +
	+ 3*b*(R^2 + b^3)*x^4 +
	- R*b^2*x^5 +
	+ (R^2 + b^3)*x^6


#######################

### Variant:
### Mixed-Order: 2+1

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
#######################

### Variant 2:
### Mixed-Order: 2+1

### x[i]^2*x[i+1] + b*x[i+2]

# x^2*y + b*z = R
# y^2*z + b*x = R
# z^2*x + b*y = R

### Solution:
# - P[12] & Trivial solution: x == y == z;
# - Shortcut: using Eq. E21a from the HtSymmetric Mixed system;

### Eq 1:
# Sum + Rotation =>
E3*S^3 - 3*b*E3*S + 9*E3^2 + 9*R*E3 - 6*E3*E2*S +
	+ E2^3 + b*E2*S^2 - 3*R*E2*S + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Eq 2:
# Sum(z*...) =>
x*y*z*(x+y+z) + b*(x^2 + y^2 + z^2) - R*S # = 0
E3*S + b*S^2 - 2*b*E2 - R*S # = 0

### Eq 3:
# Sum(y*...) =>
(x^2*y^2 + y^2*z^2 + z^2*x^2) + b*E2 - R*S # = 0
E2^2 - 2*E3*S + b*E2 - R*S # = 0

### Auxiliary:
E3Subst = 108*R*b^3 - 3*R*S^2*b^2 - 11*R*S^4*b + 6*R^2*S^3 - 36*S*b^4 - 7*S^3*b^3 + 5*S^5*b^2;
E3Div   = 7*S^4*b - 6*R*S^3 - 69*S^2*b^2 + 54*R*S*b + 108*b^3;
# E3 = - E3Subst / E3Div;

### Eq S: P[4]
(b^2*S^4 - 2*b*R*S^3 + (R^2 - 3*b^3)*S^2 + 2*b^2*R*S + 4*b^4)
# * (S^3 + 9*b*S - 27*R)^2 * S^3 * P[8];

# P[8]: is a FALSE solution
b*S^8 - R*S^7 - 94*b^2*S^6 + 117*b*R*S^5 + 945*b^3*S^4 - 27*R^2*S^4 - 810*b^2*R*S^3 +
	- 3456*b^4*S^2 + 1944*b^3*R*S + 3888*b^5 # FALSE


### Solver:
solve.S3Ht.L21z = function(R, b, b.ext=0, debug=TRUE) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0);
	# if(R[1] == 0) {}
	coeff = coeff.S3Ht.L21z(R, b, b.ext);
	S = roots(coeff);
	if(debug) print(S);
	# Step 2:
	R = R - b.ext[1]*S - b.ext[2]*S^2;
	E3Subst = 108*R*b^3 - 3*R*S^2*b^2 - 11*R*S^4*b + 6*R^2*S^3 - 36*S*b^4 - 7*S^3*b^3 + 5*S^5*b^2;
	E3Div   = 7*b*S^4 - 6*R*S^3 - 69*b^2*S^2 + 54*b*R*S + 108*b^3;
	E3 = - E3Subst / E3Div;
	E2 = (E3 + b*S - R)*S / (2*b);
	#
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	x = as.vector(x);
	S = rep(S, each=3); E3 = rep(E3, each=3); R = rep(R, each=3);
	yz.s = S - x;
	y = (R - b*x)*x / E3;
	z = yz.s - y;
	sol = cbind(x=x, y=y, z=z, S=S);
	return(sol);
}
coeff.S3Ht.L21z = function(R, b, b.ext=0) {
	coeff = c(b^2, - 2*b*R, (R^2 - 3*b^3), 2*b^2*R, 4*b^4);
	# Extensions:
	if(any(b.ext != 0)) {
		coeff = coeff + c(b.ext[1]^2 + 2*b*b.ext[1], -2*R*b.ext[1], - 2*b^2*b.ext[1],0,0);
		# TODO: b.ext[2];
	}
	return(coeff);
}
### Test:
test.S3Ht.L21z = function(sol, b, b.ext=0, R=NULL) {
	if(length(b.ext) < 2) b.ext = c(b.ext, 0);
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	S = x+y+z; ext1 = b.ext[1]*S; ext2 = b.ext[2]*S^2;
	err1 = x^2*y + b*z + ext1 + ext2;
	err2 = y^2*z + b*x + ext1 + ext2;
	err3 = z^2*x + b*y + ext1 + ext2;
	err  = cbind(err1, err2, err3);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = -1
b.ext = 0
sol = solve.S3Ht.L21z(R, b, b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

test.S3Ht.L21z(sol, b, b.ext=b.ext)


### Ex 2:
R = 5
b = -2
sol = solve.S3Ht.L21z(R, b)

test.S3Ht.L21z(sol, b)


### Ex 3:
R = 3
b = -2
b.ext = -5
sol = solve.S3Ht.L21z(R, b, b.ext)

test.S3Ht.L21z(sol, b, b.ext=b.ext)


### Ex 4: R = 0
R = 0
b = 5
b.ext = -1
sol = solve.S3Ht.L21z(R, b, b.ext)

test.S3Ht.L21z(sol, b, b.ext=b.ext)


### Test
S = x+y+z; ext1 = b.ext[1]*S; ext2 = b.ext[2]*S^2;
x^2*y + b*z + ext1 + ext2 # - R
y^2*z + b*x + ext1 + ext2 # - R
z^2*x + b*y + ext1 + ext2 # - R

### Classic Polynomial:
round0.p(poly.calc(x))

### (x^3 + b*x - R) * P[12]
b^2*x^12 - 2*b*R*x^11 + (R^2 - b^3)*x^10 + 3*R*b^2*x^9 - b*(3*R^2 - b^3)*x^8 +
	+ R*(R^2 - 4*b^3)*x^7 + (2*R^2*b^2 - b^5)*x^6 + 5*b^4*R*x^5 +
	+ (R^4 - b^3*R^2 + b^6)*x^4 - b^2*R*(3*R^2 + 3*b^3)*x^3 +
	- b*(R^4 - 3*R^2*b^3 + b^6)*x^2 + (2*R^3*b^3 - R*b^6)*x +
	- R^2*b^5 + b^8 # = 0


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


########################
########################

########################
### Mixed-Order: 3+1 ###
########################

### x[i]^3*x[i+1] + b*x[i+2]

# x^3*y + b*z = R
# y^3*z + b*x = R
# z^3*x + b*y = R

### Solution:
# R1 = 3*R - b*S;
# - shortcut: using eq. from the HtSymmetric Mixed system;

### Sum + Rotation =>
E3*S^5 - 5*E2*E3*S^3 + (7*E3^2 - R1*E2)*S^2 + (E2^2*E3 + R1*E3)*S +
	+ R1^2 + 2*R1*E2^2 + E2^4 # = 0
# R1 = 3*R - b*S;
E3*S^5 - 5*E2*E3*S^3 + (7*E3^2 - (3*R - b*S)*E2)*S^2 + (E2^2*E3 + (3*R - b*S)*E3)*S +
	+ (3*R - b*S)^2 + 2*(3*R - b*S)*E2^2 + E2^4 # = 0
E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + b*E2*S^3 - 3*R*E2*S^2 + E2^2*E3*S + 3*R*E3*S - b*S*E3*S +
	- 2*b*S*E2^2 + E2^4 + 6*R*E2^2 + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0

### Sum(z*...) =>
x*y*z*(x^2 + y^2 + z^2) + b*(x^2 + y^2 + z^2) - R*S # = 0
E3*(S^2 - 2*E2) + b*S^2 - 2*b*E2 - R*S # = 0
E3*S^2 - 2*E2*E3 + b*S^2 - 2*b*E2 - R*S # = 0

### Sum(y*...) =>
(x^3*y^2 + y^3*z^2 + z^3*x^2) + b*E2 - R*S # = 0
### Rotation/Inversion =>
# - shortcut: using eq. from the HtSymmetric Mixed system;
# R31 = R*S - b*E2;
E3^2*S^4 + (2*R31*E3 + E2*E3^2)*S^2 - (R31*E2^2 + 5*E2^3*E3)*S +
	+ R31^2 + E2^5 + 7*E2^2*E3^2 + R31*E2*E3 # = 0
E3^2*S^4 + (2*(R*S - b*E2)*E3 + E2*E3^2)*S^2 - ((R*S - b*E2)*E2^2 + 5*E2^3*E3)*S +
	+ (R*S - b*E2)^2 + E2^5 + 7*E2^2*E3^2 + (R*S - b*E2)*E2*E3 # = 0
E3^2*S^4 + 2*R*E3*S^3 + E2*E3^2*S^2 - 2*b*E2*E3*S^2 - R*E2^2*S^2 +
	- 5*E2^3*E3*S + b*E2^3*S + R*E2*E3*S - 2*b*R*E2*S +
	+ E2^5 + 7*E2^2*E3^2 - b*E2^2*E3 + b^2*E2^2 +R^2*S^2 # = 0

### System:
E3*S^2 - 2*E2*E3 + b*S^2 - 2*b*E2 - R*S # = 0
E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + b*E2*S^3 - 3*R*E2*S^2 + E2^2*E3*S + 3*R*E3*S - b*S*E3*S +
	- 2*b*S*E2^2 + E2^4 + 6*R*E2^2 + b^2*S^2 - 6*b*R*S + 9*R^2 # = 0
E3^2*S^4 + 2*R*E3*S^3 + E2*E3^2*S^2 - 2*b*E2*E3*S^2 - R*E2^2*S^2 +
	- 5*E2^3*E3*S + b*E2^3*S + R*E2*E3*S - 2*b*R*E2*S +
	+ E2^5 + 7*E2^2*E3^2 - b*E2^2*E3 + b^2*E2^2 +R^2*S^2 # = 0

### Eq:
b1^3*S^11 - 3*R*b1^2*S^10 + 3*R^2*b1*S^9 - (R^3 + 9*b1^4)*S^8 + 20*R*b1^3*S^7 - 15*R^2*b1^2*S^6 +
	+ (2*R^3*b1 + 72*b1^5)*S^5 + (R^4 - 208*R*b1^4)*S^4 + 264*R^2*b1^3*S^3 +
	- (188*R^3*b1^2 + 64*b1^6)*S^2 + (80*R^4*b1 + 96*R*b1^5)*S - 16*R^5 - 48*R^2*b1^4

### Auxiliary Eqs:
### TODO: E3;
# E2 = (E3*S^2 + b*S^2 - R*S) / (2*(E3  + b));


### Test
x^3*y + b*z # - R
y^3*z + b*x # - R
z^3*x + b*y # - R


### P[33]
b1 = b[1];
(- R^3*b1^10 + b1^14) +
(- R^2*b1^11 + 3*R^5*b1^7)*x +
(- R*b1^12 + 3*R^4*b1^8 - 3*R^7*b1^4)*x^2 +
(4*R^3*b1^9 - 3*R^6*b1^5 + R^9*b1 - b1^13)*x^3 +
(- 7*R^2*b1^10 - 6*R^5*b1^6 + R^8*b1^2)*x^4 +
(2*R*b1^11 - 3*R^4*b1^7 + 4*R^7*b1^3)*x^5 +
(R^6*b1^4 - R^9 + b1^12)*x^6 +
(15*R^2*b1^9 - 2*R^5*b1^5)*x^7 +
(- 3*R*b1^10 + 17*R^4*b1^6 + R^7*b1^2)*x^8 +
(- 11*R^3*b1^7 + 5*R^6*b1^3 - b1^11)*x^9 +
(- 15*R^2*b1^8 + 6*R^5*b1^4 - R^8)*x^10 +
(4*R*b1^9 - 36*R^4*b1^5 - R^7*b1)*x^11 +
(21*R^3*b1^6 - 10*R^6*b1^2 + b1^10)*x^12 +
(10*R^2*b1^7 + 5*R^5*b1^3)*x^13 +
(- 5*R*b1^8 + 29*R^4*b1^4 - R^7)*x^14 +
(- 25*R^3*b1^5 + 7*R^6*b1 - b1^9)*x^15 +
(- 4*R^2*b1^6 - 12*R^5*b1^2)*x^16 +
(8*R*b1^7 - 7*R^4*b1^3)*x^17 +
(22*R^3*b1^4 - R^6 - b1^8)*x^18 +
(- 3*R^2*b1^5 + 6*R^5*b1)*x^19 +
(- 7*R*b1^6 - 6*R^4*b1^2)*x^20 +
(2*R^3*b1^3 + b1^7)*x^21 +
(- 6*R^2*b1^4 - R^5)*x^22 +
(6*R*b1^5 + 5*R^4*b1)*x^23 +
(- 10*R^3*b1^2 - b1^6)*x^24 +
(10*R^2*b1^3)*x^25 +
(- 5*R*b1^4 - R^4)*x^26 +
(4*R^3*b1 + b1^5)*x^27 +
(- 6*R^2*b1^2)*x^28 +
(4*R*b1^3)*x^29 +
(- R^3 - b1^4)*x^30 +
(3*R^2*b1)*x^31 +
(- 3*R*b1^2)*x^32 +
(b1^3)*x^33

### Debug
R = 2; b = -1;
x =  0.4910569800 - 0.8076735584i;
y = -1.9917091152 + 0.0333588950i;
z = -0.3198882225 + 0.0862258573i;
S = x+y+z; E3 = x*y*z; E2 = x*(y+z)+y*z;

coeff.f = function(R, b) {
	b1 = b[1];
	coeff = c(b1^3, - 3*R*b1^2, 3*R^2*b1, - R^3 - b1^4, 4*R*b1^3, - 6*R^2*b1^2, 4*R^3*b1 + b1^5,
		- 5*R*b1^4 - R^4, 10*R^2*b1^3, - 10*R^3*b1^2 - b1^6, 6*R*b1^5 + 5*R^4*b1, - 6*R^2*b1^4 - R^5,
		2*R^3*b1^3 + b1^7, - 7*R*b1^6 - 6*R^4*b1^2, - 3*R^2*b1^5 + 6*R^5*b1, 22*R^3*b1^4 - R^6 - b1^8,
		8*R*b1^7 - 7*R^4*b1^3, - 4*R^2*b1^6 - 12*R^5*b1^2, - 25*R^3*b1^5 + 7*R^6*b1 - b1^9,
		- 5*R*b1^8 + 29*R^4*b1^4 - R^7, 10*R^2*b1^7 + 5*R^5*b1^3, 21*R^3*b1^6 - 10*R^6*b1^2 + b1^10,
		4*R*b1^9 - 36*R^4*b1^5 - R^7*b1, - 15*R^2*b1^8 + 6*R^5*b1^4 - R^8,
		- 11*R^3*b1^7 + 5*R^6*b1^3 - b1^11, - 3*R*b1^10 + 17*R^4*b1^6 + R^7*b1^2,
		15*R^2*b1^9 - 2*R^5*b1^5, R^6*b1^4 - R^9 + b1^12, 2*R*b1^11 - 3*R^4*b1^7 + 4*R^7*b1^3,
		- 7*R^2*b1^10 - 6*R^5*b1^6 + R^8*b1^2, 4*R^3*b1^9 - 3*R^6*b1^5 + R^9*b1 - b1^13,
		- R*b1^12 + 3*R^4*b1^8 - 3*R^7*b1^4, - R^2*b1^11 + 3*R^5*b1^7, - R^3*b1^10 + b1^14);
	coeff;
}
solve.y = function(x, R, b) {
	b1 = b[1];
	y0 = - 2*R^2*b1^5*x^18 + 10*R^2*b1^6*x^15 + 2*R^3*b1^4*x^17 - 10*R^3*b1^5*x^14 +
	8*R^4*b1^2*x^19 - 4*R^4*b1^6*x^7 + 8*R^4*b1^7*x^4 + 2*R^4*b1^8*x - 8*R^5*b1*x^18 + 2*R^6*b1^3*x^8 - 8*R^6*b1^4*x^5 + 2*R^8*x^9 - 2*R^2*b1^11;
	y.div = - 2*R*b1^6*x^18 + 2*R*b1^12 + 2*R^2*b1^5*x^17 - 6*R^3*b1^2*x^22 - 2*R^3*b1^6*x^10 +
	10*R^3*b1^7*x^7 - 14*R^3*b1^8*x^4 + 6*R^4*b1*x^21 - 14*R^5*b1^3*x^11 + 20*R^5*b1^4*x^8 - 2*R^7*x^12

	y = - y0 / y.div;
}
solve.classic = function(R, b, debug=TRUE) {
	coeff = coeff.f(R, b);
	x = roots(coeff);
	isZero = (round0(x) == 0);
	if(debug) print(sum(isZero));
	x = x[ ! isZero]
	y = solve.y(x, R, b);
	z = (R - x^3*y) / b[1];
	return(cbind(x=x, y=y, z=z));
}
unique.S = function(sol, digits=5) {
	S = apply(sol, 1, sum);
	S = S[order(abs(S))];
	S = S[ ! duplicated(round(S, digits=digits))]
	return(S)
}

R = -2; b = 3;
sol = solve.classic(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];
S = unique.S(sol); # S = unique.S(sol[-c(1,2,3),])
round0.p(poly.calc(S) * b^3)


#######################
#######################

#####################
### Bi-Mixt-Order ###
#####################

### x[i]^2*x[j] + a*x[i]*x[j]^2

# x^2*y + a*x*y^2 = R
# y^2*z + a*y*z^2 = R
# z^2*x + a*z*x^2 = R

### Solution:

### Sum =>
(x^2*y + y^2*z + z^2*x) + a*(x*y^2 + y*z^2 + z*x^2) - 3*R # = 0
(x^2*y + y^2*z + z^2*x + x*y^2 + y*z^2 + z*x^2) +
	+ (a-1)*(x*y^2 + y*z^2 + z*x^2) - 3*R # = 0
(E2*S - 3*E3 - 3*R) +
	+ (a-1)*(x*y^2 + y*z^2 + z*x^2) # = 0 # Eq 1-bis
### * (x^2*y + y^2*z + z^2*x) =>
(E2*S - 3*E3 - 3*R) * (x^2*y + y^2*z + z^2*x) +
	+ (a-1)*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) # = 0 # Eq 2-bis
### (...)*Eq 1-bis + (a-1)*Eq 2-bis =>
(a-1)*(E2*S - 3*E3 - 3*R) * (x^2*y + y^2*z + z^2*x + x*y^2 + y*z^2 + z*x^2) +
	+ (a-1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (E2*S - 3*E3 - 3*R)^2 # = 0
(a-1)*(E2*S - 3*E3 - 3*R) * (E2*S - 3*E3) +
	+ (a-1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (E2*S - 3*E3 - 3*R)^2 # = 0
(a-1)^2*E2^3 + 9*(a+1)*E3*R + 9*E3^2 + 9*a*(a-1)*E3^2 +
	- 6*E2*E3*S - 6*a*(a-1)*E2*E3*S - 3*(a+1)*R*E2*S +
	+ E2^2*S^2 + (a-1)*E2^2*S^2 + (a-1)^2*E3*S^3 + 9*R^2

### Sum(z*...) =>
(a+1)*E3*S - R*S # = 0
# S = 0, OR
# (a+1)*E3 = R;

### Sum(z^2*...) =>
(a+1)*E3*E2 - R*(x^2 + y^2 + z^2) # = 0
(a+1)*E3*E2 - R*(S^2 - 2*E2) # = 0
(a+1)*E3*E2 - R*S^2 + 2*R*E2 # = 0

### Alternative: Diff;

### Case 1: S = 0
(a+1)*E3*E2 + 2*R*E2 # = 0
E2*((a+1)*E3 + 2*R) # = 0
# (a+1)*E3 = - 2*R, OR
# E2 = 0;

### Case 1.1: E2 != 0:
(a-1)^2*E2^3 + 9*(a+1)*E3*R + 9*E3^2 + 9*a*(a-1)*E3^2 + 9*R^2 # = 0
# =>
(a+1)^2*E2^3 + 27*R^2 # = 0

### Case 1.2: E2 = 0
E3^2 + a*(a-1)*E3^2 + (a+1)*E3*R + R^2 # = 0

### Case 2:
# TODO


### Solver:
solve.LeadA.S3P21 = function(R, a, b.ext=0, debug=TRUE) {
	ap = a[1] + 1;
	an = a[1] - 1;
	### Case 1: S = 0
	S = 0;
	### Case 1.1: E2 != 0
	E3 = - 2*R / ap;
	# currently NO correct solutions: TODO: examine & update;
	E2 = roots(c(ap^2, 0, 0, 27*R[1]^2));
	if(debug) print(E2);
	len = length(E2);
	# S = rep(S, len); E3 = rep(E3, len);
	
	### Case 1.2: E2 = 0
	E3 = roots(c(1 + a*an, ap*R[1], R[1]^2))
	if(debug) print(E3);
	len = length(E3);
	S = rep(S, len); E2 = 0; E2 = rep(E2, len);
	### Case 2: TODO !!!
	#
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	#
	rep.m = function(x) matrix(x, ncol=len, nrow=3, byrow=TRUE)
	S = rep.m(S); E3 = rep.m(E3);
	R1 = R[1];
	yz.s = S - x;
	yz = E2 - x*yz.s;
	yz.sa = R1 / yz;
	y = (a*yz.s - yz.sa) / an;
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol)
}

### Examples:
R = 2
a = -2
b.ext = c(0, 0)
sol = solve.LeadA.S3P21(R, a)
x = sol[,1]; y = sol[,2]; z = sol[,3];


### Test
S = x+y+z; ext1 = b.ext[1]*S; ext2 = b.ext[2]*S^2;
x^2*y + a*x*y^2 + ext1 + ext2 # - R
y^2*z + a*y*z^2 + ext1 + ext2 # - R
z^2*x + a*z*x^2 + ext1 + ext2 # - R

### Classic Polynomial:
round0.p(poly.calc(x)) * 7
4 - 2*x^3 + 7*x^6


### Debug:
x =  0.8161669286 - 0.4045961913i;
y = -0.0576928844 + 0.9091193896i;
z = -0.7584740442 - 0.5045231983i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;



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
