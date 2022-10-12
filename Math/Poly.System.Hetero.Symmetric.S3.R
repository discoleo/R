########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogeneous Symmetric
###
### draft v.0.4g-clean5


### Hetero-Symmetric
### Polynomial Systems: 3 Variables

### Example:
x^n + P(x, y, z) = R
y^n + P(y, z, x) = R
z^n + P(z, x, y) = R

####################

###############
### History ###
###############


### draft v.0.4g:
# - robust solutions for S3P3 Simple & S3P2-Asymmetric Sum;
### draft v.0.4f-clean 1 & 2 & 3:
# - more cleaning;
# - S3P3 MixedSideChain: factorized P[12] => P[8]
#   => 3*8 = 24 true solutions; [v.0.4f-TrueSol]
### draft v.0.4e:
# - cleanup: moved to new file
#   Poly.System.Hetero.Symmetric.S3.Derivation.R;
# - work on x^3 + b2*y^2 + b1*z = R;
### draft v.0.4c - v.0.4d:
# - solved Order 3:
#   x^3 + b*y*z = R;
# - TODO: Case x == y, but != z;
# - moved to new file:
#   Poly.System.Hetero.Symmetric.S3.YZ.R; [draft v.0.4d]
### draft v.0.4b - v.0.4b-full:
# - [full] solution for:
#   x^3 + y^2 = R; [partial solution]
#   x^3 + b2*y^2 + b1*y = R; [v.0.4b-full]
### draft v.0.4a:
# - moved to new file:
#   systems with Composite Leading Term: e.g. x^m*y^n;
# - file: Poly.System.Hetero.Symmetric.S3.Leading.R;
### draft v.0.3h - v.0.3h-fix:
# - Order 3: (partially solved + partial bug fix)
#   x^3 + a1*y^3 + a2*z^3 = R;
### draft v.0.3g - v.0.3g-fix:
# - solved (with extensions of type A1):
#   x^2 + a1*y^2 + a2*z^2 = R;
# - more exploration of this system; [v.0.3g-expl]
# - fixed / full equation; [v.0.3g-fix]
### draft v.0.3f:
# - Structural Extension:
#   a2*(x*y*z)^2 + a1*x*y*z + x*y + b1*y = R;
### draft v.0.3e:
# - initial work on: x^3 + b2*y^2 + b1*y = R;
### draft v.0.3d - v.0.3d-simple:
# - solved: x*y + b*y = R;
# - Note: only A1-type extensions have distinct solutions;
# - simplification of the correct solution (for A1-extensions);
### draft v.0.3c - v.0.3c-poly-shift:
# - solved: x^2 + y^2 + b1*y = R;
# - added Extensions of type A1; [v.0.3c-ext]
# - added Classical polynomial: P[6]; [v.0.3c-poly & fixed minor bug]
# - special Case P[6]: b.ext[2] = -3;
#   e.g. 11 + 2*x + 5*x^2 - 2*x^3 + x^6 = 0;
# - special Case with shifted P[6]: b5 = 0;
### draft v.0.3b - v.0.3b-P2ext:
# - solved: x^3 + b*y*z = R;
# - reordering of sections & better comments; [v.0.3b-ord]
# - fix of the wrong roots (in an older Shift-x system); [v.0.3b-fix]
# - extension A1 for the P2 system;
# - TODO: cleanup;
### draft v.0.3a-ext:
# - extensions of type A1 for Ht S3P2;
# - simplification of the base Eq for Ht S3P2;
### draft v.0.3a-pre:
# - moved Difference types to new file:
#   Poly.System.Hetero.Symmetric.S3.Diff.R;
### draft v.0.2b - v.0.2b-ext: [MOVED]
# - first concepts / solved [v.0.2b-sol]:
#   x^2 - y^2 + b*x*y = R;
# - solved extension A1 (Pow 1): + b[2]*(x+y+z); [v.0.2b-ext]
### draft v.0.2a - v.0.2a-poly:
# - solved: Ht[3, 3, 1];
# - simplified solution: from P11 to P8; [v.0.2a-simplify]
# - classic polynomial: P24; [v.0.2a-poly]
# - TODO: find robust solution;
### draft v.0.1d - v.0.1d-poly:
# - minor fix: in Ht[3, 2, 1];
# - classic polynomial (P6) for Ht[3, 2, 1]; [v.0.1d-poly]
### draft v.0.1c-move:
# - moved Hetero-Mixt (v.0.1c) to separate file;
### draft v.0.1c-pre-alpha - v.0.1c-exact:
# - first look & solved + exact/robust solution: [v.0.1c-exact]
#   x*y^2 + y*z^2 + z*x^2 = R1;
#   [moved to file: Poly.System.Hetero.Symmetric.S3.Mixt.R]
### draft v.0.1b - v.0.1b-fix:
# - solved: x[i]^2 + b2*x[j] + b1*x[k];
# - classical Polynomial (P8) (v.0.1b-clP; fixed in v.0.1b-fix);
#  (done: P8 = P2*P6 in v.0.1b-fix)
# - extension: x[i]^2 + s*x + b3*x*y*z + b2*x[j] + b1*x[k]; (v.0.1b-ext)
# - TODO: find/correct (precision) bug vs correct roots;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Hetero.Symmetric.R;
########################
### former branch v.0.2:
### in Poly.System.Hetero.Symmetric.R
### draft v.0.2.d:
# - more work on x[i]^3 + b*x[i+1] = R;
# - added: x^2*y*z + b*y = R;
# - various formatting improvements;
### draft v.0.2.c:
# - solved: x[i]^2 + s*x[i] + b*x[i+1] = R;
# - TODO: correct various bugs [DONE];
### draft v.0.2b:
# - some exploration of systems with x*y*z terms;
### branch v.0.2a:
# - more work on systems with 3 variables:
#   "proper" implementation of: x[i]^2 + b*x[k];
# - TODO: robust removal of set of wrong solutions;
#   [or avoid getting superfluous solutions ???]
### branch v.0.2a-pre-a:
# - initial work on systems with 3 variables;
# - the simple cases are less rewarding;


####################

####################
### Introduction ###
####################

### 3 Variables:
# x^n + P(x, y, z) = R
# y^n + P(y, z, x) = R
# z^n + P(z, x, y) = R

### Very Simple System:
# x^n + b*y = R
# y^n + b*z = R
# z^n + b*x = R

### Solution:
# - decomposing system into a lower order system;
# - the trivial solution x = y = z will not be covered here;

### Sum =>
# (x^n + y^n + z^n) + b*S = 3*R

### Sum(x[i] * P[i]) =>
# (x^(n+1) + y^(n+1) + z^(n+1)) + b*E2 = R*S

###  Eq 3:
# - can be based on:
E21a = x^2*y + y^2*z + z^2*x;
E12a = x*y^2 + y*z^2 + z*x^2;
#
E21a + E12a - E2*S + 3*E3 # = 0
E21a * E12a - E2^3 + 6*E2*E3*S - E3*S^3 - 9*E3^2 # = 0
# =>
E21a^2 - (E2*S - 3*E3)*E21a + E2^3 - 6*E2*E3*S + E3*S^3 + 9*E3^2 # = 0

### Eq 3: [alternative]
# - can be based on: E[n,1]a;
En1a + b*(S^2 - 2*E2) - R*S # = 0
E1na + b*E2 - R*S # = 0
# =>
Sn*S - S(n+1) + b*S^2 - b*E2 - 2*R*S # = 0

### Eq 3: [alternative]
# x^n = R - b*y
# x^(2*n) = R^2 + b^2*y^2 - 2*b*R*y
### [3a] Sum =>
# x^(2*n) + y^(2*n) + z^(2*n) = 3*R^2 + b^2*S^2 - 2*b^2*E2 - 2*b*R*S;


# Note:
# Sn = (x^n + y^n + z^n) as well as the (n+1) variant:
# - can be decomposed as polynomials of: S, E2, E3;,
#   where S = x + y + z;
# - S, E2, E3 = elementary polynomials; 
# n:   (x^n + y^n + z^n) = Decomp1(S^n, E2, E3);
# n+1: (x^m + y^m + z^m) = Decomp2(S^(n+1), E2, E3), where m = n+1;

# Alternative to Eq3: for very low orders;
### Diff =>
# x^n - y^n = b*(z - y)
### Prod =>
# (x^n - y^n)*(y^n - z^n)*(z^n - x^n) = b^3*(z-y)*(x-z)*(y-x)

### Eq 3-variant: [but NO benefit]
# x^n - b*x = R - b*y - b*x
# x^(2*n) - 2*b*x^(n+1) + b^2*x^2 = R^2 + b^2*y^2 + b^2*x^2 - 2*b*R*y - 2*b*R*x + 2*b^2*x*y
### [3b] Sum =>
# x^(2*n) + y^(2*n) + z^(2*n) - 2*b*(x^(n+1) + y^(n+1) + z^(n+1)) + b^2*S^2 - 2*b^2*E2 =
#  = 3*R^2 + 2*b^2*(S^2 - 2*E2) - 4*b*R*S + 2*b^2*E2
### Sum[3b] - Sum[3a] + 2*b*Sum[2] =>
# b^2*S^2 - 2*b^2*E2 + 2*b^2*E2 = 2*b*R*S + 3*R^2 + 2*b^2*(S^2 - 2*E2) - 4*b*R*S + 2*b^2*E2 - (3*R^2 + b^2*S^2 - 2*b^2*E2 - 2*b*R*S)
# 0 == 0; # [equations are correlated]


### Complexity:
# - initial System: => order P[n^3];
#  -- polynomial can be decomposed = P[n]*P[n^3 - n];
# - Decomposed System:
#  -- P[3] o P[(n^3 - n)/3];
#  -- D(S, E2, E3): orders of E2 & E3 are usually much lower;


####################
####################

### Helper Functions

# library(polynom)
# library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

source("Polynomials.Helper.R")

### Other

test.S3Ht.Simple = function(sol, b, R=NULL, n) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	# Extensions:
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: power 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A1: power 2;
	S  = (x+y+z); ext1 = b2*S; ext2 = b3*S^2;
	ext  = ext1 + ext2;
	err1 = x^n + b[1]*y + ext;
	err2 = y^n + b[1]*z + ext;
	err3 = z^n + b[1]*x + ext;
	err = rbind(err1, err2, err3);
	err = round0(err);
	return(err);
}
### Asymmetric Sum:
test.S3Ht.SumYZ = function(sol, b, b.ext=0, R=NULL, n) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	# Extensions:
	b2 = b.ext[1];
	b3 = if(length(b.ext) > 1) b.ext[2] else 0; # Ext A1: power 2;
	S  = (x+y+z); ext1 = b2*S; ext2 = b3*S^2;
	ext  = ext1 + ext2;
	err1 = x^n + b[1]*y + b[2]*z + ext;
	err2 = y^n + b[1]*z + b[2]*x + ext;
	err3 = z^n + b[1]*x + b[2]*y + ext;
	err = rbind(err1, err2, err3);
	err = round0(err);
	return(err);
}
test.S3Ht.Product = function(sol, b, b.ext = 0, R=NULL, n) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	# Extensions:
	b2 = b.ext[1];
	b3 = if(length(b.ext) >= 2) b.ext[2] else 0; # Ext A1: power 2;
	S  = (x+y+z); ext1 = b2*S; ext2 = b3*S^2;
	ext  = ext1 + ext2;
	err1 = x^n + b[1]*y*z + ext;
	err2 = y^n + b[1]*z*x + ext;
	err3 = z^n + b[1]*x*y + ext;
	err = rbind(err1, err2, err3);
	err = round0(err);
	return(err);
}
test.S3Ht.Asym3ExtE3 = function(sol, b, s, b3, R=NULL, n) {
	x = sol[,1]; y = sol[,2]; z = sol[,3];
	b1 = b[1]; b2 = b[2];
	# Quasi-Extension:
	ext  = b3 * x*y*z;
	err1 = x^n + s*x + b2*y + b1*z + ext;
	err2 = y^n + s*y + b2*z + b1*x + ext;
	err3 = z^n + s*z + b2*x + b1*y + ext;
	err = rbind(err1, err2, err3);
	err = round0(err);
	return(err);
}


################################
################################

##############
### Simple ###
##############

### x^n + b*y = R

###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1]

# x^2 + b1*y = R
# y^2 + b1*z = R
# z^2 + b1*x = R

### Extensions:

### Ext A1: power 1
# x^2 + b1*y + b2*(x + y + z) = R
### Ext A1: power 2
# x^2 + b1*y + b2*(x + y + z) + b3*(x + y + z)^2 = R


### Solution

# Trivial solution: x = y = z;

### Case: (x, y, z) distinct

### Method 1:
### Eq 1: Sum =>
S^2 - 2*E2 + b1*S - 3*R # = 0

### Eq 2: Sum(x[i]*...) =>
6*E3 - S^3 - 2*b1*S^2 + 7*R*S + b1^2*S - 3*b1*R

### Eq 3: Sum(x[i+1]^2*...) =>
S^4 + 2*b1*S^3 - (10*R + b1^2)*S^2 + 6*(b1*R + b1^3)*S - 18*b1^2*R + 9*R^2

### Eq 3: [alternative]
# - based on E[n,1]a;
E2*S - 3*E3 + b1*S^2 - b1*E2 - 2*R*S # = 0

### Eq:
(S^2 + 3*b1*S - 9*R)*(S^2 - b1*S - R + 2*b1^2)

### Auxiliary Eqs:
# E2 = (S^2 + b[1]*S - 3*R1) / 2;
# E3 = - (S^3 - 3*E2*S + b[1]*E2 - R1*S) / 3;

### Solver:
solve.S3Ht.P2 = function(R, b, debug=TRUE) {
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: power 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A1: power 2;
	# coeff = c(1, 2*b[1], - (10*R[1] + b[1]^2), 6*(b[1]*R[1] + b[1]^3), - 18*b[1]^2*R[1] + 9*R[1]^2)
	# if(b2 != 0) coeff = coeff + c(0, 10*b2, -6*b[1]*b2 + 9*b2^2, 18*b[1]^2*b2 - 18*R[1]*b2, 0)
	coeff = c(1 + b3, b2 - b[1], - R[1] + 2*b[1]^2)
	S = roots(coeff)
	if(debug) print(S);
	len = length(S);
	if(len == 0) stop("NO solutions!")
	R1 = R[1] - b2*S - b3*S^2;
	# [REMOVED] remove x == y == z = S / 3;
	# as it causes numerical instability due to multiplicity of roots;
	# isEq = round0(S^2 + 3*b[1]*S - 9*R1) == 0
	# if(any(isZero)) print("Warning: f(S) == 0!")
	# S = S[ ! isEq]; R1 = R1[ ! isEq];
	E2 = round0(S^2 + b[1]*S - 3*R1) / 2;
	E3 = - (S^3 - 3*E2*S + b[1]*E2 - R1*S) / 3;
	E3 = round0(E3, tol=1E-10); # improve numerics when E3 == 0;
	#
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	x = as.vector(x);
	R1 = rep(R1, each=3);
	y = (R1 - x^2)/b[1];
	z = (R1 - y^2)/b[1];
	sol = cbind(x, y, z);
	return(sol);
}
test.S3Ht.P2 = function(sol, b, R=NULL) {
	err = test.S3Ht.Simple(sol, b=b, R=R, n=2);
	return(err);
}

### Examples:
R = 3
b = -1
#
sol = solve.S3Ht.P2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test.S3Ht.P2(sol, b=b)

### Classic Polynomial: P8 or P6 (when S == 0)
round0.p(poly.calc(sol[,1]))


### Classic Polynomial:
- R*b[1]^4 + 2*R^2*b[1]^2 - R^3 + b[1]^6 + (2*R*b[1]^3 - R^2*b[1] - b[1]^5)*x + (3*R^2 - 3*R*b[1]^2 + b[1]^4)*x^2 +
	+ (2*R*b[1] - b[1]^3)*x^3 - (3*R - b[1]^2)*x^4 - b[1]*x^5 + x^6

### Classic Polynomial: Shifted
# b10 = b1 / 6; x0 = x - b1;
(- 963*R*b1^4 + 69*R^2*b1^2 - R^3 + 39991*b1^6) +
	(240*R*b1^3 - 5712*b1^5)*x + (- 90*R*b1^2 + 3*R^2 + 819*b1^4)*x^2 +
	(- 112*b1^3)*x^3 + (- 3*R + 21*b1^2)*x^4 + x^6


#########
### Ex 2:
R = 3
b = 3
#
sol = solve.S3Ht.P2(R, b=b)

### Test
test.S3Ht.P2(sol, b=b)

round0.p(poly.calc(sol[,1]))

621 - 108*x + 27*x^2 - 9*x^3 - 3*x^5 + x^6


###############
### Extensions:
R = 1
b = c(1, 1)
#
sol = solve.S3Ht.P2(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
test.S3Ht.P2(sol, b=b)

round0.p(poly.calc(sol[,1]))


### Ex 2:
R = -1
b = c(1, 1)
#
sol = solve.S3Ht.P2(R, b=b)

### Test
test.S3Ht.P2(sol, b=b)

x = sol[,1]
round0.p(poly.calc(sol[,1]))
err = 25 + 12*x^2 - 2*x^3 + 3*x^4 + x^6
round0(err)


### Ex 3:
R = 2
b = c(2, 0, 1)
#
sol = solve.S3Ht.P2(R, b=b)

### Test
test.S3Ht.P2(sol, b=b)

x = sol[,1]
round0.p(poly.calc(sol[,1]))
err = 115 + 39*x + 44*x^2 - x^3 - 12*x^4 - x^5 + x^6
round0(err)


##################
##################

### Shifted roots:
### x[i]^2 + s*x[i] + b*x[i+1]

# x^2 + s*x + b1*y = R
# y^2 + s*y + b1*z = R
# z^2 + s*z + b1*x = R

### Solution

### Fast Solution:
# - shift all roots "back" & use solution to simple system: x^2 + b1*y = Rs,
#   where Rs = R + s^2/4 + b1*s/2;
(x + s/2)^2 + b1*(y + s/2) = R + s^2/4 + b1*s/2
(y + s/2)^2 + b1*(z + s/2) = R + s^2/4 + b1*s/2
(z + s/2)^2 + b1*(x + s/2) = R + s^2/4 + b1*s/2

### Long Solution:
# - moved to file:
#   Poly.System.Hetero.Symmetric.S3.Derivation.R;

### Eq:
(2*S + 3*s + 3*b1) * (S^2 - (b1 - 3*s)*S + (- R - 2*b1*s + 2*b1^2 + 2*s^2))


### Solver:
solve.S3Ht.P2Shift = function(R, b, s, debug=TRUE, all=FALSE) {
	coeff = c(1, - (b[1] - 3*s), (- R[1] - 2*b[1]*s + 2*b[1]^2 + 2*s^2))
	x.sum = roots(coeff);
	if(debug) print(x.sum);
	E3 = (x.sum^3 + (3*s + 2*b[1])*x.sum^2 +
		(2*s-b[1])*(s+b[1])*x.sum - 7*R*x.sum - 6*s*R + 3*b[1]*R) / 6;
	E2 = (x.sum^2 + (s+b[1])*x.sum - 3*R)/2
	x = as.vector(sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id]))))
	y = (R - x^2 - s*x)/b[1]
	z = (R - y^2 - s*y)/b[1]
	sol = cbind(x, y, z)
	if(all) {
		sol = rbind(sol, sol[,c(2,3,1)], sol[,c(3,1,2)]);
	}
	return(sol);
}

### Example
R = 1
b = 3
s = 1
#
sol = solve.S3Ht.P2Shift(R, b, s)
x = sol[,1]; y = sol[,2]; z = sol[,3]
sol

### Test
x^2 + s*x + b[1]*y
y^2 + s*y + b[1]*z
z^2 + s*z + b[1]*x

### Classic polynomial
round0.p(poly.calc(x[1:6]))


########################
########################

### Variant:
### x[i]^2 + b*(x[j] + x[k])

# x^2 + b1*(y+z) = R
# y^2 + b1*(x+z) = R
# z^2 + b1*(x+y) = R
# [a trivial system]


### Solution:
# Trivial solution: x = y = z;

### Diff =>
# x^2 - y^2 = b1*(x-y)
# x^2 - z^2 = b1*(x-z)
# y^2 - z^2 = b1*(y-z)

# if x != y != z
# x + y = b1
# x + z = b1
# y + z = b1
# => x = y = z, which violates assumption;


# Case: x = y
# x^2 + b1*(x+z) = R
# z^2 + 2*b1*x = R

# Case x != z
# x + z = b1 # => (Eq 1) =>
# x^2 + b1^2 = R

### Solver:

solve.S3Ht.P2ChsYZ = function(R, b, debug=TRUE, all=FALSE) {
	x = sqrt(R - b[1]^2 + 0i);
	x = c(x, -x);
	if(debug) print(x);
	y = x;
	z = b[1] - x;
	sol = cbind(x, y, z);
	if(all) sol = rbind(sol, sol[, c(2,3,1)], sol[, c(3,1,2)]);
	return(sol);
}

### Example:
R = 1
b = 3
#
sol = solve.S3Ht.P2ChsYZ(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + b[1]*(y+z)
y^2 + b[1]*(x+z)
z^2 + b[1]*(x+y)


####################

### Shifted Variant:
### (x[i] - s)^2 + b * sum(x[-i])

# (x-s)^2 + b1*(y+z) = R
# (y-s)^2 + b1*(x+z) = R
# (z-s)^2 + b1*(x+y) = R


### Solution

### Diff =>
# (x-s)^2 - (y-s)^2 = b1*(x - y)
# (x - y)*(x + y - 2*s) = b1*(x - y)
# (x - y)*(x + y - 2*s - b1) = 0
# (x - z)*(x + z - 2*s - b1) = 0
# (y - z)*(y + z - 2*s - b1) = 0

### Case x != y != z
# is NOT solvable (except in very special conditions);

### Case x = y && x != z
# x + z - 2*s - b1 = 0
# x + z = b1 + 2*s;
# =>
# (x-s)^2 + b1*(x+z) - R = 0
# x^2 - 2*s*x + s^2 + b1*(b1 + 2*s) - R
# x^2 - 2*s*x + b1^2 + 2*b1*s + s^2 - R

### Example
R = 1
b = 3
s = -1
#
x = roots(c(1, - 2*s, b[1]^2 + 2*b[1]*s + s^2 - R))
y = x
z = b[1] + 2*s - x
sol = cbind(x, y, z)
sol


### Test
(x-s)^2 + b[1]*(y+z)
(y-s)^2 + b[1]*(x+z)
(z-s)^2 + b[1]*(x+y)


####################
####################

####################
### Product-type ###
####################

### x[i]^2 + b * prod(x[-i])

# x^2 + b1*y*z = R
# y^2 + b1*x*z = R
# z^2 + b1*x*y = R

### Solution

# - degenerate P4;

### Special Case: b1 = 2
# S^2 = 3*R

### Diff =>
# x^2 - y^2 = b1*z*(x-y)
# x^2 - z^2 = b1*y*(x-z)
# y^2 - z^2 = b1*x*(y-z)
# => if x != y != z
  x + y - b1*z # = 0
  x - b1*y + z # = 0
- b1*x + y + z # = 0
# => x = y = z = 0; # Contradiction !!!
# - applies in both cases: (b1 != 2) & (b1 == 2);

### Case: x == y & x != z =>
# x + z = b1*x; # =>
# z = (b1 - 1)*x; # =>
(b1^2 - b1 + 1)*x^2 - R # = 0

### Case: x = y
# x^2 + b1*x*z = R
# z^2 + b1*x^2 = R
# =>
# b1*z = R/x - x
# b1^2*z^2 + b1^3*x^2 = b1^2*R
x^2 - 2*R + R^2/x^2 + b1^3*x^2 - b1^2*R # = 0
### Eq:
(b1^3+1)*x^4 - R*(b1^2 + 2)*x^2 + R^2 # = 0
((b1^2 - b1 + 1)*x^2 - R) * ((b1+1)*x^2 - R)


### Solver:

solve.S3Ht.P2ChpYZ = function(R, b, b.ext=0, debug=TRUE, all=TRUE) {
	coeff1 = c((b[1]^2 - b[1] + 1), 0, -R[1]); # x == y;
	coeff2 = c((b[1] + 1), 0, -R[1]); # x == y == z;
	# computes x, NOT S;
	# simple extensions: adapted to x;
	if(b.ext[1] != 0) {
		coeff1 = coeff1 + c(0, b.ext[1] * (b[1]+1), 0); # x == y
		coeff2 = coeff2 + c(0, 3*b.ext[1], 0);
	}
	sol1 = roots(coeff1);
	if(debug) print(sol1);
	if(all) {
		sol2 = roots(coeff2);
		x = c(sol1, sol2);
	}
	y = x;
	# z = (b[1] - 1)*x; # assumes x == y & x != z;
	z = (R[1] - x^2 - 2*x*b.ext[1]) / (b[1]*x + b.ext[1]);
	sol = round0(cbind(x=x, y=y, z=z))
}
test.S3Ht.P2ChpYZ = function(sol, b, b.ext = 0, R=NULL) {
	err = test.S3Ht.Product(sol, b=b, b.ext=b.ext, R=R, n=2);
	return(err)
}

### Example:
R = 2
b = 1
#
sol = solve.S3Ht.P2ChpYZ(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
test.S3Ht.P2ChpYZ(sol, b=b)


### Extension A1:
R = 2
b = 1
b.ext = c(1)
#
sol = solve.S3Ht.P2ChpYZ(R, b, b.ext)

### Test
test.S3Ht.P2ChpYZ(sol, b=b, b.ext=b.ext)


######################
######################

### Prod-Type: Order 3
### x[i]^3 + b * prod(x[-i])

# x^3 + b1*y*z = R
# y^3 + b1*x*z = R
# z^3 + b1*x*y = R

### Solution

### Case 1: x = y = z
# x^3 + b1*x^2 - R = 0
### Case 2: x != y != z
# x^3 + b1*x^2 + b1^2*x + b1^3 - R = 0
### Case 3: x = y != z
# - remaining cases;
# - unfortunately I haven't found a simpler solution yet!

### Sum =>
x^3 + y^3 + z^3 + b1*E2 - 3*R # = 0
S^3 - 3*E2*S + 3*E3 + b1*E2 - 3*R

### Sum(x*...) =>
x^4 + y^4 + z^4 + 3*b1*E3 - R*S # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 + 3*b1*E3 - R*S

### Sum(y*z*...) =>
E3*(x^2 + y^2 + z^2) + b1*((x*y)^2 + (x*z)^2 + (y*z)^2) - R*E2 # = 0
E3*(S^2 - 2*E2) + b1*(E2^2 - 2*E3*S) - R*E2 # = 0
E3*S^2 - 2*E2*E3 + b1*E2^2 - 2*b1*E3*S - R*E2 # = 0

###
E3 = (54*R + 31*S*b1^2 + 12*S^2*b1 - 30*S^3 - 15*b1^3) /
     (27*R*S*b1 - 72*R*S^2 + 45*R*b1^2 - 15*S^3*b1^2 + S^4*b1 + 12*S^5)
E2 = (S^3 + 3*E3 - 3*R) / (3*S - b1)


### Eq: (S^3 + 3*b1*S^2 - 27*R) * (6*S - 5*b1) * (S + b1) * P[6]
(5*R*b1^3 + 27*R^2) +
(18*R*b1^2)*S^1 +
(- 27*R*b1)*S^2 +
(- 5*b1^3)*S^3 +
(7*b1^2)*S^4 +
(- 3*b1)*S^5 +
(1)*S^6


### Solver:

solve.S3Ht.P3ChpYZ = function(R, b, debug=TRUE, all=FALSE, tol=1E-8) {
	b1 = b[1]; R = R[1];
	### Case: 2 equal
	coeff = c(1, - 3*b1, 7*b1^2, -5*b1^3, - 27*R*b1, 18*R*b1^2, (5*R*b1^3 + 27*R^2))
	S = roots(coeff)
	# Case: x != y != z
	# S = c(S, -b1);
	if(debug) print(S);
	#
	E2.x0  = 27*R*S*b1 - 72*R*S^2 + 45*R*b1^2 - 15*S^3*b1^2 + S^4*b1 + 12*S^5;
	E2.div = (54*R + 31*S*b1^2 + 12*S^2*b1 - 30*S^3 - 15*b1^3);
	E2 = - E2.x0 / E2.div;
	E3 = - (S^3 - 3*E2*S + b1*E2 - 3*R) / 3;
	# Robust complex solutions for x:
	# - but still some numerical instability!
	E2 = round0(E2, tol=tol); E3 = round0(E3, tol=tol);
	# Robust: P[2] instead of P[3];
	# x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	len = length(S);
	x = sapply(seq(len), function(id) {
		S = S[id]; E2 = E2[id]; E3 = E3[id];
		roots(c(S^2 - E2, E3 - E2*S - R, E3*b1 + E3*S));
	})
	len = if(is.matrix(x)) nrow(x) else 1;
	x  = as.vector(x);
	solve.x = function(len) {
		S  = rep(S, each=len);
		E3 = rep(E3, each=len);
		# yz = E3 / x;
		yz = (R - x^3)/b1; # robust ???
		yz.s = S - x;
		#
		yz.d = sqrt(yz.s^2 - 4*yz)
		y = (yz.s + yz.d)/2;
		z = (yz.s - yz.d)/2;
		sol = cbind(x, y, z);
	}
	sol = solve.x(len);
	if(all) {
		# only one permutation is missing;
		sol = rbind(sol, sol[, c(2,3,1)]);
	}
	### Case: all distinct
	S = - b1; E2 = b1^2; E3 = (R - b1^3);
	x = roots(c(1, -S, E2, -E3));
	sol2 = solve.x(len=3);
	# All permutations: x != y != z;
	# - actually only (z, y) is missing (which is actually valid);
	if(all) sol2 = rbind(sol2, sol2[, c(1,3,2)]);
	sol = rbind(sol, sol2);
	return(sol)
}
test.S3Ht.P3ChpYZ = function(sol, b, b.ext = 0, R=NULL) {
	err = test.S3Ht.Product(sol, b=b, b.ext=b.ext, R=R, n=3);
	return(err)
}

### Examples:

# - still some numerical instability!

R = 1;
b = 2;
#
sol = solve.S3Ht.P3ChpYZ(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

test.S3Ht.P3ChpYZ(sol, b)


### Ex 2:
R = -5;
b = 3;
#
sol = solve.S3Ht.P3ChpYZ(R, b)

test.S3Ht.P3ChpYZ(sol, b)


### Test
x^3 + b*y*z # - R
y^3 + b*x*z # - R
z^3 + b*x*y # - R


### Debug
R = 2; b = 3;
x = -0.5886981678 + 1.2138678310i
y = -0.2089400966 - 0.4200165444i
z = -0.2089400966 - 0.4200165444i
S = (x+y+z); E2 = x*(y+z) + y*z; E3 = x*y*z;


#########################
#########################

#########################
### "Asymmetric" Variant:
### x[i]^2 + b1*x[j] + b2*x[k]

# x^2 + b1*y + b2*z = R
# y^2 + b1*z + b2*x = R
# z^2 + b1*x + b2*y = R

### Solution:

### Sum =>
S^2 - 2*E2 + (b1+b2)*S - 3*R # = 0
# E2 = (S^2 + (b1+b2)*S - 3*R)/2;

### Sum(x[i]*P[i]) =>
-S^3 - 2*(b1+b2)*S^2 + 6*E3 + 7*R*S + (b1+b2)*((b1+b2)*S - 3*R) # = 0
# 6*E3 = S^3 + 2*(b1+b2)*S^2 - (b1+b2)^2*S - 7*R*S + 3*(b1+b2)*R

### Eq 3:
# Squaring =>
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 # =
	b1*b2*S^2 - (b1^2 + b2^2 - b1*b2)*(b1+b2)*S - 2*b2*R*S - 2*b1*R*S + 3*R^2 + 3*(b1^2 + b2^2 - b1*b2)*R;

### Eq 3: Alternative
# [but redundant with Eq 3]
# Sum(y*...) & Sum(x^2*...) =>
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - R*S^2 - (b1 - b2)*(b1*S^2 - 2*b1*E2 + b2*E2 - R*S) +
	+ b2*E2*S - 3*b2*E3 + 2*R*E2 # = 0

### Diff =>
# [NOT used]


### Eq S:
(S^2 + 3*(b1+b2)*S - 9*R) * (S^2 - (b1 + b2)*S + 2*(b1 + b2)^2 - 6*b2*b1 - R)


### Solver:

solve.S3Ht.P2Ch2s = function(b, R, debug=TRUE) {
	b.sum = b[1] + b[2];
	# Excluded: Equal roots;
	coeff = c(1, - b.sum, 2*b.sum^2 - 6*b[1]*b[2] - R);
	S = roots(coeff);
	if(debug) print(S);
	E2 = (S^2 + b.sum*S - 3*R)/2;
	E3 = - (S^3 - 3*E2*S - R*S + b.sum*E2)/3;
	# E3 = (x.sum^3 + 2*b.sum*x.sum^2 - b.sum^2*x.sum - 7*R*x.sum + 3*b.sum*R) / 6
	# Robust: solve (x, y, z)
	x = sapply(1:length(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])) );
	x = as.vector(x);
	S = rep(S, each=3);
	yz.sum = S - x;
	yz.b.sum = R - x^2;
	if(b[1] == b[2]) {
		E3 = rep(E3, each=3)
		yz = E3 / x
		yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
		y = (yz.sum + yz.diff)/2
		z = (yz.sum - yz.diff)/2
	} else {
		# Note: order in derivation is b[1]*y + b[2]*z;
		y = (yz.sum * b[2] - yz.b.sum) / (b[2] - b[1])
		z = yz.sum - y
	}
	sol = cbind(x, y, z)
	return(sol)
}
test.S3Ht.P2SumYZ = function(sol, b, b.ext=0, R=NULL) {
	test.S3Ht.SumYZ(sol, b=b, b.ext=b.ext, R=R, n=2);
}

### Example:
b = c(1,-2)
R = 1
#
sol = solve.S3Ht.P2Ch2s(b, R)

test.S3Ht.P2SumYZ(sol, b)

### Classic Polynomial
round0.p(poly.calc(x))
#
171 + 36*x - 15*x^2 + 3*x^3 - 2*x^4 + x^5 + x^6
# Trivial solution:
# -1 - 3*x + 5*x^3 - 3*x^5 + x^6


### Ex 2:
b = c(2,-5)
R = -3
#
sol = solve.S3Ht.P2Ch2s(b, R)

test.S3Ht.P2SumYZ(sol, b)


### Test
x = sol[,1]; y = sol[,2]; z = sol[,3]
x^2 + b[1]*y + b[2]*z
y^2 + b[1]*z + b[2]*x
z^2 + b[1]*x + b[2]*y


### Classic Poly: P[6]
x = sol[1,]; b1 = b[1]; b2 = b[2];
x^6 - (b1 + b2)*x^5 + ((b1 + b2)^2 - 3*R)*x^4 +
	- (b1 + b2)*((b1 + b2)^2 - 2*b1*b2 - 2*R)*x^3 +
	+ ((b1 + b2)^4 - 2*b1*b2*(b1 + b2)^2 - 3*b1^2*b2^2 + 4*b1*b2*R - 3*(b1 + b2)^2*R + 3*R^2)*x^2 +
	- (b1 + b2)*((b1 + b2)^4 - 6*b1*b2*(b1 + b2)^2 + 9*b1^2*b2^2 +
		+ R^2 - 2*(b1 + b2)^2*R + 6*b1*b2*R)*x +
	+ ((b1 + b2)^6 - 7*b1*b2*(b1 + b2)^4 + 14*b1^2*b2^2*(b1 + b2)^2 - 6*b1^3*b2^3 +
		- R^3 + 2*(b1 + b2)^2*R^2 - 4*b1*b2*R^2 - (b1 + b2)^4*R + 11*b1^2*b2^2*R);


##############
### Extensions

# x^2 + s*x + b2*y + b1*z = R;
# s = shift of main variable;
# x^2 + s*x + b3*x*y*z + b2*y + b1*z = R;

##############
### Extension:
# x^2 + s*x + b3*x*y*z + b2*y + b1*z = R

### Sum =>
S^2 - 2*E2 + 3*b3*E3 + (s+b1+b2)*S - 3*R # = 0
# E2 = (S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R)/2;

### Sum(x[i]*P[i]) =>
S^3 + s*S^2 - 3*E2*S + 3*E3 + b3*E3*S - 2*s*E2 - R*S + (b1 + b2)*E2 # = 0

### Eq 3:
# Squaring =>
S^4 + 2*s*S^3 + (s^2 + b1^2 + 2*b1*b2 - b2^2 + 2*b3*E3 - 4*E2)*S^2 +
	+ (2*s*b3*E3 - 6*s*E2 + 4*b1*b3*E3 - 2*b1*R + 2*b1*E2 + 2*b2*R + 4*E3)*S +
	+ 2*E2^2 - 2*s^2*E2 + 4*s*b1*E2 - 4*b1*b2*E2 + 2*b2^2*E2 +
	+ 3*b3^2*E3^2 - 4*b3*E2*E3 + 6*s*E3 - 6*b1*E3 - 3*R^2 # = 0


### Eq S:
# (b3*S^3 + 3*S^2 + 9*(b1+b2+s)*S - 27*R) * P[3]
b1 = b[1]; b2 = b[2]; b3 = b[3];
b3*S^3 + (4*s*b3 - 2*(b1 + b2)*b3 - 1)*S^2 +
	+ (b1 + b2 - 3*s + 5*s^2*b3 - 5*s*(b1 + b2)*b3 + 3*(b1 + b2)^2*b3 - 7*b1*b2*b3)*S +
	- (b1 + b2 - s)^3*b3 - 2*(b1 + b2 - s)^2 +
	+ 3*b1*b2*(b1 + b2) + 2*s*(b1 + b2) + s^3*b3 - 6*s*b1*b2*b3 + 6*b1*b2 + R;


### Solver:

solve.S3Ht.P2Asym3ExtE3 = function(R, b, s, debug=TRUE) {
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	bs = b1 + b2; bd = bs - s;
	coeff = c(b3, (4*s*b3 - 2*bs*b3 - 1),
		(bs - 3*s + 5*s^2*b3 - 5*s*bs*b3 + 3*bs^2*b3 - 7*b1*b2*b3),
		(- bd^3*b3 - 2*bd^2 + 3*b1*b2*bs + 2*s*bs +
			+ s^3*b3 - 6*s*b1*b2*b3 + 6*b1*b2 + R)
	)
	x.sum = round0(roots(coeff));
	if(debug) {
		print(coeff);
		print(x.sum);
	}
	# solve (x, y, z)
	S = x.sum
	E3 = (-S^3 - 3*s*S^2 - 2*(b1+b2)*S^2 - 2*s^2*S - s*(b1+b2)*S + (b1+b2)^2*S + 7*R*S + 6*s*R - 3*(b1+b2)*R) /
		(7*b3*S - 3*b3*(b1 + b2) + 6*s*b3 - 6);
	E2 = (S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R)/2;
	#
	x = sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id])))
	x = as.vector(x)
	x.sum = rep(x.sum, each=3)
	E3 = rep(E3, each=3)
	#
	yz.sum = x.sum - x
	yz = E3 / x
	yz.b.sum = R - x^2 - s*x - b3*E3
	if(b[1] == b[2]) {
		yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
		y = (yz.sum + yz.diff)/2
		z = (yz.sum - yz.diff)/2
	} else {
		y = (yz.b.sum - b[1]*yz.sum) / (b[2] - b[1])
		z = yz.sum - y
	}
	sol = cbind(x, y, z)
	return(sol)
}
test.S3Ht.P2Asym3ExtE3 = function(sol, b, s, R=NULL) {
	b3 = b[3]; b = b[c(1,2)];
	test.S3Ht.Asym3ExtE3(sol, b=b, s=s, b3=b3, R=R, n=2);
}

### Example
R = 1
s = 1
b = c(-2,2,3)
#
sol = solve.S3Ht.P2Asym3ExtE3(R, b, s)

test.S3Ht.P2Asym3ExtE3(sol, b=b, s=s)


### Ex 2:
R = -5
s = 1
b = c(-2,2,3)
#
sol = solve.S3Ht.P2Asym3ExtE3(R, b, s)

test.S3Ht.P2Asym3ExtE3(sol, b=b, s=s)


### Test
b1 = b[1]; b2 = b[2]; b3 = b[3];
x = sol[,1]; y = sol[,2]; z = sol[,3];
x^2 + s*x + b[3]*x*y*z + b[2]*y + b[1]*z
y^2 + s*y + b[3]*x*y*z + b[2]*z + b[1]*x
z^2 + s*z + b[3]*x*y*z + b[2]*x + b[1]*y


#########################

#########################
### "Asymmetric" Variant:
### Order 3
### x[i]^3 + b1*x[j] + b2*x[k]

# x^3 + b1*y + b2*z = R
# y^3 + b1*z + b2*x = R
# z^3 + b1*x + b2*y = R

### Solution:

### Sum =>
(x^3 + y^3 + z^3) + (b1+b2)*S - 3*R # = 0
S^3 - 3*E2*S + 3*E3 + (b1+b2)*S - 3*R # = 0
# 3*E3 = -(S^3 - 3*E2*S + (b1+b2)*S - 3*R)

### Sum(x[i]*P(x)) =>
# (x^4 + y^4 + z^4) + (b1+b2)*E2 - R*S = 0
# S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 + (b1+b2)*E2 - R*S = 0

### Eq 3:
# (x^3 + b1*y + b1*z)^2 = (R + b1*z - b2*z)^2

### TODO

### Test
x^3 + b[1]*y + b[2]*z
y^3 + b[1]*z + b[2]*x
z^3 + b[1]*x + b[2]*y


#########################

#########################
#########################

#######################
### 2 Leading Terms ###
#######################

### x[i]^2 + x[j]^2 + b*x[j]

# x^2 + y^2 + b1*y = R
# y^2 + z^2 + b1*z = R
# z^2 + x^2 + b1*x = R

### Solution

### Sum =>
2*(x^2 + y^2 + z^2) + b1*(x+y+z) - 3*R # = 0
2*S^2 - 4*E2 + b1*S - 3*R
# 4*E2 = 2*S^2 + b1*S - 3*R;

### Sum(z*...) =>
x^2*z+y^2*z + y^2*x+z^2*x + x^2*y+z^2*y + b1*E2 - R*S # = 0
E2*S - 3*E3 + b1*E2 - R*S # = 0

### Diff =>
# x^2 - z^2 = -b1*(y - z)
# Note: excludes x == y == z;
### Prod =>
(x+y)*(x+z)*(y+z) - b1^3 # = 0
x^2*z+y^2*z + y^2*x+z^2*x + x^2*y+z^2*y + 2*x*y*z - b1^3 # = 0
E2*S - E3 - b1^3 # = 0
# E3 = E2*S - b1^3

### =>
E2*S - 3*E3 + b1*E2 - R*S # = 0
8*E2*S - 4*b1*E2 + 4*R*S - 12*b1^3 # = 0
2*(2*S^2 + b1*S - 3*R)*S - b1*(2*S^2 + b1*S - 3*R) + 4*R*S - 12*b1^3 # = 0
4*S^3 - (2*R + b1^2)*S - 12*b1^3 + 3*b1*R # = 0
### Eq:
(2*S - 3*b1)*(2*S^2 + 3*b1*S + 4*b1^2 - R)

### Alternatives:
### Redundant:
# Sum((x+y)*...), Sum(x*y*...);

### Alternative Eq:
# Sum(y^2*...) =>
(x^2*y^2+y^2*z^2+z^2*x^2) + (x^4+y^4+z^4) + b1*(x^3+y^3+z^3) - R*(x^2+y^2+z^2) # = 0


### Solver:

solve.2H.S3P2 = function(R, b, b.ext=0, debug=TRUE) {
	be1 = b.ext[1];
	be2 = if(length(b.ext) < 2) 0 else b.ext[2];
	# coeff = c(4, 0, - (2*R[1] + b[1]^2), -12*b[1]^3 + 3*b[1]*R[1])
	coeff = c(2, 3*b[1], 4*b[1]^2 - R[1])
	coeff = coeff + c(be2, be1, 0)
	S = round0(roots(coeff)) # numerical stability
	if(debug) print(S)
	R1 = R[1] - be1*S - be2*S^2;
	E2 = (2*S^2 + b[1]*S - 3*R1) / 4
	E3 = E2*S - b[1]^3;
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	### fast prototype
	# max.perm=0; sol = solve.EnAll(x, max.perm=max.perm, n=3)
	### robust
	S  = matrix(S, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	R1 = matrix(R1, ncol=len, nrow=3, byrow=T)
	yz.s = S - x; yz = E3 / x;
	y = 2*(R1 - x^2) - b[1]*x - yz.s^2 + 2*yz;
	y = y / b[1];
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol)
}
poly.2H.S3P2 = function(R, b, b.ext=0, max.leading=FALSE) {
	b1 = b[1]; b2 = b.ext[1];
	b3 = if(length(b.ext) > 1) b.ext[2] else 0;
	coeff = c((2+b3)^3, (2+b3)^2*(3*b1 + b2),
		(2+b3)*(-6*R + 3*b1^2 - 2*b1*b2 - b2^2 - 3*b3*R - 6*b1^2*b3 + 2*b1*b2*b3 - 6*b1^2*b3^2),
		(-12*b1*R + b1^3 - 4*b2*R - 11*b1^2*b2 - 5*b1*b2^2 - b2^3 - 6*b1*b3*R - 6*b1^3*b3 - 2*b2*b3*R +
			- 4*b1^2*b2*b3 + 2*b1*b2^2*b3 - 6*b1^3*b3^2 - 6*b1^2*b2*b3^2 + 2*b1^3*b3^3),
		(6*R^2 - 7*b1^2*R - 2*b1^4 + 6*b1*b2*R - 8*b1^3*b2 + b2^2*R + 2*b1^2*b2^2 + 3*b3*R^2 +
			8*b1^2*b3*R - 5*b1^4*b3 - 2*b1^3*b2*b3 + 3*b1^2*b2^2*b3 + 8*b1^2*b3^2*R + 7*b1^4*b3^2 +
			- 5*b1^3*b2*b3^2 + 9*b1^4*b3^3),
		(3*b1*R^2 - 2*b1^3*R - b1^5 + b2*R^2 + 8*b1^2*b2*R - 7*b1^4*b2 + 2*b1*b2^2*R + 5*b1^3*b2^2 + 3*b1^2*b2^3 +
			- 2*b1^3*b3*R + 5*b1^5*b3 + 2*b1^2*b2*b3*R - 6*b1^4*b2*b3 - 7*b1^3*b2^2*b3 - 2*b1^3*b3^2*R +
			5*b1^5*b3^2 + 11*b1^4*b2*b3^2 - 6*b1^5*b3^3),
		(-R^3 + 2*b1^2*R^2 + 3*b1^4*R + b1^6 - 2*b1*b2*R^2 + 4*b1^3*b2*R + b1^5*b2 - 3*b1^2*b2^2*R + 7*b1^4*b2^2 +
			- b1^3*b2^3 - 2*b1^2*b3*R^2 + 5*b1^4*b3*R + 6*b1^6*b3 + 3*b1^3*b2*b3*R - 12*b1^5*b2*b3 +
			2*b1^4*b2^2*b3 - 5*b1^4*b3^2*R + 17*b1^6*b3^2 - 3*b1^5*b2*b3^2 + b1^6*b3^3))
	if(max.leading) coeff else rev(coeff);
}

### Examples:

R = -2
b = 4
sol = solve.2H.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y # - R
y^2 + z^2 + b[1]*z # - R
z^2 + x^2 + b[1]*x # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.2H.S3P2(R, b)


### Extensions:

R = -2;
b = 1;
b.ext = c(1)
#
sol = solve.2H.S3P2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y + b.ext[1]*(x+y+z) # - R
y^2 + z^2 + b[1]*z + b.ext[1]*(x+y+z) # - R
z^2 + x^2 + b[1]*x + b.ext[1]*(x+y+z) # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.2H.S3P2(R, b, b.ext)
err = 1 + 2*x^2 + 2*x^3 + 3*x^4 + 2*x^5 + x^6
round0(err)


### Ext 2:
R = 1;
b = 1;
b.ext = c(0, 1)
#
sol = solve.2H.S3P2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
y^2 + z^2 + b[1]*z + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
z^2 + x^2 + b[1]*x + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.2H.S3P2(R, b, b.ext)
err = 1 + x^2 - x^3 - 2*x^4 + x^5 + x^6
round0(err)


### Ext 2 Ex 2:
R = 1;
b = 2;
b.ext = c(1, -1)
#
sol = solve.2H.S3P2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
y^2 + z^2 + b[1]*z + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
z^2 + x^2 + b[1]*x + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.2H.S3P2(R, b, b.ext)
err = 991 + 447*x - 88*x^2 - 89*x^3 + 7*x^5 + x^6
round0(err)


### Ext 2 Ex 3:
R = 6;
b = -1;
b.ext = c(3, -3)
#
sol = solve.2H.S3P2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
y^2 + z^2 + b[1]*z + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R
z^2 + x^2 + b[1]*x + b.ext[1]*(x+y+z) + b.ext[2]*(x+y+z)^2 # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.2H.S3P2(R, b, b.ext)
err = 11 + 2*x + 5*x^2 - 2*x^3 + x^6
round0(err)


### Special Cases:
### b.ext[2] = -3
b1 = b[1]; b2 = b.ext[1]; b3 = b.ext[2];
(3*R*b1^2*b2^2 + 5*R*b1^3*b2 + 57*R*b1^4 + 2*R^2*b1*b2 - 8*R^2*b1^2 + R^3 + b1^3*b2^3 +
		- b1^4*b2^2 - 10*b1^5*b2 - 109*b1^6) +
	(- 2*R*b1*b2^2 - 2*R*b1^2*b2 + 14*R*b1^3 - 3*R^2*b1 - R^2*b2 - 3*b1^2*b2^3 - 26*b1^3*b2^2 +
		- 110*b1^4*b2 - 191*b1^5)*x +
	(- 6*R*b1*b2 - 41*R*b1^2 - R*b2^2 + 3*R^2 + 7*b1^2*b2^2 + 47*b1^3*b2 + 167*b1^4)*x^2 +
	(- 6*R*b1 - 2*R*b2 + 11*b1*b2^2 + 53*b1^2*b2 + 89*b1^3 + b2^3)*x^3 +
	(3*R - 8*b1*b2 - 33*b1^2 - b2^2)*x^4 +
	- (3*b1 + b2)*x^5 + x^6
### shifted
b1 = b[1] / 2; b2 = b.ext[1] / 6; x = x - b1 - b2;
(- 324*R*b1*b2^3 - 158*R*b1^2*b2^2 - 132*R*b1^3*b2 + 851*R*b1^4 - 45*R*b2^4 + 18*R^2*b1*b2 +
		- 35*R^2*b1^2 - 3*R^2*b2^2 + R^3 + 1170*b1*b2^5 + 1905*b1^2*b2^4 + 1692*b1^3*b2^3 +
		- 4975*b1^4*b2^2 - 8238*b1^5*b2 - 9841*b1^6 + 175*b2^6) +
	16*(- 27*R*b1*b2^2 - 37*R*b1^2*b2 - 15*R*b1^3 - 6*R*b2^3 + 171*b1*b2^4 + 393*b1^2*b2^3 +
		504*b1^3*b2^2 + 331*b1^4*b2 + 51*b1^5 + 30*b2^5)*x + # (2*b2 - 3*b1)*(...)
	(- 108*R*b1*b2 - 182*R*b1^2 - 54*R*b2^2 + 3*R^2 + 1836*b1*b2^3 + 4770*b1^2*b2^2 +
		5868*b1^3*b2 + 3971*b1^4 + 387*b2^4)*x^2 +
	16*(9*b1*b2^2 + 15*b1^2*b2 + 9*b1^3 + 2*b2^3)*x^3 + # (2*b2 - 3*b1)*(...)
	3*(R - 42*b1*b2 - 49*b1^2 - 17*b2^2)*x^4 + x^6


### Debug
R = 2
b = 3
x = -3.0643873807 + 0.5677216544i
y = -0.7500000000 + 2.3196254315i
z =  1.5643873807 + 0.5677216544i
S = x+y+z; E2 = x*(y+z)+y*z; E3 = x*y*z;


#########################

### Variant:

### x[i]^2 + x[j]^2 + b*(x[i] + x[j])

# x^2 + y^2 + b1*(x+y) = R
# y^2 + z^2 + b1*(y+z) = R
# x^2 + z^2 + b1*(x+z) = R

### Solution

### trivial solution: x = y = z;

### Diff =>
# x^2 - z^2 = -b1*(x - z)
# (x-z)*(x + z + b1) = 0

### Case: x != y != z
# - has NO solutions;

### Case: x = y && x != z
# - trivial;
# x^2 + b1*x = R/2
# z^2 + b1*z = R/2 # x & z are the conjugate roots;
# also:
# x + z = - b1
# x*z = - R / 2

### Example
b = 3
R = 1
#
x = roots(c(1, b[1], -R/2))
y = x
z = x[c(2,1)]
sol = cbind(x,y,z)
sol

### Test
x^2 + y^2 + b[1]*(x+y)
y^2 + z^2 + b[1]*(y+z)
x^2 + z^2 + b[1]*(x+z)

### TODO:
# - extensions;



########################
########################

### High-Power Terms: 2
### Variant

### x[i]^2 + x[j]^2 + b*x[k]

# x^2 + y^2 + b1*z = R
# y^2 + z^2 + b1*x = R
# x^2 + z^2 + b1*y = R

# - trivial solution: x = y = z;
# - trivial system;
# - equivalent to the previous variant: + b1*(x+y);

### Solution

### Diff =>
# x^2 - z^2 = b1*(x - z)
# (x-z)*(x + z - b1) = 0

### Case: x != y != z
# - has NO solutions;

### Case x = y && x != z
# z = -x + b1;
# =>
# 2*x^2 + b1*z - R = 0
# 2*x^2 - b1*x + b1^2 - R

### Example

b = 3
R = 1
#
x = roots(c(2, - b[1], b[1]^2 - R))
y = x
z = -x + b[1]
sol = cbind(x, y, z)
sol

### Test
x^2 + y^2 + b[1]*z
y^2 + z^2 + b[1]*x
x^2 + z^2 + b[1]*y


########################
########################

########################
### Leading Terms: 3 ###

# x^2 + a1*y^2 + a2*z^2 = R
# y^2 + a1*z^2 + a2*x^2 = R
# z^2 + a1*x^2 + a2*y^2 = R

### Solution:
# - complicated solution;

### Sum =>
(a1 + a2 + 1)*(x^2 + y^2 + z^2) - 3*R # = 0
(a1 + a2 + 1)*(S^2 - 2*E2) - 3*R # = 0
# 2*(a1 + a2 + 1)*E2 = (a1 + a2 + 1)*S^2 - 3*R;

### Sum(x^2*...) =>
(x^4 + y^4 + z^4) + (a1 + a2)*((x*y)^2 + (x*z)^2 + (y*z)^2) - R*(x^2 + y^2 + z^2) # = 0
(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (a1 + a2)*(E2^2 - 2*E3*S) - R*(S^2 - 2*E2)
S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 - 2*(a1 + a2 - 2)*E3*S + 2*R*E2
# 2*(a1 + a2 - 2)*E3*S =
#   S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 + 2*R*E2

### Sum(z*...) =>
a2*(x^3 + y^3 + z^3) + (x^2*z + y^2*x + z^2*y) + a1*(y^2*z + z^2*x + x^2*y) - R*S # = 0
a2*(S^3 - 3*E2*S + 3*E3) + (E2*S - 3*E3) + (a1 - 1)*(y^2*z + z^2*x + x^2*y) - R*S # = 0
(a1 - 1)*(y^2*z + z^2*x + x^2*y) + a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S # = 0 # Eq 3a
# *(x^2*z + y^2*x + z^2*y) =>
(a1 - 1)*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*(x^2*z + y^2*x + z^2*y)
### Sum: Eq 3a + Eq 3b =>
(a1 - 1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*
	((a1 - 1)*(E2*S - 3*E3) + (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S))
(a1 - 1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*
	(a2*S^3 - 3*a2*E2*S + a1*E2*S + 3*a2*E3 - 3*a1*E3 - R*S)



### Eq:
((a1 + a2 + 1)*S^2 - 9*R)^2 * ((a1 + a2 + 1)*S^2 - R)^2 # * P0;
### P[0]
(4 - 8*a1 - 8*a2 + 6*a1*a2 - 3*a1^2*a2 - 3*a1*a2^2 + a1*a2^3 + a1^3*a2 +
	+ 9*a1^2 + 9*a2^2 - 5*a1^3 - 5*a2^3 + a1^4 + a2^4) 


### Q:
# - Do A1-type extensions have additional roots?
# - It seems NO additional roots possible!
### Technique: Sequential factorization
# - NO additional factors in this case;


### Solver:
solve.HP3.S3P2 = function(R, a, b.ext=0, debug=TRUE) {
	a.s = (a[1] + a[2] + 1);
	coeff = c(a.s, 0, -R) # only Non-equal roots!
	len = max(length(coeff), length(b.ext) + 1)
	coeff = c(rep(0, len - length(coeff)), coeff)
	b.all = c(rep(0, len - length(b.ext) - 1), rev(b.ext), 0)
	coeff = coeff + b.all;
	S = roots(coeff);
	if(debug) print(S);
	#
	pow = seq(length(b.ext));
	R1 = R[1] - sapply(S, function(S) sum(b.ext * (S^pow)));
	E2 = (a.s*S^2 - 3*R1) / a.s / 2;
	E3 = (S^4 - R1*S^2 - 4*E2*S^2 + (a.s + 1)*E2^2 + 2*R1*E2) /
		(2*(a.s - 3)*S) # TODO: a.s == 3
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	#
	sol = solve.EnAll(x, n = 3)
	return(sol)
}
test.HP3.S3P2 = function(sol, R, a, b.ext=0) {
	a = c(1, a) # TODO: include 1 in a;
	S = apply(sol, 1, sum)
	len = length(b.ext)
	ext = sapply(S, function(S) sum(b.ext * S^seq(len)));
	err = apply(sol, 1, function(sol) sum(a*sol^2))
	err = err + ext;
	return(round0(err))
}

### Examples:

R = -2
a = c(2, 3)
#
sol = solve.HP3.S3P2(R, a)
test.HP3.S3P2(sol, R, a)


### Ex 2:
R = -2
a = c(2, 3)
b.ext = c(-1)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 3:
R = -2
a = c(2, 3)
b.ext = c(-1, -1)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 4:
R = -5
a = c(2, 3)
b.ext = c(0, 10)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 5:
R = -2
a = c(2, 3)
b.ext = c(-1, -1, 2)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)



### Test
x^2 + a[2]*y^2 + a[3]*z^2 # - R
y^2 + a[2]*z^2 + a[3]*x^2 # - R
z^2 + a[2]*x^2 + a[3]*y^2 # - R

perm.gen = function(x) {
	len = length(x)
	id = seq(len)
	id.m = outer(id, id, function(i, j) ((i+j+1) %% len + 1))
	p.m = x[id.m]
	dim(p.m) = dim(id.m)
	p.m
}

R = 1;
a = c(1,2,3)
a1 = a[2]; a2 = a[3];
p.m = perm.gen(a)
d = det(p.m)

sol = solve(p.m, rep(R, 3))
x = sqrt(sol[1]); y = -x; z = -x;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


### Sum(y^2*...) =>
a1*(x^4 + y^4 + z^4) + (a2 + 1)*((x*y)^2 + (x*z)^2 + (y*z)^2) - R*(x^2 + y^2 + z^2) # = 0
### Diff: Eq2 - Eq3 =>
(a1 - 1)*((x*y)^2 + (x*z)^2 + (y*z)^2) - (a1 - 1)*(x^4 + y^4 + z^4) # = 0
### Case: a1 != 1
((x*y)^2 + (x*z)^2 + (y*z)^2) - (x^4 + y^4 + z^4) # = 0
(E2^2 - 2*E3*S) - (S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) # = 0
E2^2 + 6*E3*S + S^4 - 4*E2*S^2
# =>
(a1 + a2 - 2)*E2^2 + 3*(S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 + 2*R*E2) +
	+ (a1 + a2 - 2)*S^4 - 4*(a1 + a2 - 2)*E2*S^2
4*(a1 + a2 + 1)*E2^2 - 4*(a1 + a2 + 1)*E2*S^2 +
	+ (a1 + a2 + 1)*S^4 + 6*R*E2 - 3*R*S^2
- 2*(a1 + a2 + 1)*E2*S^2 + (a1 + a2 + 1)*S^4 - 3*R*S^2 # redundant


########################
########################

###############
### Order 3 ###
###############

#######################
### 3 Leading Terms ###

# x^3 + a1*y^3 + a2*z^3 = R
# y^3 + a1*z^3 + a2*x^3 = R
# z^3 + a1*x^3 + a2*y^3 = R

### Solution

### Sum =>
(a1 + a2 + 1)*(x^3 + y^3 + z^3) - 3*R # = 0
(a1 + a2 + 1)*(S^3 - 3*E2*S + 3*E3) - 3*R

### Sum(x^3*...) =>
(x^6 + y^6 + z^6) + (a1 + a2)*(x^3*y^3 + x^3*z^3 + y^3*z^3) - R*(x^3 + y^3 + z^3) # = 0
(- 2*E2^3 + 3*E3^2 - 12*E2*E3*S + 9*E2^2*S^2 + 6*E3*S^3 - 6*E2*S^4 + S^6) +
	(a1 + a2)*(E2^3 - 3*E3*E2*S + 3*E3^2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
(a1 + a2 - 2)*E2^3 + 3*(a1 + a2 + 1)*E3^2 - 3*(a1 + a2 + 4)*E2*E3*S + 9*E2^2*S^2 - 6*E2*S^4 +
	+ 3*R*E2*S + 6*E3*S^3 - 3*R*E3 + S^6 - R*S^3 # = 0

### Sum(y*...) =>
a1*(x^4 + y^4 + z^4) + (x^3*y + y^3*z + z^3*x) + a2*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (x^3*y + y^3*z + z^3*x) + a2*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (E2*(S^2 - 2*E2) - E3*S) + (a2 - 1)*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S +
	+ (a2 - 1)*(x*y^3 + y*z^3 + z*x^3) # = 0 ### Eq 3a
### * (x^3*y + y^3*z + z^3*x) =>
(a2 - 1)*(E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4) +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S) *
	(x^3*y + y^3*z + z^3*x) ### Eq 3b
### Eq 3a + Eq 3b =>
(a2 - 1)^2*(E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4) +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S)^2 +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S) *
	(a2 - 1)*(E2*S^2 - 2*E2^2 - E3*S)


### TODO:

### Eq:
S * ((a1 + a2 + 1)*S^3 - 27*R) * ((a1 + a2 + 1)^2 * S^6 + 27*R^2)

### Solver:
solve.3HT.S3P3 = function(R, a, b.ext=0, max.perm=1, debug=TRUE) {
	a.s = (a[1] + a[2] + 1);
	coeff = c(a.s^2, 0,0,0,0,0, 27*R^2) # only Non-equal roots!
	len = max(length(coeff), 2*length(b.ext) + 1)
	coeff = c(rep(0, len - length(coeff)), coeff)
	b.all = c(rep(0, len - length(b.ext) - 1), -2*27*R*rev(b.ext), 0)
	id2 = len - rev(2 * seq(1, length(b.ext)))
	b.all[id2] = b.all[id2] + 27*b.ext^2; # TODO: cross-products
	coeff = coeff + b.all;
	S = roots(coeff); S = c(S, 0)
	if(debug) print(S);
	#
	pow = seq(length(b.ext));
	R1 = R[1] - sapply(S, function(S) sum(b.ext * (S^pow)));
	E2 = round0(calc.E2(S, R1, a)) # !!! # TODO: has 6000 monoms & overflows !!!
	# E2[non-linearized] = P[E2^4, S^8];
	E3 = (3*R1 - a.s*(S^3 - 3*E2*S)) / 3 / a.s;
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	#
	sol = solve.EnAll(x, n = 3, max.perm=max.perm)
	return(sol)
}
calc.E2 = function(S, R, a) {
	if(a[1] != -2 || a[2] != 2) stop("works only with a = c(-2, 2)!");
	# possible coefficients: BUT possible NOT exact!
	# TODO: find exact coefficients;
	E2Subst = 108*R*S^2 - 3240*R*S^6 + 11700*R*S^10 + 56610*R^2*S^7 - 112860*R^3*S^4 - 18*S^5 +
		+ 540*S^9 - 3000*S^13;
	E2Div = - 135*R + 4050*R*S^4 - 450*R*S^8 - 1404*R^2*S - 71055*R^2*S^5 + 10395*R^3*S^2 +
		+ 45*S^3 - 1350*S^7 + 7600*S^11;
	- E2Subst / E2Div;
}

### Examples:
R = -1
a = c(-2, 2)
b.ext = c(1)
sol = solve.3HT.S3P3(R, a, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
ext = b.ext[1] * (x+y+z)
x^3 + a[1]*y^3 + a[2]*z^3 + ext # - R
y^3 + a[1]*z^3 + a[2]*x^3 + ext # - R
z^3 + a[1]*x^3 + a[2]*y^3 + ext # - R


### Debug
R = -1
a = c(-2, 2)
m = complex(re=cos(2*pi/3), im=sin(2*pi/3))
x = rootn(R / (a[1] + a[2] + 1), 3)
y = x*m; z = x*m^2;
a1 = a[1]; a2 = a[2];
S = x + y + z; E2 = x*(y+z) + y*z; E3 = x*y*z;

poly.calc(c(x+y+y, x+z+z, x+x+y, x+x+z, y+y+z, y+z+z)) * (a[1]+a[2]+1)^2


########################
########################
########################

###############
### Order 3 ###
###############

### Problems:
# - Difference works well for systems with 2 variables;
# - but it does NOT work well in systems with 3 variables (and higher power);


### x[i]^3 + b*x[i+1]

# x^3 + b*y = R
# y^3 + b*z = R
# z^3 + b*x = R

### Solution:
# Trivial solution: x = y = z;

### Case: (x, y, z) distinct

### Sum =>
S^3 - 3*E2*S + 3*E3 + b*S - 3*R # = 0;
# 6*E3 = 6*E2*S - 2*S^3 - 2*b*S + 6*R

### Sum(x[i]*...) =>
S^4 - 4*E2*S^2 + 4*E3 * S + 2*E2^2 + b*E2 - R*S # = 0;

### Sum(x[i+1]^3*...) =>
b*S^4 - R*S^3 - 4*b*E2*S^2 + 3*R*E2*S + 4*b*E3*S - 3*E3*E2*S - 3*R*E3 + E2^3 + 2*b*E2^2 + 3*E3^2

### Auxiliary Eqs:
E2x0 = 6*S^6 + 13*b*S^4 - 30*R*S^3 + 4*b^2*S^2 - 9*b*R*S + 12*b^3;
E2Div = 14*S^4 + 14*b*S^2 - 18*R*S + 3*b^2;
E2 = E2x0 / E2Div;

#######
### Eq:
# "P12" = S * P[3] * P[8]
# S * (S^3 + 9*b*x - 27*R) * P[8]
### P[8]
S^8 - 3*b*S^6 + 18*b^2*S^4 - 27*R*b*S^3 + (27*R^2 + 4*b^3)*S^2 - 27*R*b^2*S + 9*b^4


### Solver:
solve.sysHt33 = function(R, b, debug=TRUE, do.S=FALSE) {
	# only S8 used to compute the roots:
	coeff = c(1, 0, - 3*b[1], 0, 18*b[1]^2, - 27*R[1]*b[1], (27*R[1]^2 + 4*b[1]^3), - 27*R[1]*b[1]^2, 9*b[1]^4)
	S = roots(coeff);
	isEq = round0(S^3 + 9*b[1]*S - 27*R) == 0
	S = S[ ! isEq]
	if(debug) print(S);
	E2x0 = 6*S^6 + 13*b*S^4 - 30*R*S^3 + 4*b^2*S^2 - 9*b*R*S + 12*b^3;
	E2Div = 14*S^4 + 14*b*S^2 - 18*R*S + 3*b^2;
	E2 = E2x0 / E2Div;
	E3 = (6*E2*S - 2*S^3 - 2*b[1]*S + 6*R[1]) / 6;
	x = sapply(1:length(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	y = (R - x^3) / b[1]
	z = (R - y^3) / b[1]
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	if(do.S) sol = cbind(sol, S=as.vector(S));
	return(sol);
}

###########
### Example
R = 1
b = 2
#
sol = solve.sysHt33(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b[1]*y
y^3 + b[1]*z
z^3 + b[1]*x


### Classic Polynomial: P[24]
(- R[1]^2*b[1]^9 + 3*R[1]^4*b[1]^6 - 3*R[1]^6*b[1]^3 + R[1]^8 + b[1]^12) +
(- R[1]*b[1]^10 + 3*R[1]^3*b[1]^7 - 3*R[1]^5*b[1]^4 + R[1]^7*b[1])*x +
(3*R[1]^2*b[1]^8 - 3*R[1]^4*b[1]^5 + R[1]^6*b[1]^2 - b[1]^11)*x^2 +
(2*R[1]*b[1]^9 - 9*R[1]^3*b[1]^6 + 16*R[1]^5*b[1]^3 - 8*R[1]^7)*x^3 +
(- 6*R[1]^2*b[1]^7 + 13*R[1]^4*b[1]^4 - 7*R[1]^6*b[1] + b[1]^10)*x^4 +
(- 3*R[1]*b[1]^8 + 10*R[1]^3*b[1]^5 - 6*R[1]^5*b[1]^2)*x^5 +
(10*R[1]^2*b[1]^6 - 35*R[1]^4*b[1]^3 + 28*R[1]^6 - b[1]^9)*x^6 +
(4*R[1]*b[1]^7 - 22*R[1]^3*b[1]^4 + 21*R[1]^5*b[1])*x^7 +
(- 12*R[1]^2*b[1]^5 + 15*R[1]^4*b[1]^2 + b[1]^8)*x^8 +
(- 5*R[1]*b[1]^6 + 40*R[1]^3*b[1]^3 - 56*R[1]^5)*x^9 +
(18*R[1]^2*b[1]^4 - 35*R[1]^4*b[1] - b[1]^7)*x^10 +
(6*R[1]*b[1]^5 - 20*R[1]^3*b[1]^2)*x^11 +
(- 25*R[1]^2*b[1]^3 + 70*R[1]^4 + b[1]^6)*x^12 +
(- 7*R[1]*b[1]^4 + 35*R[1]^3*b[1])*x^13 +
+ (15*R[1]^2*b[1]^2 - b[1]^5)*x^14 +
+ (8*R[1]*b[1]^3 - 56*R[1]^3)*x^15 +
+ (- 21*R[1]^2*b[1] + b[1]^4)*x^16 +
+ (- 6*R[1]*b[1]^2)*x^17 +
+ (28*R[1]^2 - b[1]^3)*x^18 +
+ 7*R[1]*b[1]*x^19 +
+ b[1]^2*x^20 +
- 8*R[1]*x^21 +
- b[1]*x^22 + x^24


############################

############################
### x[i]^3 + b*(x[j] + x[k])

# x^3 + b1*(y+z) = R
# y^3 + b1*(x+z) = R
# z^3 + b1*(x+y) = R

# Trivial solution: x = y = z;
# Trivial system;

### Solution

### Diff =>
# x^3 - y^3 = b1*(x-y)
# x^3 - z^3 = b1*(x-z)
# y^3 - z^3 = b1*(y-z)

# (x-y)*(x^2 + y^2 + x*y - b1) = 0

# Case: x != y != z
# x^2 + y^2 + x*y - b1 = 0
# x^2 + z^2 + x*z - b1 = 0
# y^2 + z^2 + y*z - b1 = 0
### Sum =>
# E2 = 2/3*S^2 - b1
### Diff =>
# y^2 - z^2 + x*(y - z) = 0
# Case: x != y != z =>
# x + y + z = 0
# => E2 = -b1

### Sum =>
# x^3 + y^3 + z^3 + 2*b1*S = 3*R, where S = 0
# x^3 + y^3 + z^3 = 3*R

# x^3 + y^3 + z^3 = S^3 - 3*E2*S + 3*E3
# =>
# E3 = R


### Example 1:
b = 3
R = 1
#
x = roots(c(1,0, -b[1], -R[1]))
y = as.vector(sapply(x, function(x) roots(c(1, x, x^2 - b[1]))))
x = rep(x, each=2)
z = -x-y
sol = cbind(x,y,z)
sol

### Test
x^3 + b[1]*(y+z)
y^3 + b[1]*(x+z)
z^3 + b[1]*(x+y)

# trivial:
(x^3 - b*x - R)^2


##########################
##########################

### TODO:
# - move to appropriate section;

### Shifted
### x[i]^2 * (x[i] - shift) + b*Sum

# x^2*(x - s) + b1*(x+y+z) = R
# y^2*(y - s) + b1*(x+y+z) = R
# z^2*(z - s) + b1*(x+y+z) = R

### Solution

# - trivial solution: x = y = z;

### Diff =>
# x^2*(x-s) = y^2*(y-s)
# y^2*(y-s) = z^2*(z-s)
# z^2*(z-s) = x^2*(x-s)
# =>
# x^3 - y^3 - s*(x^2 - y^2) = 0
# (x-y)*(x^2 + y^2 + x*y - s*(x+y)) = 0

# Case: x != y != z
# x^2 + y^2 + x*y - s*(x+y) = 0
# x^2 + z^2 + x*z - s*(x+z) = 0
# y^2 + z^2 + y*z - s*(y+z) = 0
# =>
# y^2 + x*y - s*(x+y) = z^2 + x*z - s*(x+z)
# y^2 - z^2 + x*(y-z) - s*(y-z) = 0
# (y-z)*(x + y + z - s) = 0
# x + y + z = s

### Sum =>
# x^3 + y^3 + z^3 - s*(x^2 + y^2 + z^2) + 3*b1*S = 3*R
# x^3 + y^3 + z^3 = S^3 - 3*E2*S + 3*E3 =>
# S^3 - 3*E2*S + 3*E3 - s*(S^2 - 2*E2) + 3*b1*S - 3*R = 0
# S = s =>
# - s*E2 + 3*b1*s - 3*R + 3*E3 = 0
# E3 =  s*E2/3 - b1*s + R

# E2 = 2/3*S^2 - 2/3*s*S
# E2 = 0
# E3 = s*E2/3 - b*s + R
# E3 = - b1*s + R

### Example:

b = 3
R = 1
s = 2
#
x = roots(c(1, -s, 0, b[1]*s - R))
yz.sum = s - x
yz = yz.sum^2 - s*yz.sum
yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
y = (yz.sum + yz.diff)/2
z = (yz.sum - yz.diff)/2
sol = cbind(x,y,z)
sol

### Test
x^2*(x - s) + b[1]*(x+y+z)
y^2*(y - s) + b[1]*(x+y+z)
z^2*(z - s) + b[1]*(x+y+z)

### Classical Polynomial
# x => trivial polynomial: P3;
# y, z => (P3)^2

round0.p(poly.calc(sol[,2:3]))
round0.p(poly.calc(sol[,2:3] - s))

# (x^3 - 2*x^2 + 5)^2
err = 25 - 20*x^2 + 10*x^3 + 4*x^4 - 4*x^5 + x^6
round0(err)
x = x - s # shift back;
err = 25 + 40*x + 56*x^2 + 42*x^3 + 24*x^4 + 8*x^5 + x^6
round0(err)



########################

###############
### Order 3 ###

### x^3 + b3*x*y*z + b[j]*Sum[j]

# x^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R
# y^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R
# z^3 + b3*x*y*z + b2*(x^2+y^2+z^2) + b1*(x+y+z) = R

### TODO: extend with + b4*(x*y*z)^2;

m3 = unity(3, all=F)

### Solution:

# Trivial system;

### Diff =>
# x^3 - y^3 = 0
# x^3 - z^3 = 0
# Case: x != y != z =>
# y = x*m3
# z = x*m3^2
# =>
# S = 0;
# x^2 + y^2 + z^2 = x^2 *(1 + m + m^2) = 0;

### =>
# x^3 + b3*x^3 = R
# x^3 = R / (b3 + 1)

### TODO


###########################
###########################

###########################
### x^3 + b2*x^2*y*z + b1*x

# x^3 + b2*x^2*y*z + b1*x = R
# y^3 + b2*x*y^2*z + b1*y = R
# z^3 + b2*x*y*z^2 + b1*z = R

### Solution:

### Diff =>
# x^3 - y^3 + b2*x*y*z*(x-y) + b1*(x-y) = 0
(x-y)*(x^2 + y^2 + x*y + b2*x*y*z + b1) # = 0
(y-z)*(y^2 + z^2 + y*z + b2*x*y*z + b1) # = 0
(z-x)*(x^2 + z^2 + x*z + b2*x*y*z + b1) # = 0
# Case: x != y != z =>
x^2 + y^2 + x*y + b2*x*y*z + b1 # = 0
y^2 + z^2 + y*z + b2*x*y*z + b1 # = 0
x^2 + z^2 + x*z + b2*x*y*z + b1 # = 0
# => Diff =>
# x + y + z = 0

### Sum =>
# x^3 + y^3 + z^3 + b2*x*y*z*(x+y+z) + b1*(x+y+z) = 3*R
# x^3 + y^3 + z^3 = 3*R
# (x+y+z)^3 - E2*(x+y+z) - 3*x*y*z = 3*R
# x*y*z = -R

### Sum(x^2) =>
# 2*(x^2 + y^2 + z^2) + E2 + 3*b2*E3 + 3*b1 = 0
# 2*(x+y+z)^2 - 4*E2 + E2 - 3*b2*R + 3*b1 = 0
# E2 + b2*R - b1 = 0
# E2 = -b2*R + b1

### Example
b = c(1, 3)
R = 1
#
x = roots(c(1,0, -b[2]*R + b[1], R))
yz = -R / x;
yz.sum = -x;
yz.diff = sqrt(yz.sum^2 - 4*yz + 0i)
y = (yz.sum + yz.diff)/2
z = (yz.sum - yz.diff)/2
sol = cbind(x, y, z)
sol = rbind(sol, sol[,c(2,3,1)], sol[,c(3,1,2)]) # add all variants

### Test
x^3 + b[2]*x^2*y*z + b[1]*x
y^3 + b[2]*x*y^2*z + b[1]*y
z^3 + b[2]*x*y*z^2 + b[1]*z

### Classic Polynomial

# Trivial: P9 = (P3)^3

### TODO


#############################
#############################

#############################
### Univariate Side-Chain ###
### Higher Powers         ###
#############################

### x[i]^3 + b2*x[j]^2 + b1*x[j]

# x^3 + b2*y^2 + b1*y = R
# y^3 + b2*z^2 + b1*z = R
# z^3 + b2*x^2 + b1*x = R


### Solution:

### Sum =>
S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 + 3*E3 - 3*R # = 0
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)

### Sum( (b[1]*x[i] + b[2]*x[i]^2) * ...) =>
b2*S^5 + b1*S^4 - 5*b2*E2*S^3 + 5*b2*E3*S^2 - 5*b2*E2*E3 + 4*b1*E3*S - 2*b2^2*E3*S +
	- 4*b1*E2*S^2 + 5*b2*E2^2*S - 2*b1*b2*E2*S +
	+ 2*b1*E2^2 + b2^2*E2^2 + b1^2*E2 - 2*b1*b2^2*E2 +
	+ b1*b2*S^3 + b1*b2^2*S^2 + b1^2*b2*S +
	- b2*R*S^2 + 2*b2*R*E2 - b1*R*S - 3*b1*b2*R # = 0

### Sum(y^3*...) =>
E2^3 - 3*E3*E2*S + 3*E3^2 +
	+ b2*S^5 - 5*b2*E2*S^3 + 5*b2*E3*S^2 + 5*b2*E2^2*S - 5*b2*E2*E3 +
	+ b1*S^4 - 4*b1*E2*S^2 + 4*b1*E3*S + 2*b1*E2^2 +
	- R*S^3 + 3*R*E2*S - 3*R*E3 # = 0


### Eq: P[3] * P[8]
# [the original eq. had 730 monomials]
S^8 - b2*S^7 + (2*b2^2 - 3*b1)*S^6 + (8*b2^3 + 12*b2*b1)*S^5 +
	+ (-b2^4 + 27*b2^2*b1 + 18*b1^2)*S^4 +
	- (27*R*b2^2 + 5*b2^5 + 27*R*b1 + 9*b2^3*b1 - 24*b2*b1^2)*S^3 +
	+ (27*R^2 + 14*R*b2^3 + 7*b2^6 - 36*R*b2*b1 - 18*b2^4*b1 - 10*b2^2*b1^2 + 4*b1^3)*S^2 +
	+ (-27*R^2*b2 - 2*R*b2^4 + 27*R*b2^2*b1 + 21*b2^5*b1 - 27*R*b1^2 - 26*b2^3*b1^2 + 2*b2*b1^3)*S +
	(27*R^2*b2^2 - 7*R*b2^5 - 18*R*b2^3*b1 + 27*R*b2*b1^2 + 14*b2^4*b1^2 - 20*b2^2*b1^3 + 9*b1^4)

### Special Case:
# for b2 == 1; b1 == 0;
R*(27*R - 7) - R*(27*R + 2)*S + (27*R^2 + 14*R + 7)*S^2 - (27*R + 5)*S^3 +
	- S^4 + 8*S^5 + 2*S^6 - S^7 + S^8
# TODO:
# - simplified solution: for b2 == 1; b1 == 0; (if it is useful?)


### Solver:
solve.Y2Y1.S3P3 = function(R, b=c(1,0), debug=TRUE) {
	b2 = b[1]; b1 = if(length(b) > 1) b[2] else 0;
	coeff = c(1, - b2, (2*b2^2 - 3*b1), (8*b2^3 + 12*b2*b1),
		(-b2^4 + 27*b2^2*b1 + 18*b1^2),
		- (27*R*b2^2 + 5*b2^5 + 27*R*b1 + 9*b2^3*b1 - 24*b2*b1^2),
		(27*R^2 + 14*R*b2^3 + 7*b2^6 - 36*R*b2*b1 - 18*b2^4*b1 - 10*b2^2*b1^2 + 4*b1^3),
		(-27*R^2*b2 - 2*R*b2^4 + 27*R*b2^2*b1 + 21*b2^5*b1 - 27*R*b1^2 - 26*b2^3*b1^2 + 2*b2*b1^3),
		(27*R^2*b2^2 - 7*R*b2^5 - 18*R*b2^3*b1 + 27*R*b2*b1^2 + 14*b2^4*b1^2 - 20*b2^2*b1^3 + 9*b1^4));
	# [OLD] simple coeffs
	# coeff = c(1, -1, 2, 8, -1, - (27*R + 5),
	#	(27*R^2 + 14*R + 7), - R*(27*R + 2), R*(27*R - 7))
	S = roots(coeff)
	if(debug) print(S);
	R = R[1] - 0*S; # extensions
	E2 = E2.S3P3.UnivSCh(S, R, b);
	E3 = - (S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R) / 3;
	# solve: x
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	R = rep(R, each=3); S = rep(S, each=3); E3 = rep(E3, each=3);
	# robust
	yz.s = S - x; yz = E3 / x;
	y_x  = b1*x^3 - b2^2*(yz.s^3 - 3*yz.s*yz) - b2^3*x^2 - b1*b2^2*x - (b1 - b2^2)*R;
	ydiv =  b2*(x^3 - R) - b1^2; # TODO: DIV0!
	y = y_x / ydiv;
	z = yz.s - y;
	return(cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z)));
}
E2.S3P3.UnivSCh = function(S, R, b) {
	# TODO: re-order b!
	b1 = b[2]; b2 = b[1]; R = R[1];
	E2Subst = 1701*R*S*b1*b2^4 - 243*R*S*b1^3 - 756*R*S*b2^6 + 918*R*S^2*b1*b2^3 - 243*R*S^2*b1^2*b2 +
		+ 288*R*S^2*b2^5 - 1782*R*S^3*b1*b2^2 + 648*R*S^3*b1^2 + 450*R*S^3*b2^4 - 486*R*S^4*b1*b2 +
		- 945*R*S^4*b2^3 - 702*R*S^5*b2^2 - 1134*R*b1*b2^5 + 1620*R*b1^2*b2^3 - 729*R*b1^3*b2 +
		+ 729*R^2*S*b1*b2 - 486*R^2*S*b2^3 + 972*R^2*S^2*b2^2 - 729*R^2*b1*b2^2 + 378*S*b1^2*b2^5 +
		- 540*S*b1^3*b2^3 + 243*S*b1^4*b2 + 630*S^2*b1*b2^6 - 1233*S^2*b1^2*b2^4 + 504*S^2*b1^3*b2^2 +
		- 537*S^3*b1*b2^5 - 180*S^3*b1^2*b2^3 + 90*S^3*b1^3*b2 + 252*S^3*b2^7 - 516*S^4*b1*b2^4 +
		+ 540*S^4*b1^2*b2^2 - 189*S^4*b1^3 + 30*S^4*b2^6 + 576*S^5*b1*b2^3 - 27*S^5*b1^2*b2 - 144*S^5*b2^5 +
		+ 513*S^6*b1*b2^2 - 108*S^6*b1^2 + 186*S^6*b2^4 + 45*S^7*b1*b2 + 297*S^7*b2^3 + 90*S^8*b2^2;
	E2Div = 1539*R*S*b1*b2^2 - 486*R*S*b1^2 - 459*R*S*b2^4 + 648*R*S^2*b1*b2 + 702*R*S^2*b2^3 + 810*R*S^3*b2^2 +
		- 567*R*b1*b2^3 + 189*R*b2^5 - 729*R^2*b2^2 + 819*S*b1*b2^5 - 45*S*b1^2*b2^3 - 108*S*b1^3*b2 +
		- 504*S*b2^7 + 954*S^2*b1*b2^4 - 621*S^2*b1^2*b2^2 + 108*S^2*b1^3 - 123*S^2*b2^6 - 810*S^3*b1*b2^3 +
		+ 162*S^3*b1^2*b2 + 312*S^3*b2^5 - 945*S^4*b1*b2^2 + 270*S^4*b1^2 - 369*S^4*b2^4 - 108*S^5*b1*b2 +
		- 684*S^5*b2^3 - 225*S^6*b2^2 - 756*b1*b2^6 + 1458*b1^2*b2^4 - 1026*b1^3*b2^2 + 243*b1^4;
	E2 = - E2Subst / E2Div;
	return(E2);
}

### Examples
R = 2
b = c(-1, 3)
sol = solve.Y2Y1.S3P3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b[1]*y^2 + b[2]*y # - R
y^3 + b[1]*z^2 + b[2]*z # - R
z^3 + b[1]*x^2 + b[2]*x # - R


### Classic Polynomial
# - see file: Poly.System.S3P3.Mega.ClassicPoly.R;


########################

########################
### Mixed Side-Chain ###
########################

### x[i]^3 + b2*x[j]^2 + b1*x[k]

# x^3 + b2*y^2 + b1*z = R
# y^3 + b2*z^2 + b1*x = R
# z^3 + b2*x^2 + b1*y = R

### Solution:

### Sum =>
S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 + 3*E3 - 3*R # = 0
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)

### Eq2:
# Sum(Tr(Eq)^2) =>
- 3*E3^2 + 12*E2*E3*S - 6*E3*S^3 + 4*b2^2*E3*S + 2*b1*E3*S - 6*b1*b2*E3 +
	+ 2*E2^3 - 9*E2^2*S^2 + 2*b2^2*E2^2 + 4*b1*E2^2 +
	+ 6*E2*S^4 - 4*b2^2*E2*S^2 - 2*b1*E2*S^2 + 6*b1*b2*E2*S + 4*b2*R*E2 +
	- S^6 + b2^2*S^4 - 2*b1*b2*S^3 - 2*b2*R*S^2 - b1^2*S^2  + 2*b1*R*S + 3*R^2 # = 0

### Eq3:
# Sum(Tr(Eq)^2) & Diff(Eq 2 - Eq 3) =>
b2*E2*E3 + 2*b2*E3*S^2 - 2*b2^2*E3*S - b1*E3*S +
	- b2*E2^2*S - 3*b2^2*E2^2 - 2*b1*E2^2 +
	+ 4*b2^2*E2*S^2 + b1*E2*S^2 - b1^2*E2 - 4*b2*R*E2 +
	- b2^2*S^4 + b1^2*S^2 + 2*b2*R*S^2 - 2*b1*R*S # = 0

### Eq S:
S^8 - b2*S^7 - 3*b1*S^6 + 2*b2^2*S^6 - 9*b1*b2*S^5 + 8*b2^3*S^5 + (18*b1^2 - b2^4)*S^4 +
	- (5*b2^5 + b1*b2^3 + 27*b1*R + 27*b2^2*R)*S^3 +
	+ (7*b2^6 + 11*b1*b2^4 + 4*b1^3 + 14*b2^3*R - 7*b1^2*b2^2 + 45*b1*b2*R + 27*R^2)*S^2 +
	- (7*b1*b2^5 + 2*b2^4*R + b1^3*b2 + 9*b1^2*b2^3 + 27*b1^2*R + 18*b1*b2^2*R + 27*b2*R^2)*S +
	- 7*b2^5*R + 9*b1^4 + 14*b1^2*b2^4 + 19*b1^3*b2^2 - 27*b1*b2^3*R + 27*b2^2*R^2

### Auxiliary Eqs:
# E2 = ...; see file PP.S3P3.MixedSCh.R;
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)


### Solver:
solve.S3P3.MixedSCh = function(R, b, debug=TRUE) {
	coeff = coeff.S3P3.MixedSideChain(R, b);
	S = roots(coeff);
	if(debug) print(S);
	# function E2.S3P3.MixedSideChain(S, R, b):
	# is defined in file: PP.S3P3.MixedSCh.R;
	E2 = E2.S3P3.MixedSideChain(S, R, b);
	b1 = b[1]; b2 = b[2]; R = R[1];
	E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R) / 3;
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	yz.s = S - x; yz = E3 / x;
	# robust
	yz.d = R - x^3 + b2*x^2 - b1*x - b2*(yz.s^2 - 2*yz);
	yz.d = yz.d / (yz.s^2 - yz - b1);
	y = (yz.s + yz.d)/2;
	z = (yz.s - y);
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol);
}
coeff.S3P3.MixedSideChain = function(R, b) {
	if(length(b) < 2) {
		print("Warning: Missing b!");
		b = c(b, 0);
	}
	b1 = b[1]; b2 = b[2]; R = R[1];
	# formula for coeff's is in reverse order;
	coeff = c(
		- 7*b2^5*R + 9*b1^4 + 14*b1^2*b2^4 + 19*b1^3*b2^2 - 27*b1*b2^3*R + 27*b2^2*R^2,
		- 7*b1*b2^5 - 2*b2^4*R - b1^3*b2 - 9*b1^2*b2^3 - 27*b1^2*R - 18*b1*b2^2*R - 27*b2*R^2,
		7*b2^6 + 11*b1*b2^4 + 4*b1^3 + 14*b2^3*R - 7*b1^2*b2^2 + 45*b1*b2*R + 27*R^2,
		- 5*b2^5 - b1*b2^3 - 27*b1*R - 27*b2^2*R,
		18*b1^2 - b2^4,
		- 9*b1*b2 + 8*b2^3,
		- 3*b1 + 2*b2^2, - b2, 1
	);
	return(rev(coeff));
}

### Examples:

R = -1;
b = c(3,-1);
sol = solve.S3P3.MixedSCh(R, b);
x = sol[,1]; y = sol[,2]; z= sol[,3];

### Test
x^3 + b[2]*y^2 + b[1]*z # - R
y^3 + b[2]*z^2 + b[1]*x # - R
z^3 + b[2]*x^2 + b[1]*y # - R


##########################
##########################

##########################
### Difference Types   ###

# quasi-[Negative Correlations]

# x^2 - y^2 + b*x*y = R
# y^2 - z^2 + b*y*z = R
# z^2 - x^2 + b*x*z = R

### Solution:
# - moved to file:
#   Poly.System.Hetero.Symmetric.S3.Diff.R;


########################
########################
########################

#########################
### Type: x^n + b*x*z ###
#########################

###############
### Order 2 ###
###############

# x^2 + b*x*z = R
# y^2 + b*y*z = R
# z^2 + b*x*y = R

### Solution:

### Diff: (1) - (2) =>
x^2 - y^2 + b*z*(x - y) # = 0
(x - y)*(x + y + b*z) # = 0
### Cases:
# x == y *OR*
# x + y + b*z = 0;
x + y + b*z # = 0

### Sum => [not used]
x^2 + y^2 + z^2 + b*E2 - 3*R # = 0
S^2 + (b - 2)*E2 - 3*R # = 0
# (b - 2)*E2 = - (S^2 - 3*R);

### Sum: (1) + (2) =>
x^2 + y^2 + b*z*(x+y) - 2*R # = 0
(x+y)^2 - 2*x*y + b*z*(-b*z) - 2*R # = 0
(-b*z)^2 - 2*x*y + b*z*(-b*z) - 2*R # = 0
x*y + R # = 0
# =>
# z^2 = R*(b + 1);
# =>
# x + y = - b*z;

### Solution:
R = 1
b = 2
#
z = sqrt(R[1]*(b[1] + 1) + 0i)
z = c(z, -z);
xy.s = - b[1] * z;
xy = -R[1];
xy.d = sqrt(xy.s^2 - 4*xy)
xy.d = c(xy.d, -xy.d)
xy.s = rep(xy.s, 2); z = rep(z, 2)
x = (xy.s + xy.d) / 2;
y = (xy.s - xy.d) / 2;

### Test
x^2 + b[1]*x*z # - R
y^2 + b[1]*y*z # - R
z^2 + b[1]*x*y # - R


########################
########################
########################

