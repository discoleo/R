
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogenous Symmetric
###
### draft v.0.4a


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
# - decomposing system into lower order system;
# - the trivial solution x = y = z will not be covered here;

### Sum =>
# (x^n + y^n + z^n) + b*S = 3*R

### Sum(x[i] * P[i]) =>
# (x^(n+1) + y^(n+1) + z^(n+1)) + b*E2 = R*S

###  Eq3:
# x^n = R - b*y
# x^(2*n) = R^2 + b^2*y^2 - 2*b*R*y
### [3a] Sum =>
# x^(2*n) + y^(2*n) + z^(2*n) = 3*R^2 + b^2*S^2 - 2*b^2*E2 - 2*b*R*S;

### Eq3-variant:
# x^n - b*x = R - b*y - b*x
# x^(2*n) - 2*b*x^(n+1) + b^2*x^2 = R^2 + b^2*y^2 + b^2*x^2 - 2*b*R*y - 2*b*R*x + 2*b^2*x*y
### [3b] Sum =>
# x^(2*n) + y^(2*n) + z^(2*n) - 2*b*(x^(n+1) + y^(n+1) + z^(n+1)) + b^2*S^2 - 2*b^2*E2 =
#  = 3*R^2 + 2*b^2*(S^2 - 2*E2) - 4*b*R*S + 2*b^2*E2
### Sum[3b] - Sum[3a] + 2*b*Sum[2] =>
# b^2*S^2 - 2*b^2*E2 + 2*b^2*E2 = 2*b*R*S + 3*R^2 + 2*b^2*(S^2 - 2*E2) - 4*b*R*S + 2*b^2*E2 - (3*R^2 + b^2*S^2 - 2*b^2*E2 - 2*b*R*S)
# 0 == 0; # [equations are correlated]


# Note:
# (x^n + y^n + z^n) as well as the (n+1) variant:
# - can be decomposed as polynomials of: S, E2, E3;,
#   where S = x + y + z;
# - S, E2, E3 = elementary polynomials; 
# n:   (x^n + y^n + z^n) = D1(S^n, E2, E3);
# n+1: (x^n + y^n + z^n) = D2(S^(n+1), E2, E3);

# Alternative to Eq3: for very low orders;
### Diff =>
# x^n - y^n = b*(z - y)
### Prod =>
# (x^n - y^n)*(y^n - z^n)*(z^n - x^n) = b^3*(z-y)*(x-z)*(y-x)

### Complexity:
# - initial System: => order P[n^3];
#  -- polynomial can be decomposed = P[n]*P[n^3 - n];
# - decomposed system:
#  -- D(S, E2, E3): orders of E2 & E3 are usually much lower;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R


################################
################################


###############
### Order 2 ###

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

### Method 1:
### Sum =>
# x^2 + y^2 + z^2 + b1*(x+y+z) = 3*R
# S^2 - 2*E2 + b1*S - 3*R = 0;
# 2*E2 = S^2 + b1*S - 3*R;

### Sum(x[i]*...) =>
# x^3 + y^3 + z^3 + b1*E2 = R*S
# S^3 - 3*E2*S + 3*E3 + b1*E2 = R*S
# 2*S^3 - 6*E2*S + 6*E3 + 2*b1*E2 - 2*R*S
# 2*S^3 - 3*(S^2 + b1*S - 3*R)*S + 6*E3 + b1*(S^2 + b1*S - 3*R) - 2*R*S
6*E3 -S^3 - 2*b1*S^2 + 7*R*S + b1^2*S - 3*b1*R
# 6*E3 = S^3 + 2*b1*S^2 - 7*R*S - b1^2*S + 3*b1*R

### Sum(x[i+1]^2*...) =>
E2^2 - 2*E3*S + b1*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2)
(S^2 + b1*S - 3*R)^2 - 8*E3*S + b1*(4*S^3 - 6*(S^2 + b1*S - 3*R)*S + 12*E3) - 4*R*(S^2 - (S^2 + b1*S - 3*R))
(S^2 + b1*S - 3*R)^2 - 8*E3*S + b1*(-2*S^3 - 6*b1*S^2 + 22*R*S + 12*E3) - 12*R^2
S^4 - 6*R*S^2 - 5*b1^2*S^2 + 16*b1*R*S - 8*E3*S + 12*b1*E3 - 3*R^2
S^4 + 2*b1*S^3 - (10*R + b1^2)*S^2 + 6*(b1*R + b1^3)*S - 18*b1^2*R + 9*R^2
### Eq:
(S^2 + 3*b1*S - 9*R)*(S^2 - b1*S - R + 2*b1^2)


### [old/unstable]
### Diff =>
# x^2 - y^2 = b1*(z-y)
# y^2 - z^2 = b1*(x-z)
# z^2 - x^2 = b1*(y-x)
# Prod =>
# (x+y)*(x+z)*(y+z) = (-1)*b1^3;
# S^3 - (x^3 + y^3 + z^3) + 3*b1^3 = 0;

# [Prod] =>
# (x^2 + x*y + x*z + y*z)*(y+z) = - b1^3
# S[x^2*y] + 2*x*y*z + b1^3 = 0;
# E2*S - 3*E3 + 2*E3 + b1^3 = 0;
# E2*S - E3 + b1^3 = 0

### Alternative:
### Method 2: classic
# b1*y = R - x^2
# b1^3*z = b1^2*R - (R - x^2)^2
# b1^3*z = b1^2*R - R^2 - x^4 + 2*R*x^2
# =>
# x^8 - 4*R*x^6 + (6*R^2 - 2*b1^2*R)*x^4 + 4*R^2*(b1^2 - R)*x^2 + b[1]^7*x + (b1^2*R - R^2)^2 - b[1]^6*R = 0
# (x^2 + b1*x - R) * P6;
- R*b[1]^4 + 2*R^2*b[1]^2 - R^3 + b[1]^6 + (2*R*b[1]^3 - R^2*b[1] - b[1]^5)*x + (3*R^2 - 3*R*b[1]^2 + b[1]^4)*x^2 +
	+ (2*R*b[1] - b[1]^3)*x^3 - (3*R - b[1]^2)*x^4 - b[1]*x^5 + x^6

### Classic Polynomial: Shifted
# b10 = b1 / 6; x0 = x - b1;
(- 963*R*b1^4 + 69*R^2*b1^2 - R^3 + 39991*b1^6) +
	(240*R*b1^3 - 5712*b1^5)*x + (- 90*R*b1^2 + 3*R^2 + 819*b1^4)*x^2 +
	(- 112*b1^3)*x^3 + (- 3*R + 21*b1^2)*x^4 + x^6


### Solution:
solve.sysHt32 = function(R, b, doPrint=TRUE) {
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1: power 1;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A1: power 2;
	# coeff = c(1, 2*b[1], - (10*R[1] + b[1]^2), 6*(b[1]*R[1] + b[1]^3), - 18*b[1]^2*R[1] + 9*R[1]^2)
	# if(b2 != 0) coeff = coeff + c(0, 10*b2, -6*b[1]*b2 + 9*b2^2, 18*b[1]^2*b2 - 18*R[1]*b2, 0)
	coeff = c(1 + b3, b2 - b[1], - R[1] + 2*b[1]^2)
	S = roots(coeff)
	if(doPrint) print(S)
	len = length(S);
	if(len == 0) stop("NO solutions!")
	R1 = R[1] - b2*S - b3*S^2;
	# [REMOVED] remove x == y == z = S / 3;
	# as it causes numerical instability due to roots multiplicity;
	# isEq = round0(S^2 + 3*b[1]*S - 9*R1) == 0
	# if(any(isZero)) print("Warning: f(S) == 0!")
	# S = S[ ! isEq]; R1 = R1[ ! isEq];
	E2 = round0(S^2 + b[1]*S - 3*R1)/2
	E3 = - (S^3 - 3*E2*S + b[1]*E2 - R1*S) / 3
	E3 = round0(E3, tol=1E-10); # improve numerics when E3 == 0;
	x = sapply(1:length(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	R1 = matrix(R1, ncol=len, nrow=3, byrow=TRUE)
	y = (R1 - x^2)/b[1]
	z = (R1 - y^2)/b[1]
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	sol
}

### Example:
R = 3
b = -1
#
sol = solve.sysHt32(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^2 + b[1]*y
y^2 + b[1]*z
z^2 + b[1]*x

### Classic Polynomial: P8 or P6 (when S == 0)
round0.p(poly.calc(sol[,1]))


#########
### Ex 2:
R = 3
b = 3
#
sol = solve.sysHt32(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^2 + b[1]*y 
y^2 + b[1]*z
z^2 + b[1]*x

round0.p(poly.calc(sol[,1]))

621 - 108*x + 27*x^2 - 9*x^3 - 3*x^5 + x^6


### alternative / classic
coeff = c(1,0, - 4*R,0, (6*R^2 - 2*b[1]^2*R), 0, 4*R^2*(b[1]^2 - R), b[1]^7, (b[1]^2*R - R^2)^2 - b[1]^6*R)
x = roots(coeff)
y = (R - x^2)/b[1]
z = (R - y^2)/b[1]
sol = cbind(x, y, z)
sol

### Test
x^2 + b[1]*y
y^2 + b[1]*z
z^2 + b[1]*x


###############
### Extensions:
R = 1
b = c(1, 1)
#
sol = solve.sysHt32(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^2 + b[1]*y + b[2]*(x+y+z)
y^2 + b[1]*z + b[2]*(x+y+z)
z^2 + b[1]*x + b[2]*(x+y+z)

round0.p(poly.calc(sol[,1]))


### Ex 2:
R = -1
b = c(1, 1)
#
sol = solve.sysHt32(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^2 + b[1]*y + b[2]*(x+y+z)
y^2 + b[1]*z + b[2]*(x+y+z)
z^2 + b[1]*x + b[2]*(x+y+z)

round0.p(poly.calc(sol[,1]))
err = 25 + 12*x^2 - 2*x^3 + 3*x^4 + x^6
round0(err)


### Ex 3:
R = 2
b = c(2, 0, 1)
#
sol = solve.sysHt32(R, b=b)
x = sol[,1]; y = sol[,2]; z = sol[,3]

### Test
x^2 + b[1]*y + b[2]*(x+y+z) + b[3]*(x+y+z)^2
y^2 + b[1]*z + b[2]*(x+y+z) + b[3]*(x+y+z)^2
z^2 + b[1]*x + b[2]*(x+y+z) + b[3]*(x+y+z)^2

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
# - shift all roots "back" & use solution to simple system: x^2 + b1*y = R;
(x + s/2)^2 + b1*(y + s/2) = R + s^2/4 + b1*s/2
(y + s/2)^2 + b1*(z + s/2) = R + s^2/4 + b1*s/2
(z + s/2)^2 + b1*(x + s/2) = R + s^2/4 + b1*s/2

# Trivial solution: x = y = z;

### Diff =>
# (x-y)*(x+y+s) = -b1*(y-z)
# (y-z)*(y+z+s) = b1*(x-z)
# (x-z)*(x+z+s) = b1*(x-y)
### Prod =>
# (x+y+s)*(x+z+s)*(y+z+s) = -b1^3

### Sum =>
# S^2 - 2*E2 + (s+b1)*S = 3*R
# 2*E2 = S^2 + (s+b1)*S - 3*R;

### Sum(x[-i] * ...) =>
# S[x^2*y] + 2*s*E2 + b1*(x^2+y^2+z^2) + b1*E2 = 2*R*S
# S*E2 - 3*E3 + 2*s*E2 + b1*(S^2 - 2*E2) + b1*E2 - 2*R*S = 0
# 3*E3 = S*E2 + 2*s*E2 + b1*(S^2 - 2*E2) + b1*E2 - 2*R*S
# 3*E3 = S*E2 + 2*s*E2 + b1*(S^2 - E2) - 2*R*S
# 3*E3 = b1*S^2 + E2*S + 2*s*E2 - b1*E2 - 2*R*S
# 6*E3 = 2*b1*S^2 + 2*E2*S + 4*s*E2 - 2*b1*E2 - 4*R*S
# 6*E3 = 2*b1*S^2 + S*(S^2 + (s+b1)*S - 3*R) + 4*s*E2 - 2*b1*E2 - 4*R*S
# 6*E3 = S^3 + (3*s + 2*b1)*S^2 + (2*s-b1)*(s+b1)*S - 7*R*S - 6*s*R + 3*b1*R

### Sum(x[i]*...) =>
# x^3 + y^3 + z^3 + s*(x^2 + y^2 + z^2) + b1*E2 = R*S
# S^3 - 3*E2*S + 3*E3 + s*(S^2 - 2*E2) + b1*E2 - R*S = 0
# 2*S^3 - 6*E2*S + 6*E3 + 2*s*(S^2 - 2*E2) + 2*b1*E2 - 2*R*S
# -S^3 - (3*s + 2*b1)*S^2 + 6*E3 + (b1 - 2*s)*(s+b1)*S + 7*R*S + 6*s*R - 3*b1*R
# -S^3 - (3*s + 2*b1)*S^2 + (S^3 + (3*s + 2*b1)*S^2 + (2*s-b1)*(s+b1)*S - 7*R*S - 6*s*R + 3*b1*R) +
#  + (b1 - 2*s)*(s+b1)*S + 7*R*S + 6*s*R - 3*b1*R
# 0 == 0 [redundant]

### Prod =>
# (x+y+s)*(x+z+s)*(y+z+s) = -b1^3
# (x^2 + x*y + x*z + y*z + s*(2*x+y+z) + s^2)*(y+z+s) + b1^3 = 0
# (x^2 + x*y + x*z + y*z)*(y+z+s) + (s*(2*x+y+z) + s^2)*(y+z+s) + b1^3 = 0
# (x^2 + x*y + x*z + y*z)*(y+z) + s*(x^2 + x*y + x*z + y*z) + (s*(2*x+y+z) + s^2)*(y+z+s) + b1^3 = 0
# S*E2 - E3 + s*x^2 + s*E2 + s*(2*x+y+z)*(y+z+s) + s^2*(y+z+s) + b1^3
# S*E2 - E3 + s*x^2 + s*E2 + s*(2*x+y+z)*(y+z) + s^2*(2*x+y+z) + s^2*(y+z+s) + b1^3
# S*E2 - E3 + s*x^2 + s*E2 + s*(2*x+y+z)*(y+z) + 2*s^2*(x+y+z) + s^3 + b1^3
# S*E2 - E3 + s*(x^2+y^2+z^2) + s*E2 + 2*s*E2 + 2*s^2*S + s^3 + b1^3
# S*E2 - E3 + s*S^2 + s*E2 + 2*s^2*S + s^3 + b1^3
# 2*S*E2 - 2*E3 + 2*s*S^2 + 2*s*E2 + 4*s^2*S + 2*s^3 + 2*b1^3
# S*(S^2 + (s+b1)*S - 3*R) - 2*E3 + 2*s*S^2 + s*(S^2 + (s+b1)*S - 3*R) + 4*s^2*S + 2*s^3 + 2*b1^3
# S^3 + (4*s+b1)*S^2 - 2*E3 - 3*R*S + 5*s^2*S + s*b1*S + 2*s^3 + 2*b1^3 - 3*s*R
# 3*S^3 + 3*(4*s+b1)*S^2 - 6*E3 - 9*R*S + 15*s^2*S + 3*s*b1*S + 6*s^3 + 6*b1^3 - 9*s*R
# 2*S^3 + (9*s+b1)*S^2 + b1^2*S - 2*R*S + 13*s^2*S + 2*s*b1*S + 6*s^3 + 6*b1^3 - 3*s*R - 3*b1*R
### Eq:
(2*S + 3*s + 3*b1)*(S^2 - (b1 - 3*s)*S + (- R - 2*b1*s + 2*b1^2 + 2*s^2))

# Note:
# Bug corrected: + 6*b1^3;

### Solver:

solve.htShX.S3P2 = function(R, b, s) {
	coeff = c(1, - (b[1] - 3*s), (- R[1] - 2*b[1]*s + 2*b[1]^2 + 2*s^2))
	x.sum = roots(coeff)
	E3 = (x.sum^3 + (3*s + 2*b[1])*x.sum^2 +
		(2*s-b[1])*(s+b[1])*x.sum - 7*R*x.sum - 6*s*R + 3*b[1]*R) / 6;
	E2 = (x.sum^2 + (s+b[1])*x.sum - 3*R)/2
	x = as.vector(sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id]))))
	y = (R - x^2 - s*x)/b[1]
	z = (R - y^2 - s*y)/b[1]
	sol = cbind(x, y, z)
	sol = rbind(sol, sol[,c(2,3,1)], sol[,c(3,1,2)])
	return(sol);
}

### Example
R = 1
b = 3
s = 1
#
sol = solve.htShX.S3P2(R, b, s)
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

# Trivial solution: x = y = z;

### Solution

### Diff =>
# x^2 - y^2 = b1*(x-y)
# x^2 - z^2 = b1*(x-z)
# y^2 - z^2 = b1*(y-z)

# if x != y != z
# x + y = b1
# x + z = b1
# y + z = b1
# => x = y = z, which violates assumption;

# =>
# x = y
# x^2 + b1*(x+z) = R
# z^2 + 2*b1*x = R

# Case x != z
# x + z = b1
# x^2 + b1^2 = R

### Example:
b = 3
R = 1
#
x = sqrt(R - b[1]^2 + 0i)
x = c(x, -x)
y = x
z = b[1] - x
sol = cbind(x, y, z)
sol

### Test
x^2 + b[1]*(y+z)
y^2 + b[1]*(x+z)
z^2 + b[1]*(x+y)

####################

### Shifted Variant:
### (x[i] - s)^2 + b*(x[j] + x[k])

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
b = 3
s = -1
R = 1
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

### x[i]^2 + b*x[-i]

# x^2 + b1*y*z = R
# y^2 + b1*x*z = R
# z^2 + b1*x*y = R

### Solution

# - degenerate P4;

### Special Case: b1 = 2
# Z^2 = 3*R

### Diff =>
# x^2 - y^2 = b1*z*(x-y)
# x^2 - z^2 = b1*y*(x-z)
# y^2 - z^2 = b1*x*(y-z)
# => if x != y != z
#   x + y - b1*z = 0
#   x - b1*y + z = 0
# - b1*x + y + z = 0
# => x = y = z = 0; # Contradiction !!!
# x == y & x != z =>
# x + z = b1*y;
# z = (b1 - 1)*x;

### Case: x = y
# x^2 + b1*x*z = R
# z^2 + b1*x^2 = R
# =>
# b1*z = R/x - x
# b1^2*z^2 + b1^3*x^2 = b1^2*R
# x^2 - 2*R + R^2/x^2 + b1^3*x^2 - b1^2*R = 0
### Eq:
(b1^3+1)*x^4 - R*(b1^2 + 2)*x^2 + R^2 # = 0
((b1^2 - b1 + 1)*x^2 - R) * ((b1+1)*x^2 - R)


### Solver:

solve.htYZ.S3P2 = function(R, b, b.ext=0) {
	coeff1 = c((b[1]^2 - b[1] + 1), 0, -R[1]); # x == y;
	coeff2 = c((b[1] + 1), 0, -R[1]); # x == y == z;
	# computes x, NOT S;
	# simple extensions: adapted to x;
	if(b.ext[1] != 0) {
		coeff1 = coeff1 + c(0, b.ext[1] * (b[1]+1), 0); # x == y
		coeff2 = coeff2 + c(0, 3*b.ext[1], 0);
	}
	sol1 = roots(coeff1)
	sol2 = roots(coeff2)
	x = c(sol1, sol2);
	y = x
	# z = (b[1] - 1)*x; # assumes x == y & x != z;
	z = (R[1] - x^2 - 2*x*b.ext[1]) / (b[1]*x + b.ext[1]);
	sol = round0(cbind(x=x, y=y, z=z))
}

### Example:
R = 2
b = 1
#
sol = solve.htYZ.S3P2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + b[1]*y*z
y^2 + b[1]*x*z
z^2 + b[1]*x*y


### Extension A1:
R = 2
b = 1
b.ext = c(1)
#
sol = solve.htYZ.S3P2(R, b, b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + b[1]*y*z + b.ext[1]*(x+y+z)
y^2 + b[1]*x*z + b.ext[1]*(x+y+z)
z^2 + b[1]*x*y + b.ext[1]*(x+y+z)


######################
######################

### Prod-Type: Order 3
### x[i]^3 + b*x[-i]

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

solve.htyz.S3P3 = function(R, b, do2Eq = FALSE, tol=1E-8) {
	b1 = b[1]; R = R[1];
	coeff = c(1, - 3*b1, 7*b1^2, -5*b1^3, - 27*R*b1, 18*R*b1^2, (5*R*b1^3 + 27*R^2))
	S = roots(coeff)
	# Case: x != y != z
	S = c(S, -b1);
	#
	E2.part = 27*R*S*b1 - 72*R*S^2 + 45*R*b1^2 - 15*S^3*b1^2 + S^4*b1 + 12*S^5;
	divE2 = (54*R + 31*S*b1^2 + 12*S^2*b1 - 30*S^3 - 15*b1^3);
	E2 = - E2.part / divE2
	E3 = - (S^3 - 3*E2*S + b1*E2 - 3*R) / 3
	# robust complex solutions for x;
	E2 = round0(E2, tol=tol); E3 = round0(E3, tol=tol);
	#
	x = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	# yz = E3 / x;
	yz = (R - x^3)/b1; # robust ???
	yz.s = S - x;
	#
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d)/2;
	z = (yz.s - yz.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	# all permutations: x != y != z
	sol2 = sol[nrow(sol) - 0:2, c(1,3,2)]
	sol = rbind(sol, sol2)
	### y == z: ??? more robust ???
	if(do2Eq) {
		isEq = round0((yz.s)^2 - 4*yz, tol=1E-7) == 0
		sol = sol[ ! as.vector(isEq) , ]
		y = z = as.vector(yz.s [isEq]) / 2;
		sol2 = cbind(x=as.vector(x[isEq]), y=y, z=z)
		sol = rbind(sol, sol2)
	}
	return(sol)
}

### Examples:

R = 1;
b = 2;
#
sol = solve.htyz.S3P3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b*y*z # - R
y^3 + b*x*z # - R
z^3 + b*x*y # - R


### Debug
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
# (x^3+y^3+z^3) + (b1 + b2)*E2 = R*S
# S^3 - 3*E2*S + 3*E3 - R*S + (b1 + b2)*E2 = 0
# 2*S^3 - 2*3*E2*S + 6*E3 - 2*R*S + 2*(b1 + b2)*E2 = 0
# 2*S^3 - 3*(S^2 + (b1+b2)*S - 3*R)*S + 6*E3 - 2*R*S + (b1 + b2)*(S^2 + (b1+b2)*S - 3*R) = 0
# -S^3 - 2*(b1+b2)*S^2 + 6*E3 + 7*R*S + (b1+b2)*((b1+b2)*S - 3*R) = 0
# 6*E3 = S^3 + 2*(b1+b2)*S^2 - (b1+b2)^2*S - 7*R*S + 3*(b1+b2)*R

### Eq3:
# x^2 + b1*y + b1*z = R - (b2-b1)*z
# x^4 + b1^2*(y+z)^2 + 2*b1*x^2*(y+z) = R^2 + (b2-b1)^2*z^2 - 2*(b2-b1)*R*z
### Sum =>
# (x^4+y^4+z^4) + 2*b1^2*(x^2+y^2+z^2 + E2) + 2*b1*(E2*S - 3*E3) = 3*R^2 + (b2-b1)^2*(S^2 - 2*E2) - 2*(b2-b1)*R*S
# (x^4+y^4+z^4) + 2*b1^2*(S^2 - E2) + 2*b1*E2*S - 6*b1*E3 = 3*R^2 + (b2-b1)^2*S^2 - 2*(b2-b1)^2*E2 - 2*(b2-b1)*R*S
# (x^4+y^4+z^4) = 3*R^2 + (-b1^2 + b2^2 - 2*b1*b2)*S^2 - 2*b1*E2*S - 2*(b2^2 - 2*b1*b2)*E2 + 6*b1*E3 - 2*(b2-b1)*R*S
# S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 =
#  = b1*b2*S^2 - (b1^2 + b2^2 - b1*b2)*(b1+b2)*S - 2*b2*R*S - 2*b1*R*S + 3*R^2 + 3*(b1^2 + b2^2 - b1*b2)*R
# - 2*S^4 - 4*(b1+b2)*S^3 + 8*R*S^2 - 4*(b1^2 + b2^2)*S^2 - 14*b1*b2*S^2 + 3*(S^2 + (b1+b2)*S - 3*R)^2 + 24*(b1+b2)*R*S =
#  = - 6*(b1^2 + b2^2 - b1*b2)*(b1+b2)*S + 18*R^2 + 18*(b1^2 + b2^2 - b1*b2)*R
# S^4 + 2*(b1+b2)*S^3 - 10*R*S^2 - (b1^2 + b2^2 + 8*b1*b2)*S^2 + 6*(b1+b2)*R*S +
#  + 6*(b1^2 + b2^2 - b1*b2)*(b1+b2)*S + 9*R^2 - 18*(b1^2 + b2^2 - b1*b2)*R  =  0;


### Diff =>
# [NOT used]
# x^2 - y^2 = b2*(x - z) - b1*(y - z)
# x^2 - z^2 = b2*(y - z) - b1*(y - x)
# y^2 - z^2 = b2*(y - x) - b1*(z - x)


solve.htS3L2 = function(b, R) {
	b.sum = b[1] + b[2]
	coeff = c(1, 2*b.sum, - 10*R - (b.sum^2 + 6*b[1]*b[2]), 6*b.sum*R + 6*(b.sum^2 - 3*b[1]*b[2])*b.sum,
		9*R^2 - 18*(b.sum^2 - 3*b[1]*b[2])*R)
	x.sum = roots(coeff)
	E2 = (x.sum^2 + b.sum*x.sum - 3*R)/2
	E3 = - (x.sum^3 - 3*E2*x.sum - R*x.sum + b.sum*E2)/3
	# E3 = (x.sum^3 + 2*b.sum*x.sum^2 - b.sum^2*x.sum - 7*R*x.sum + 3*b.sum*R) / 6
	# solve (x, y, z)
	x = as.vector(sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id]))) )
	x.sum = rep(x.sum, each=3)
	yz.sum = x.sum - x
	yz.b.sum = R - x^2
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

### Example:
b = c(1,-2)
R = 1
#
sol = solve.htS3L2(b, R)
x = sol[,1]; y = sol[,2]; z = sol[,3]
sol
# TODO:
# - find and correct bugs;

### Test
x^2 + b[1]*y + b[2]*z
y^2 + b[1]*z + b[2]*x
z^2 + b[1]*x + b[2]*y

### Classic Polynomial
# Note: set of valid roots may change!
round0.p(poly.calc(x[c(4:9)]))
round0.p(poly.calc(x[-c(4:9)]))
#
171 + 36*x - 15*x^2 + 3*x^3 - 2*x^4 + x^5 + x^6
-1 - 3*x + 5*x^3 - 3*x^5 + x^6

### Derivation:

# b2*z = -x^2 - b1*y + R
# =>
# b2*y^2 - b1^2*y - b1*x^2 + b2^2*x + (b1 - b2)*R # = 0
# &
# (x^2 + b1*y - R)^2 + b1*b2^2*x + b2^3*y - b2^2*R # = 0
# x^4 + b1^2*y^2 + R^2 + 2*b1*x^2*y - 2*R*x^2 - 2*b1*R*y + b1*b2^2*x + b2^3*y - b2^2*R # = 0
# 4*b2^2*x^4 - 8*b2^2*R*x^2 + b1^2*(4*b2^2*y^2) + (4*b1*b2*x^2 + 2*b2^4 - 4*b1*b2*R)*2*b2*y + 4*b1*b2^4*x - 4*b2^4*R + 4*b2^2*R^2 # = 0
# Sol Eq 2 =>
# y = (b1^2 + sqrt(b1^4 + 4*b1*b2*x^2 - 4*b2^3*x - 4*b2*(b1 - b2)*R)) / (2*b2)
# y = (b1^2 - sqrt(b1^4 + 4*b1*b2*x^2 - 4*b2^3*x - 4*b2*(b1 - b2)*R)) / (2*b2)
# =>
4*b2^2*x^4 - 8*b2^2*R*x^2 +
	+ 4*b1^3*b2*x^2 - 4*b1^2*b2^3*x - 4*b1^3*b2*R + 4*b1^2*b2^2*R +
	+ 4*b1^3*b2*x^2 + 2*b1^2*b2^4 - 4*b1^3*b2*R + 4*b1*b2^4*x - 4*b2^4*R + 4*b2^2*R^2 + 2*b1^6 + # b1^4 +
	(2*b1^4 + 2*b2^4 + 4*b1*b2*x^2 - 4*b1*b2*R)*sqrt(b1^4 + 4*b1*b2*x^2 - 4*b2^3*x - 4*b2*(b1 - b2)*R)
###
(16*R*b1*b2^9 - 16*R*b1^2*b2^8 - 16*R*b1^4*b2^6 + 16*R*b1^5*b2^5 - 16*R*b1^6*b2^4 - 16*R*b2^10 + 64*R^2*b1*b2^7 - 80*R^2*b1^2*b2^6 + 64*R^2*b1^3*b2^5 + 16*R^2*b1^4*b2^4 + 16*R^2*b2^8 - 32*R^3*b1^2*b2^4 - 32*R^3*b2^6 + 16*R^4*b2^4) +
(- 96*R*b1*b2^8 + 32*R*b1^2*b2^7 + 32*R*b1^3*b2^6 - 96*R*b1^4*b2^5 + 32*R^2*b1*b2^6 + 32*R^2*b1^2*b2^5 + 16*b1^3*b2^8 + 16*b1^4*b2^7 + 16*b1^7*b2^4 + 16*b2^11)*x^1 +
(- 64*R*b1*b2^7 + 96*R*b1^2*b2^6 - 64*R*b1^3*b2^5 + 64*R^2*b1^2*b2^4 + 64*R^2*b2^6 - 64*R^3*b2^4 - 16*b1*b2^9 + 16*b1^2*b2^8 - 32*b1^3*b2^7 + 16*b1^4*b2^6 - 16*b1^5*b2^5)*x^2 +
(- 64*R*b1*b2^6 - 64*R*b1^2*b2^5 + 64*b1*b2^8 + 64*b1^4*b2^5)*x^3 +
(- 32*R*b1^2*b2^4 - 32*R*b2^6 + 96*R^2*b2^4 - 48*b1^2*b2^6)*x^4 +
(32*b1*b2^6 + 32*b1^2*b2^5)*x^5 +
(- 64*R*b2^4)*x^6 +
(16*b2^4)*x^8

b1 = b[1]; b2 = b[2];
coeff = c(
	(16*R*b1*b2^9 - 16*R*b1^2*b2^8 - 16*R*b1^4*b2^6 + 16*R*b1^5*b2^5 - 16*R*b1^6*b2^4 - 16*R*b2^10 + 64*R^2*b1*b2^7 - 80*R^2*b1^2*b2^6 +
	 64*R^2*b1^3*b2^5 + 16*R^2*b1^4*b2^4 + 16*R^2*b2^8 - 32*R^3*b1^2*b2^4 - 32*R^3*b2^6 + 16*R^4*b2^4),
	(- 96*R*b1*b2^8 + 32*R*b1^2*b2^7 + 32*R*b1^3*b2^6 - 96*R*b1^4*b2^5 + 32*R^2*b1*b2^6 + 32*R^2*b1^2*b2^5 + 16*b1^3*b2^8 + 16*b1^4*b2^7 + 16*b1^7*b2^4 + 16*b2^11),
	(- 64*R*b1*b2^7 + 96*R*b1^2*b2^6 - 64*R*b1^3*b2^5 + 64*R^2*b1^2*b2^4 + 64*R^2*b2^6 - 64*R^3*b2^4 - 16*b1*b2^9 + 16*b1^2*b2^8 - 32*b1^3*b2^7 + 16*b1^4*b2^6 - 16*b1^5*b2^5),
	(- 64*R*b1*b2^6 - 64*R*b1^2*b2^5 + 64*b1*b2^8 + 64*b1^4*b2^5),
	(- 32*R*b1^2*b2^4 - 32*R*b2^6 + 96*R^2*b2^4 - 48*b1^2*b2^6),
	(32*b1*b2^6 + 32*b1^2*b2^5),
	(- 64*R*b2^4),
	0, (16*b2^4)  ### 0*x^7, x^8
)
coeff = rev(coeff)
coeff

### P[8] is decomposable:
# P[8] = P[2]*P[6]
# P[2] = x^2 + (b1+b2)*x - R;
### P[6]
(- 64*R*b1*b2^7 + 80*R*b1^2*b2^6 - 64*R*b1^3*b2^5 - 16*R*b1^4*b2^4 - 16*R*b2^8 + 32*R^2*b1^2*b2^4 + 32*R^2*b2^6 - 16*R^3*b2^4 +
 - 16*b1*b2^9 + 16*b1^2*b2^8 + 16*b1^4*b2^6 - 16*b1^5*b2^5 + 16*b1^6*b2^4 + 16*b2^10) +
(32*R*b1^3*b2^4 + 32*R*b2^7 - 16*R^2*b1*b2^4 - 16*R^2*b2^5 + 16*b1*b2^8 - 16*b1^2*b2^7 - 16*b1^3*b2^6 + 16*b1^4*b2^5 - 16*b1^5*b2^4 - 16*b2^9)*x^1 +
(- 32*R*b1*b2^5 - 48*R*b1^2*b2^4 - 48*R*b2^6 + 48*R^2*b2^4 + 32*b1*b2^7 - 16*b1^2*b2^6 + 32*b1^3*b2^5 + 16*b1^4*b2^4 + 16*b2^8)*x^2 +
(32*R*b1*b2^4 + 32*R*b2^5 - 16*b1*b2^6 - 16*b1^2*b2^5 - 16*b1^3*b2^4 - 16*b2^7)*x^3 +
(- 48*R*b2^4 + 32*b1*b2^5 + 16*b1^2*b2^4 + 16*b2^6)*x^4 +
(- 16*b1*b2^4 - 16*b2^5)*x^5 +
(16*b2^4)*x^6

b1 = b[1]; b2 = b[2]
coeff = c(
	(- 64*R*b1*b2^7 + 80*R*b1^2*b2^6 - 64*R*b1^3*b2^5 - 16*R*b1^4*b2^4 - 16*R*b2^8 + 32*R^2*b1^2*b2^4 + 32*R^2*b2^6 - 16*R^3*b2^4 +
	- 16*b1*b2^9 + 16*b1^2*b2^8 + 16*b1^4*b2^6 - 16*b1^5*b2^5 + 16*b1^6*b2^4 + 16*b2^10),
	(32*R*b1^3*b2^4 + 32*R*b2^7 - 16*R^2*b1*b2^4 - 16*R^2*b2^5 + 16*b1*b2^8 - 16*b1^2*b2^7 - 16*b1^3*b2^6 + 16*b1^4*b2^5 - 16*b1^5*b2^4 - 16*b2^9),
	(- 32*R*b1*b2^5 - 48*R*b1^2*b2^4 - 48*R*b2^6 + 48*R^2*b2^4 + 32*b1*b2^7 - 16*b1^2*b2^6 + 32*b1^3*b2^5 + 16*b1^4*b2^4 + 16*b2^8),
	(32*R*b1*b2^4 + 32*R*b2^5 - 16*b1*b2^6 - 16*b1^2*b2^5 - 16*b1^3*b2^4 - 16*b2^7),
	(- 48*R*b2^4 + 32*b1*b2^5 + 16*b1^2*b2^4 + 16*b2^6),
	(- 16*b1*b2^4 - 16*b2^5), (16*b2^4)
)
coeff = rev(coeff)
coeff

x = roots(coeff)
yz.bsum = R - x^2
y = as.vector(sapply(1:length(x), function(id) roots(c(1, -b[1]^2/b[2], -R + b[2]*x[id] + b[1]/b[2]*yz.bsum[id]))))
yz.bsum = rep(yz.bsum, each=2); x = rep(x, each=2);
z = (yz.bsum - b[1]*y) / b[2]

### Test
x^2 + b[1]*y + b[2]*z
y^2 + b[1]*z + b[2]*x
z^2 + b[1]*x + b[2]*y

# sapply(x, function(x) sum(x^(8:0) * coeff))
# poly.calc(x[c(7,8)])
# poly.calc(x[-c(7,8)])


##############
### Extensions

# x^2 + s*x + b2*y + b1*z = R;
# s = shift of main variable;
# x^2 + s*x + b3*x*y*z + b2*y + b1*z = R;

##############
### Extension:
# x^2 + s*x + b3*x*y*z + b2*y + b1*z = R

### Sum =>
# S^2 - 2*E2 + 3*b3*E3 + (s+b1+b2)*S = 3*R
# E2 = (S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R)/2;

### Sum(x[i]*P[i]) =>
# (x^3+y^3+z^3) + s*(S^2 - 2*E2) + b3*E3*S + (b1 + b2)*E2 = R*S
# S^3 + s*S^2 - 3*E2*S + 3*E3 + b3*E3*S - 2*s*E2 - R*S + (b1 + b2)*E2 = 0
# 2*S^3 + 2*s*S^2 - 2*3*E2*S + 6*E3 + 2*b3*E3*S - 4*s*E2 - 2*R*S + 2*(b1 + b2)*E2 = 0
# 2*S^3 + 2*s*S^2 - 3*S*(S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R) + 6*E3 + 2*b3*E3*S +
#  - 2*s*(S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R) - 2*R*S +
#  + (b1 + b2)*(S^2 + (s+b1+b2)*S + 3*b3*E3 - 3*R) = 0
# -S^3 - 3*s*S^2 - 2*(b1+b2)*S^2 - 2*s^2*S - s*(b1+b2)*S + (b1+b2)^2*S + 7*R*S + 6*s*R - 3*(b1+b2)*R =
#  = (7*b3*S - 3*b3*(b1 + b2) + 6*s*b3 - 6)*E3

### Eq3:
# x^2 + s*x + b3*x*y*z + b1*y + b1*z = R - (b2-b1)*z
# (x^2 + s*x + b3*x*y*z + b1*y + b1*z)^2 = (R + b1*z - b2*z)^2

2*E1_2*E3*b3 + 2*E1_2*b1*b2 + E1_2*b1^2 - E1_2*b2^2 + E1_2*s^2 + 2*E1_3*s + E1_4 + 4*E2*b1*s + 2*E2*b1^2 + 2*E2c1*b1 + 4*E3*S*b1*b3 + 2*E3*S*b3*s + 3*E3^2*b3^2 - 2*R*S*b1 + 2*R*S*b2 - 3*R^2

- 4*E2*E3*b3 + 2*E2*S*b1 - 6*E2*S*s - 4*E2*S^2 - 4*E2*b1*b2 + 4*E2*b1*s + 2*E2*b2^2 - 2*E2*s^2 + 2*E2^2 + 4*E3*S + 4*E3*S*b1*b3 + 2*E3*S*b3*s + 2*E3*S^2*b3 - 6*E3*b1 + 6*E3*s + 3*E3^2*b3^2 - 2*R*S*b1 + 2*R*S*b2 - 3*R^2 + 2*S^2*b1*b2 + S^2*b1^2 - S^2*b2^2 + S^2*s^2 + 2*S^3*s + S^4

- 6*E3*R*b3 + 8*E3*S + 16*E3*S*b1*b3 + 2*E3*S*b2*b3 - 12*E3*S*b3*s - 6*E3*S^2*b3 - 12*E3*b1 - 12*E3*b1*b2*b3 + 12*E3*b1*b3*s + 6*E3*b2^2*b3 - 6*E3*b3*s^2 + 12*E3*s + 3*E3^2*b3^2 - 16*R*S*b1 - 2*R*S*b2 + 12*R*S*s + 6*R*S^2 + 12*R*b1*b2 - 12*R*b1*s - 6*R*b2^2 + 6*R*s^2 + 3*R^2 - 2*S*b1*b2^2 + 2*S*b1*s^2 - 4*S*b1^2*b2 + 4*S*b1^2*s - 2*S*b2*s^2 + 2*S*b2^2*s + 2*S*b2^3 - 2*S*s^3 + 4*S^2*b1*b2 + 2*S^2*b1*s + 5*S^2*b1^2 - 4*S^2*b2*s + S^2*b2^2 - 5*S^2*s^2 - 2*S^3*b2 - 4*S^3*s - S^4

(216*R*b1*b2 - 324*R*b1*b3*s^2 + 216*R*b1*s - 216*R*b1^2 + 324*R*b1^2*b3*s - 108*R*b1^3*b3 - 324*R*b2*b3*s^2 + 216*R*b2*s - 216*R*b2^2 + 324*R*b2^2*b3*s - 108*R*b2^3*b3 + 216*R*b3*s^3 - 216*R*s^2 + 108*R^2) +
(72*R*b1 - 108*R*b1*b2*b3 - 540*R*b1*b3*s + 324*R*b1^2*b3 + 72*R*b2 - 540*R*b2*b3*s + 324*R*b2^2*b3 + 540*R*b3*s^2 - 360*R*s + 216*b1*b2*b3*s^2 - 216*b1*b2*s - 108*b1*b2^2*b3*s + 36*b1*b2^3*b3 + 36*b1*b3*s^3 - 108*b1^2*b2*b3*s + 72*b1^3 + 36*b1^3*b2*b3 - 72*b1^3*b3*s + 36*b1^4*b3 + 36*b2*b3*s^3 + 72*b2^3 - 72*b2^3*b3*s + 36*b2^4*b3 - 72*b3*s^4 + 72*s^3)*S^1 +
(- 120*R - 216*R*b1*b3 - 216*R*b2*b3 + 432*R*b3*s - 96*b1*b2 + 396*b1*b2*b3*s - 72*b1*b2^2*b3 + 36*b1*b3*s^2 + 48*b1*s - 12*b1^2 - 72*b1^2*b2*b3 + 36*b1^2*b3*s - 96*b1^3*b3 + 36*b2*b3*s^2 + 48*b2*s - 12*b2^2 + 36*b2^2*b3*s - 96*b2^3*b3 - 204*b3*s^3 + 132*s^2)*S^2 +
(104*R*b3 + 24*b1 + 148*b1*b2*b3 - 20*b1*b3*s + 12*b1*b3^2*s^2 + 44*b1^2*b3 - 12*b1^2*b3^2*s + 4*b1^3*b3^2 + 24*b2 - 20*b2*b3*s + 12*b2*b3^2*s^2 + 44*b2^2*b3 - 12*b2^2*b3^2*s + 4*b2^3*b3^2 - 196*b3*s^2 - 8*b3^2*s^3 + 72*s)*S^3 +
(12 + 4*b1*b2*b3^2 - 16*b1*b3 + 20*b1*b3^2*s - 12*b1^2*b3^2 - 16*b2*b3 + 20*b2*b3^2*s - 12*b2^2*b3^2 - 72*b3*s - 20*b3^2*s^2)*S^4 +
(8*b1*b3^2 + 8*b2*b3^2 - 8*b3 - 16*b3^2*s)*S^5 +
(- 4*b3^2)*S^6

###############

solve.htS3X2Ext = function(b, R, s) {
	b1 = b[1]; b2 = b[2]; b3 = b[3]
	coeff = c(
	(216*R*b1*b2 - 324*R*b1*b3*s^2 + 216*R*b1*s - 216*R*b1^2 + 324*R*b1^2*b3*s - 108*R*b1^3*b3 - 324*R*b2*b3*s^2 +
	216*R*b2*s - 216*R*b2^2 + 324*R*b2^2*b3*s - 108*R*b2^3*b3 + 216*R*b3*s^3 - 216*R*s^2 + 108*R^2),
	(72*R*b1 - 108*R*b1*b2*b3 - 540*R*b1*b3*s + 324*R*b1^2*b3 + 72*R*b2 - 540*R*b2*b3*s + 324*R*b2^2*b3 +
	540*R*b3*s^2 - 360*R*s + 216*b1*b2*b3*s^2 - 216*b1*b2*s - 108*b1*b2^2*b3*s + 36*b1*b2^3*b3 +
	36*b1*b3*s^3 - 108*b1^2*b2*b3*s + 72*b1^3 + 36*b1^3*b2*b3 - 72*b1^3*b3*s + 36*b1^4*b3 + 36*b2*b3*s^3 +
	72*b2^3 - 72*b2^3*b3*s + 36*b2^4*b3 - 72*b3*s^4 + 72*s^3),
	(- 120*R - 216*R*b1*b3 - 216*R*b2*b3 + 432*R*b3*s - 96*b1*b2 + 396*b1*b2*b3*s - 72*b1*b2^2*b3 +
	36*b1*b3*s^2 + 48*b1*s - 12*b1^2 - 72*b1^2*b2*b3 + 36*b1^2*b3*s - 96*b1^3*b3 + 36*b2*b3*s^2 +
	48*b2*s - 12*b2^2 + 36*b2^2*b3*s - 96*b2^3*b3 - 204*b3*s^3 + 132*s^2),
	(104*R*b3 + 24*b1 + 148*b1*b2*b3 - 20*b1*b3*s + 12*b1*b3^2*s^2 + 44*b1^2*b3 - 12*b1^2*b3^2*s +
	4*b1^3*b3^2 + 24*b2 - 20*b2*b3*s + 12*b2*b3^2*s^2 + 44*b2^2*b3 - 12*b2^2*b3^2*s + 4*b2^3*b3^2 - 196*b3*s^2 - 8*b3^2*s^3 + 72*s),
	(12 + 4*b1*b2*b3^2 - 16*b1*b3 + 20*b1*b3^2*s - 12*b1^2*b3^2 - 16*b2*b3 + 20*b2*b3^2*s - 12*b2^2*b3^2 - 72*b3*s - 20*b3^2*s^2),
	(8*b1*b3^2 + 8*b2*b3^2 - 8*b3 - 16*b3^2*s), (- 4*b3^2)
	)
	coeff = rev(coeff)
	print(coeff)
	x.sum = round0(roots(coeff))
	x.sum = x.sum[x.sum != 0]
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

### Example
s = 1
b = c(-2,2,3)
R = 1
#
sol = solve.htS3X2Ext(b, R, s)
x = sol[,1]; y = sol[,2]; z = sol[,3];
sol

### Test
# TODO: precision error vs bug ???
# some values are correct
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

### Solution

# Trivial solution: x = y = z

# Formulas used:
# x^3 + y^3 + z^3 = S^3 - 3*E2*S + 3*E3;
# x^4 + y^4 + z^4 = S^4 - 4*E2 * S^2 + 4*E3 * S + 2*E2^2;
# x^3*y^3 + x^3*z^3 + y^3*z^3 = E2^3 + 3*E3^2 - 3*E3*E2*S
# x^4*y^2 + x^4*z^2+ x^2*y^4 + y^4*z^2 + y^2*z^4 + x^2*z^4 = E2^2*S^2 - 2*E3*S^3 + 4*E3*E2*S - 3*E3^2 - 2*E2^3

### Sum =>
# x^3 + y^3 + z^3 + b*(x+y+z) = 3*R
S^3 - 3*E2*S + 3*E3 + b*S - 3*R # = 0;
# 3*E2*S = S^3 + 3*E3 + b*S - 3*R
# =>
# 6*E3 = 6*E2*S - 2*S^3 - 2*b*S + 6*R

### Sum(x[i]*...) =>
# x^4 + y^4 + z^4 + b*E2 = R*S
S^4 - 4*E2*S^2 + 4*E3 * S + 2*E2^2 + b*E2 - R*S # = 0;

### Solve E2 & E3:
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 + b*E2 - R*S # = 0
3*S^4 - 12*E2*S^2 + 2*6*E3*S + 6*E2^2 + 3*b*E2 - 3*R*S # = 0
3*S^4 - 12*E2*S^2 + 2*(6*E2*S - 2*S^3 - 2*b*S + 6*R)*S + 6*E2^2 + 3*b*E2 - 3*R*S # = 0
# E2 =>
6*E2^2 + 3*b*E2 - (S^4 + 4*b*S^2 - 9*R*S) # = 0
# alternative: E2 =>
3*S^4 - 4*3*E2*S^2 + 12*E3*S + 2*3*E2^2 + 3*b*E2 - 3*R*S # = 0
3*S^4 - 4*(S^3 + 3*E3 + b*S - 3*R)*S + 12*E3*S + 2*3*E2^2 + 3*b*E2 - 3*R*S # = 0
S^4 + 4*b*S^2 - 9*R*S - 2*3*E2^2 - 3*b*E2 # = 0
S^5 + 3*b*S^3 - 9*R*S^2 - b^2*S - 2*(S^3 + 3*E3 + b*S - 3*R)*E2 - 3*b*E3 + 3*b*R # = 0
S^5 + 4*b*S^3 - 9*R*S^2 - 2*(S^3 + 3*E3 + 5/2*b*S - 3*R)*E2 # = 0
# ??? robust ???

### Sum(x[i+1]^3*...) =>
(x^3*y^3 + x^3*z^3 + y^3*z^3) + b*(x^4 + y^4 + z^4) - R*(x^3 + y^3 + z^3) # = 0
E2^3 + 3*E3^2 - 3*E3*E2*S + b*(S^4 - 4*E2 * S^2 + 4*E3 * S + 2*E2^2) +
	- R*(S^3 - 3*E2*S + 3*E3)
b*S^4 - R*S^3 - 4*b*E2*S^2 + 3*R*E2*S + 4*b*E3*S - 3*E3*E2*S - 3*R*E3 + E2^3 + 2*b*E2^2 + 3*E3^2

### [old] Diff =>
# x^3 - y^3 = b*(z-y)
# y^3 - z^3 = b*(x-z)
# z^3 - x^3 = b*(y-x)
### Prod =>
# (x^2 + y^2 + x*y)*(x^2 + z^2 + x*z)*(y^2 + z^2 + y*z) = -b^3
# =>
# - E2^3 + E2^2*S^2 - E3*S^3 + b^3 = 0
# E2^3 - E2^2*S^2 + E3*S^3 - b^3 = 0

### Solve: E2, E3, S
6*E2^2 + 3*b*E2 - (S^4 + 4*b*S^2 - 9*R*S)
6*E3 - (6*E2*S - 2*S^3 - 2*b*S + 6*R)
b*S^4 - R*S^3 - 4*b*E2*S^2 + 3*R*E2*S + 4*b*E3*S - 3*E3*E2*S - 3*R*E3 + E2^3 + 2*b*E2^2 + 3*E3^2
# [old]
E2^3 - E2^2*S^2 + E3*S^3 - b^3 # = 0

### E2, E3 =>
12*b*S^4 - 12*R*S^3 - 48*b*E2*S^2 + 36*R*E2*S + 8*b*(6*E2*S - 2*S^3 - 2*b*S + 6*R)*S +
	- 6*(6*E2*S - 2*S^3 - 2*b*S + 6*R)*E2*S - 6*R*(6*E2*S - 2*S^3 - 2*b*S + 6*R) +
	+ 12*E2^3 + 24*b*E2^2 + (6*E2*S - 2*S^3 - 2*b*S + 6*R)^2
4*S^6 - 12*E2*S^4 + 4*b*S^4 - 24*R*S^3 - 12*b^2*S^2 - 12*b*E2*S^2 + 36*b*R*S +
	+ 36*R*E2*S + 12*E2^3 + 24*b*E2^2  # (6*E2*S - 2*S^3 - 2*b*S + 6*R)^2
4*S^6 + 7*b*S^4 - 24*R*S^3 + 9*b*R*S - 10*E2*S^4 - 4*b*E2*S^2 + 18*R*E2*S - 9*b^2*E2
# 4*E2 = - b + 1/3 * sqrt(Det); 4*E2 = - b - 1/3 * sqrt(Det);
Det = (9*b^2 + 24*S^4 + 96*b*S^2 - 216*R*S)
# =>
4*S^6 + 7*b*S^4 - 24*R*S^3 + 9*b*R*S - (10*S^4 + 4*b*S^2 - 18*R*S + 9*b^2)*E2
16*S^6 + 28*b*S^4 - 96*R*S^3 + 36*b*R*S - (10*S^4 + 4*b*S^2 - 18*R*S + 9*b^2)*(- b + 1/3 * sqrt(Det))
48*S^6 + 114*b*S^4 + 12*b^2*S^2 - 288*R*S^3 + 54*b*R*S + 27*b^3 - (10*S^4 + 4*b*S^2 - 18*R*S + 9*b^2)*sqrt(Det)

###############
### "P12" / P11 / P8 Polynomial:
S^12 + 6*b*S^10 - 27*R*S^9 - 9*b^2*S^8 + 54*R*b*S^7 + (27*R^2 + 166*b^3)*S^6 +
	+ (- 756*R*b^2)*S^5 + (972*R^2*b + 45*b^4)*S^4 + (- 351*R*b^3 - 729*R^3)*S^3 +
	+ (729*R^2*b^2 + 81*b^5)*S^2 - 243*R*b^4*S + 0
S*(S^3 + 9*b*x - 27*R)*(S^8 - 3*b*S^6 + 18*b^2*S^4 - 27*R*b*S^3 + (27*R^2 + 4*b^3)*S^2 - 27*R*b^2*S + 9*b^4)
#
S^8 - 3*b*S^6 + 18*b^2*S^4 - 27*R*b*S^3 + (27*R^2 + 4*b^3)*S^2 - 27*R*b^2*S + 9*b^4


### Derivation
# [old version]
(x^3*y^3 + x^3*z^3 + y^3*z^3) + # = E2^3 - 6*E3^2 - 3*E3*S*(S^2 - 2*E2) + 3*E3*(S^3 - 3*E2*S + 3*E3)
	+ x^4*y^2 + x^4*z^2+ x^2*y^4 + y^4*z^2 + y^2*z^4 + x^2*z^4 + # = E2^2 *(S^2 - 2*E2) - 2*E3*(E2*S - 2*E3) - 2*E3*(S^3 - 3*E2*S + 3*E3) - E3^2
	+ 2*x*y*z*(x^2*y + x^2*z + x*y^2 + y^2*z + x*z^2 + y*z^2) + # = 2*E3*(E2*S - 3*E3)
	+ x*y*z*(x^3 + y^3 + z^3) +  3*(x*y*z)^2 + b^3
#
E2^3 + 3*E3^2 - 3*E3*E2*S + E2^2*S^2 - 2*E3*S^3 + 4*E3*E2*S - 3*E3^2 - 2*E2^3 + 
	+ 2*E3*(E2*S - 3*E3) + E3*(S^3 - 3*E2*S + 3*E3) + 3*E3^2 + b^3
#
E2^3 - E2^2*S^2 + E3*S^3 - b[1]^3

### Solution:

solve.sysHt33 = function(R, b) {
	# only S8 used to compute the roots:
	coeff = c(1, 0, - 3*b[1], 0, 18*b[1]^2, - 27*R[1]*b[1], (27*R[1]^2 + 4*b[1]^3), - 27*R[1]*b[1]^2, 9*b[1]^4)
	S = roots(coeff)
	# exclude roots: x == y == z = S/3,
	# as they create numerical instability due to root multiplicity;
	# [should be actually excluded from P8]
	isEq = round0(S^3 + 9*b[1]*S - 27*R) == 0
	S = S[ ! isEq]
	print(S)
	### TODO: find robust roots!
	# - there are 2 sets, each of 24 roots;
	# - system of the 2nd set is unknown;
	#   e.g. {S1-4: + Det, S5-8: - Det}, {S1-4: - Det, S5-8: + Det};
	Det = sqrt(9*b[1]^2 + 24*S^4 + 96*b[1]*S^2 - 216*R[1]*S)
	E2 = c(- b[1] + 1/3 * Det, - b[1] - 1/3 * Det) / 4
	S = c(S, S)
	E3 = (6*E2*S - 2*S^3 - 2*b[1]*S + 6*R[1]) / 6
	x = sapply(1:length(S), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	y = (R - x^3) / b[1]
	z = (R - y^3) / b[1]
	return(cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z), S=as.vector(S)))
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


### Classic Polynomial

### Derivation:
y = (-x^3 + R) / b[1]
z = ((x^3 - R)^3 + R*b[1]^3) / b[1]^4
# b^12*z^3 + b^13*x - R*b^12 = 0
((x^3 - R)^3 + R*b[1]^3)^3 + b[1]^13*x - R*b[1]^12


x^27 - 9*R*x^24 + 36*R^2*x^21 + 3*R*(b[1]^3 - 28*R^2)*x^18 - 18*R^2*(b[1]^3 - 7*R^2)*x^15 + 9*R^3*(5*b[1]^3 - 14*R^2)*x^12 +
 + 3*(R^2*b[1]^6 - 20*R^4*b[1]^3 + 28*R^6)*x^9 - 9*R^3*(b[1]^6 - 5*R^2*b[1]^3 + 4*R^4)*x^6 +
 + 9*R^4*(b[1]^6 - 2*R^2*b[1]^3 + R^4)*x^3 + b[1]^13*x - (R*b[1]^12 - R^3*b[1]^9 + 3*R^5*b[1]^6 - 3*R^7*b[1]^3 + R^9)


coeff = c(1,0,0, - 9*R, 0,0, 36*R^2, 0,0, 3*R*(b[1]^3 - 28*R^2), 0,0, - 18*R^2*(b[1]^3 - 7*R^2), 0,0, 9*R^3*(5*b[1]^3 - 14*R^2),
	0,0, 3*(R^2*b[1]^6 - 20*R^4*b[1]^3 + 28*R^6), 0,0, - 9*R^3*(b[1]^6 - 5*R^2*b[1]^3 + 4*R^4), 0,0,
	9*R^4*(b[1]^6 - 2*R^2*b[1]^3 + R^4), 0, b[1]^13, - (R*b[1]^12 - R^3*b[1]^9 + 3*R^5*b[1]^6 - 3*R^7*b[1]^3 + R^9))
x = roots(coeff)

# - factorization: P3 * P24;
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


#################################
#################################


#################################
### Higher Power Correlations ###


################################
### x[i]^3 + b2*x[j]^2 + b1*x[j]

# x^3 + b2*y^2 + b1*y = R
# y^3 + b2*z^2 + b1*z = R
# z^3 + b2*x^2 + b1*x = R


### TODO

### Sum =>
# (x^3 + y^3 + z^3) + b2*(x^2 + y^2 + z^2) + b1*S = 3*R
S^3 - 3*E2*S + 3*E3 + b2*(S^2 - 2*E2) + b1*S - 3*R # = 0
S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 + 3*E3 - 3*R # = 0
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)

### Sum( (b[1]*x[i] + b[2]*x[i]^2) * ...) =>
b[2]*(x^5 + y^5 + z^5) + b[1]*(x^4 + y^4 + z^4) +
	+ b[1]*b[2]*(x^2*y + x*y^2 + x^2*z + x*z^2 + y^2*z + y*z^2) + b[1]^2*E2 + b[2]^2*(E2^2 - 2*E3*S) +
	- b[2]*R*(x^2+y^2+z^2) - b[1]*R*S # = 0
b[2]*(S^5 - 5*E2*S^3 + 5*E3*S^2 + 5*E2^2*S - 5*E2*E3) + b[1]*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ b[1]*b[2]*(E2*S - 3*E3) + b[1]^2*E2 + b[2]^2*(E2^2 - 2*E3*S) +
	- b[2]*R*(S^2 - 2*E2) - b[1]*R*S # = 0
b2*S^5 + b1*S^4 - 5*b2*E2*S^3 + 5*b2*E3*S^2 - 5*b2*E2*E3 + 4*b1*E3*S - 2*b2^2*E3*S +
	- 4*b1*E2*S^2 + 5*b2*E2^2*S - 2*b1*b2*E2*S +
	+ 2*b1*E2^2 + b2^2*E2^2 + b1^2*E2 - 2*b1*b2^2*E2 +
	+ b1*b2*S^3 + b1*b2^2*S^2 + b1^2*b2*S +
	- b2*R*S^2 + 2*b2*R*E2 - b1*R*S - 3*b1*b2*R # = 0

### Sum(y^3*...) =>
(x^3*y^3 + x^3*z^3 + y^3*z^3) + b2*(x^5 + y^5 + z^5) + b1*(x^4 + y^4 + z^4) +
	- R*(x^3 + y^3 + z^3) # = 0
(E2^3 - 3*E3*(E2*S - E3)) +
	+ b2*(S^5 - 5*E2*S^3 + 5*E3*S^2 + 5*E2^2*S - 5*E2*E3) +
	+ b1*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	- R*(S^3 - 3*E2*S + 3*E3) # = 0
E2^3 - 3*E3*E2*S + 3*E3^2 +
	+ b2*S^5 - 5*b2*E2*S^3 + 5*b2*E3*S^2 + 5*b2*E2^2*S - 5*b2*E2*E3 +
	+ b1*S^4 - 4*b1*E2*S^2 + 4*b1*E3*S + 2*b1*E2^2 +
	- R*S^3 + 3*R*E2*S - 3*R*E3 # = 0

E2Subst = 1701*R*S*b1*b2^4 - 243*R*S*b1^3 - 756*R*S*b2^6 + 918*R*S^2*b1*b2^3 - 243*R*S^2*b1^2*b2 + 288*R*S^2*b2^5 - 1782*R*S^3*b1*b2^2 + 648*R*S^3*b1^2 + 450*R*S^3*b2^4 - 486*R*S^4*b1*b2 - 945*R*S^4*b2^3 - 702*R*S^5*b2^2 - 1134*R*b1*b2^5 + 1620*R*b1^2*b2^3 - 729*R*b1^3*b2 + 729*R^2*S*b1*b2 - 486*R^2*S*b2^3 + 972*R^2*S^2*b2^2 - 729*R^2*b1*b2^2 + 378*S*b1^2*b2^5 - 540*S*b1^3*b2^3 + 243*S*b1^4*b2 + 630*S^2*b1*b2^6 - 1233*S^2*b1^2*b2^4 + 504*S^2*b1^3*b2^2 - 537*S^3*b1*b2^5 - 180*S^3*b1^2*b2^3 + 90*S^3*b1^3*b2 + 252*S^3*b2^7 - 516*S^4*b1*b2^4 + 540*S^4*b1^2*b2^2 - 189*S^4*b1^3 + 30*S^4*b2^6 + 576*S^5*b1*b2^3 - 27*S^5*b1^2*b2 - 144*S^5*b2^5 + 513*S^6*b1*b2^2 - 108*S^6*b1^2 + 186*S^6*b2^4 + 45*S^7*b1*b2 + 297*S^7*b2^3 + 90*S^8*b2^2;

E2Div = 1539*R*S*b1*b2^2 - 486*R*S*b1^2 - 459*R*S*b2^4 + 648*R*S^2*b1*b2 + 702*R*S^2*b2^3 + 810*R*S^3*b2^2 - 567*R*b1*b2^3 + 189*R*b2^5 - 729*R^2*b2^2 + 819*S*b1*b2^5 - 45*S*b1^2*b2^3 - 108*S*b1^3*b2 - 504*S*b2^7 + 954*S^2*b1*b2^4 - 621*S^2*b1^2*b2^2 + 108*S^2*b1^3 - 123*S^2*b2^6 - 810*S^3*b1*b2^3 + 162*S^3*b1^2*b2 + 312*S^3*b2^5 - 945*S^4*b1*b2^2 + 270*S^4*b1^2 - 369*S^4*b2^4 - 108*S^5*b1*b2 - 684*S^5*b2^3 - 225*S^6*b2^2 - 756*b1*b2^6 + 1458*b1^2*b2^4 - 1026*b1^3*b2^2 + 243*b1^4;
# E2 = - E2Subst / E2Div;

### Eq: has 730 monoms
# TODO: find way to factorize;


### Test
x^3 + b[2]*y^2 + b[1]*y # - R
y^3 + b[2]*z^2 + b[1]*z # - R
z^3 + b[2]*x^2 + b[1]*x # - R

### Debug
R = 2; b = c(3,-1);
x = -0.8426668811 - 0.6742950178i;
y =  0.4545225198 + 0.5403469156i;
z =  1.7559415104 + 0.3460249473i;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;




################################
### x[i]^3 + b2*x[j]^2 + b1*x[k]

# x^3 + b2*y^2 + b1*z = R
# y^3 + b2*z^2 + b1*x = R
# z^3 + b2*x^2 + b1*y = R

### Solution:

### Sum =>
# (x^3 + y^3 + z^3) + b2*(x^2 + y^2 + z^2) + b1*S = 3*R
# S^3 - 3*E2*S + 3*E3 + b2*(S^2 - 2*E2) + b1*S - 3*R = 0
# S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 + 3*E3 - 3*R = 0
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)

### Eq2:
# x^3 + b1*z = R - b2*y^2
# x^3 + b1*y + b1*z = R - b2*y^2 + b1*y
# (x^3 + b1*y + b1*z)^2 = (b2*y^2 - b1*y - R)^2
# x^6 + b1^2*y^2 + b1^2*z^2 + 2*b1*x^3*y + 2*b1*x^3*z + 2*b1^2*y*z =
#  b2^2*y^4 + b1^2*y^2 + R^2 - 2*b1*b2*y^3 - 2*b2*R*y^2 + 2*b1*R*y
### Sum =>
# (x^6+y^6+z^6) + 2*b1^2*(x^2+y^2+z^2) + 2*b1*sum(x^3*y) + 2*b1^2*E2 =
#   b2^2*(x^4+y^4+z^4) - 2*b1*b2*(x^3+y^3+z^3) + (b1^2 - 2*b2*R)*(x^2+y^2+z^2) + 2*b1*R*S + 3*R^2
### sum(x^3*y) = (x^2+y^2+z^2)*E2 - 3*E3*S
# TODO: ...

### Eq3:
# x^3 + b2*y^2 = R - b1*z
# x^3 + b2*y^2 + b2*z^2 = b2*z^2 - b1*z + R
# (x^3 + b2*y^2 + b2*z^2)^2 = (b2*z^2 - b1*z + R)^2
# x^6 + 2*b2*x^3*y^2 + 2*b2*x^3*z^2 + b2^2*y^4 + b2^2*z^4 + 2*b2^2*y^2*z^2 =
#  b2^2*z^4 - 2*b1*b2*z^3 + b1^2*z^2 - 2*b2*R*z^2 - 2*b1*R*z + R^2
### Sum =>
(x^6+y^6+z^6) + 2*b2*sum(x^3*y^2 + x^3*z^2) + 2*b2^2*(x^4+y^4+z^4) + 2*b2^2*(E2^2 - 2*E3*S) =
b2^2*(x^4+y^4+z^4) - 2*b1*b2*(x^3+y^3+z^3) + (b1^2 - 2*b2*R)*(x^2+y^2+z^2) - 2*b1*R*S + 3*R^2
### sum(x^3*y^2) = (E2^2 - 2*E3*S)*S - E3*E2
# TODO: ...


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

### [Temporary]
# TODO:
# - clean & integrate where appropriate;

# x^2 + b*x*z = R
# y^2 + b*y*z = R
# z^2 + b*x*y = R

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

