
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogenous Symmetric
###
### draft v.0.3b


### Hetero-Symmetric
### Polynomial Systems: 3 Variables

### Example:
x^n + P(x, y, z) = R
y^n + P(y, z, x) = R
z^n + P(z, x, y) = R

####################

###############
### History ###

### draft v.0.3b:
# - solved: x^3 + b*y*z = R;
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


################

################
### Introduction

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

### Classic Polynomial: P8 or P6 (for S == 0)
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

### Shifted roots:
### x[i]^2 + s*x[i] + b*x[i+1]

# x^2 + s*x + b1*y = R
# y^2 + s*y + b1*z = R
# z^2 + s*z + b1*x = R

### Solution

### Fast Solution:
# - shift all roots "back" & use solution to simple system: x^2 + b1*y = R;
# (x + s/2)^2 + b1*(y + s/2) = R + s^2/4 + b1*s/2
# (y + s/2)^2 + b1*(z + s/2) = R + s^2/4 + b1*s/2
# (z + s/2)^2 + b1*(x + s/2) = R + s^2/4 + b1*s/2

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
# 0 == 0

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

### Bug corrected: + 6*b1^3;

### Example
b = 3
s = 1
R = 1
# TODO:
# - clarify reason for the 3 incorrect roots
#   & find way to remove 3 incorect roots;
coeff = c(2, (9*s+b[1]), (b[1]^2 - 2*R + 13*s^2 + 2*s*b[1]), 6*s^3 + 6*b[1]^3 - 3*s*R - 3*b[1]*R)
x.sum = roots(coeff)
E3 = (x.sum^3 + (3*s + 2*b[1])*x.sum^2 + (2*s-b[1])*(s+b[1])*x.sum - 7*R*x.sum - 6*s*R + 3*b[1]*R)/6
E2 = (x.sum^2 + (s+b[1])*x.sum - 3*R)/2
x = as.vector(sapply(1:length(x.sum), function(id) roots(c(1, -x.sum[id], E2[id], -E3[id]))))
y = (R - x^2 - s*x)/b[1]
z = (R - y^2 - s*y)/b[1]
sol = cbind(x, y, z)
sol = rbind(sol, sol[,c(2,3,1)], sol[,c(3,2,1)])
x = sol[,1]; y = sol[,2]; z = sol[,3]
sol

### Test
x^2 + s*x + b[1]*y
y^2 + s*y + b[1]*z
z^2 + s*z + b[1]*x

### Classic polynomial
round0.p(poly.calc(sol[4:9,1]))


######################
######################

####################
### x[i]^2 + b*x[-i]

# x^2 + b1*y*z = R
# y^2 + b1*x*z = R
# z^2 + b1*x*y = R

### Solution

# - degenerate P4;

# Special Case: b1 = 2
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

# Case: x = y
# x^2 + b1*x*z = R
# z^2 + b1*x^2 = R
# =>
# b1*z = R/x - x
# b1^2*z^2 + b1^3*x^2 = b1^2*R
# x^2 - 2*R + R^2/x^2 + b1^3*x^2 - b1^2*R = 0
# (b1^3+1)*x^4 - R*(b1^2 + 2)*x^2 + R^2 = 0

### Example:
b = 1
R = 2
#
x = roots(c((b[1]^3+1), 0, - R*(b[1]^2 + 2), 0, R^2))
y = x
z = (R - x^2)/y/b[1]
sol = round0(cbind(x, y, z))
sol

### Test
x^2 + b[1]*y*z
y^2 + b[1]*x*z
z^2 + b[1]*x*y


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


#########################

### Prod-Type: Order 3
### x[i]^3 + b*x[-i]

# x^3 + b1*y*z = R
# y^3 + b1*x*z = R
# z^3 + b1*x*y = R

### Solution

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

solve.htyz.S3P3 = function(R, b, doEq = FALSE, tol=1E-8) {
	b1 = b[1]; R = R[1];
	coeff = c(1, - 3*b1, 7*b1^2, -5*b1^3, - 27*R*b1, 18*R*b1^2, (5*R*b1^3 + 27*R^2))
	S = roots(coeff)
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
	yz = E3 / x; # (R - x^3)/b1 # robust ???
	yz.s = S - x;
	#
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d)/2;
	z = (yz.s - yz.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	### y == z
	if(doEq) {
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
# S^2 - 2*E2 + (b1+b2)*S = 3*R
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
# (x^3 + y^3 + z^3) + (b1+b2)*S - 3*R = 0
# S^3 - 3*E2*S + 3*E3 + (b1+b2)*S - 3*R = 0
# 3*E3 = -(S^3 - 3*E2*S + (b1+b2)*S - 3*R)

### Sum(x[i]*P(x)) =>
# (x^4 + y^4 + z^4) + (b1+b2)*E2 - R*S = 0
# S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 + (b1+b2)*E2 - R*S = 0

### Eq 3:
# (x^3 + b1*y + b1*z)^2 = (R + b1*z - b2*z)^2


### Test
x^3 + b[1]*y + b[2]*z
y^3 + b[1]*z + b[2]*x
z^3 + b[1]*x + b[2]*y


#########################

#########################
#########################

### High-Power Terms: > 1

### x[i]^2 + x[j]^2 + b*(x[i] + x[j])

# x^2 + y^2 + b1*(x+y) = R
# y^2 + z^2 + b1*(y+z) = R
# x^2 + z^2 + b1*(x+z) = R

# trivial solution: x = y = z;

### Solution

### Diff =>
# x^2 - z^2 = -b1*(x - z)
# (x-z)*(x + z + b1) = 0

# Case: x != y != z
# - has NO solutions;

# Case x = y && x != z
# trivial;
# x^2 + b1*x = R/2
# z^2 + b1*z = R/2 # x & z are the conjugate roots;

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




########################
########################

### High-Power Terms: > 1

### x[i]^2 + x[j]^2 + b*x[k]

# x^2 + y^2 + b1*z = R
# y^2 + z^2 + b1*x = R
# x^2 + z^2 + b1*y = R

# trivial solution: x = y = z;
# trivial system;

### Solution

### Diff =>
# x^2 - z^2 = b1*(x - z)
# (x-z)*(x + z - b1) = 0

# Case: x != y != z
# - has NO solutions;

# Case x = y && x != z
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

### Problems:
# - Difference works well for systems with 2 variables;
# - but it does NOT work well in systems with 3 variables;

###############
### Order 3 ###
###############

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


#######################

#######################
### Mixt-Order: 2+1 ###

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



##########################

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
# (x-y)*(x^2 + y^2 + x*y + b2*x*y*z + b1) = 0
# (y-z)*(y^2 + z^2 + y*z + b2*x*y*z + b1) = 0
# (z-x)*(x^2 + z^2 + x*z + b2*x*y*z + b1) = 0
# Case: x != y != z =>
# x^2 + y^2 + x*y + b2*x*y*z + b1 = 0
# y^2 + z^2 + y*z + b2*x*y*z + b1 = 0
# x^2 + z^2 + x*z + b2*x*y*z + b1 = 0
# =>
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



##########################

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
# S^3 - 3*E2*S + 3*E3 + b2*(S^2 - 2*E2) + b1*S - 3*R = 0
# S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 + 3*E3 - 3*R = 0
# 3*E3 = -(S^3 + b2*S^2 + b1*S - 3*E2*S - 2*b2*E2 - 3*R)

### Sum( (b[1]*x[i] + b[2]*x[i]^2) * ...) =>
b[2]*(x^5 + y^5 + z^5) + b[1]*(x^4 + y^4 + z^4) +
	+ b[1]*b[2]*(x^2*y + x*y^2 + x^2*z + x*z^2 + y^2*z + y*z^2) + b[1]^2*E2 + b[2]^2*(E2^2 - 2*E3*S) +
	- b[2]*R*(x^2+y^2+z^2) - b[1]*R*S # = 0



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

