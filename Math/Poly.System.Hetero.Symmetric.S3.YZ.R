########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogenous Symmetric: Y*Z-Type
###
### draft v.0.1d-case2


### Hetero-Symmetric
### Polynomial Systems: 3 Variables
### Side-chain: Y*Z-Type

### Example:
x^n + b*y*z = R
y^n + b*x*z = R
z^n + b*x*y = R


######################

###############
### History ###
###############


### draft v.0.1d - v.0.1d-case:
# - solved type x*y*z:
#   x^3 + b*x*y*z + b1*x = R;
# - [DONE] Case x == y, but != z; [v.0.1d-case]
### draft v.0.1b - v.0.1c:
# - [started work] / solved extension:
#   x^3 + b*y*z + b1*x = R;
# - [DONE] Case x == y, but != z; [v.0.1d-case2]
### draft v.0.1a:
# - moved to this new file from file:
#   Poly.System.Hetero.Symmetric.S3.R;
### [old file]
# - solved Order 3:
#   x^3 + b*y*z = R;
# - [DONE] Case x == y, but != z; [v.0.1c]


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

################################
################################

#########################
### Type: x^n + b*y*z ###
#########################

###############
### Order 2 ###
###############

# x^2 + b*y*z = R
# y^2 + b*x*z = R
# z^2 + b*x*y = R

### Solution:

### TODO


########################
########################

###############
### Order 3 ###
###############

# x^3 + b*y*z = R
# y^3 + b*x*z = R
# z^3 + b*x*y = R

### Solution:

### Case 1:
# x = y = z
x^3 + b*x^2 - R # = 0

### Case 2:
# x = y, but != z
x^3 + b*x*z - R # = 0
z^3 + b*x^2 - R # = 0
# TODO!

### Case 3:
# x != y != z

### Diff =>
(x-y)*(x^2 + y^2 + x*y - b*z) # = 0
(x-z)*(x^2 + z^2 + x*z - b*y) # = 0
(y-z)*(y^2 + z^2 + y*z - b*x) # = 0
### Sum =>
2*(x^2+y^2+z^2) + E2 - b*S # = 0
2*S^2 - 3*E2 - b*S # = 0

# x*... =>
x^4 + b*x*y*z - R*x # = 0
# ...
### Diff: (x[i]*Eq[i] - x[j]*Eq[j]) =>
(x^4 - y^4) - R*(x - y) # = 0
(x-y)*(x^3 + y^3 + x^2*y + x*y^2 - R) # = 0
### Sum(Diff) =>
2*(x^3 + y^3 + z^3) + (E2*S - 3*E3) - 3*R # = 0
2*(S^3 - 3*E2*S + 3*E3) + E2*S - 3*E3 - 3*R
2*S^3 - 5*E2*S + 3*E3 - 3*R

### Sum =>
(x^3 + y^3 + z^3) + b*E2 - 3*R # = 0
S^3 - 3*E2*S + 3*E3 + b*E2 - 3*R

### Relations:
# 3*E2 = 2*S^2 - b*S
# 3*E3 = - (2*S^3 - 5*E2*S - 3*R)
# 3*E3 = - (S^3 - 3*E2*S + b*E2 - 3*R)

### =>
2*S^3 - 5*E2*S - 3*R - (S^3 - 3*E2*S + b*E2 - 3*R) # = 0
S^3 - 2*E2*S - b*E2 # = 0
3*S^3 - 2*(2*S^2 - b*S)*S - b*(2*S^2 - b*S) # = 0
S^3 - b^2*S # = 0

###########
### Case 2:
# x = y, but != z
x^3 + b*x*z - R # = 0
z^3 + b*x^2 - R # = 0

### Diff =>
(x-z)*(x^2 + z^2 + x*z - b*x) # = 0
x^2 + z^2 + x*z - b*x # = 0 # * z =>
z^3 + x^2*z + x*z^2 - b*x*z # + Eq 1 =>
x^3 + z^3 + x^2*z + x*z^2 - R
s^3 - 2*x*z*s - R
### Diff: x*Eq1 - z*Eq2 => [redundant]
x^4 - z^4 - R*(x-z) # = 0
(x-z)*(x^3 + z^3 + x*z*s - R) # = 0
### Diff: z*Eq1 - x*Eq2 =>
x*z*(x^2 - z^2) - b*x*(x-z)*s + R*(x-z) # = 0
x*z*s - b*x*s + R # = 0

### [classic approach]
(x^3 - R)^3 - b^4*x^5 + b^3*x^3*R # = 0
x^9 - 3*R*x^6 - b^4*x^5 + 3*R^2*x^3 + b^3*x^3*R - R^3
# (x^3 + b*x^2 - R) * P[6]
x^6 - b*x^5 + b^2*x^4 - (b^3 + 2*R)*x^3 + b*R*x^2 + R^2


### Solver:
solve.Htxy.S3P3 = function(R, b, debug=TRUE) {
	# S = c(0, b[1], -b[1]);
	S = -b[1]; # it seems this is the ONLY solution;
	if(debug) print(S);
	E2 = (2*S^2 - b[1]*S) / 3;
	E3 = - (2*S^3 - 5*E2*S - 3*R[1]) / 3;
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	if(any(round0(x) == 0)) print("Division by 0!")
	yz.s = S - x; yz = E3 / x;
	yz.d = sqrt(yz.s^2 - 4*yz + 0i);
	y = (yz.s + yz.d) / 2;
	z = (yz.s - yz.d) / 2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = rbind(sol, sol[,c(1,3,2)])
	### Case: x == y, but != z;
	coeff = c(1, - b[1], b[1]^2, - (b[1]^3 + 2*R[1]), b[1]*R[1], 0, R[1]^2);
	x = roots(coeff);
	z = (R - x^3) / (b[1]*x);
	sol21 = cbind(x=x, y=x, z=z);
	sol = rbind(sol, sol21);
	return(sol);
}

### Examples:

R = -1
b = 2
sol = solve.Htxy.S3P3(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b*y*z # - R
y^3 + b*x*z # - R
z^3 + b*x*y # - R


### Debug
R = 2; b = 3;
x = -0.2089400963 - 0.4200165444i;
y = -0.2089400963 - 0.4200165444i;
z = -0.5886981677 + 1.2138678312i;
S = x+y+z; E2 = x*(y+z)+y*z; E3 = x*y*z;


#########################

### Asymmetric Extension:

# x^3 + b*y*z + b1*x = R
# y^3 + b*x*z + b1*y = R
# z^3 + b*x*y + b1*z = R

### Solution:

### Case 3:
# x != y != z

### Diff =>
(x-y)*(x^2 + y^2 + x*y - b*z + b1) # = 0
(x-z)*(x^2 + z^2 + x*z - b*y + b1) # = 0
(y-z)*(y^2 + z^2 + y*z - b*x + b1) # = 0
### Sum =>
2*(x^2+y^2+z^2) + E2 - b*S + 3*b1 # = 0
2*S^2 - 3*E2 - b*S + 3*b1 # = 0

# x*... =>
x^4 + b*x*y*z + b1*x^2 - R*x # = 0
# ...
### Diff: (x[i]*Eq[i] - x[j]*Eq[j]) =>
(x^4 - y^4) + b1*(x^2 - y^2) - R*(x - y) # = 0
(x-y)*(x^3 + y^3 + x^2*y + x*y^2 + b1*(x+y) - R) # = 0
### Sum(Diff) =>
2*(x^3 + y^3 + z^3) + (E2*S - 3*E3) + 2*b1*S - 3*R # = 0
2*(S^3 - 3*E2*S + 3*E3) + E2*S - 3*E3 + 2*b1*S - 3*R
2*S^3 - 5*E2*S + 3*E3 + 2*b1*S - 3*R
### Diff(Diff) =>
(y^3 - z^3) + (x^2*y + x*y^2 - x^2*z - x*z^2) + b1*(y-z) # = 0
(y-z)*(x^2 + y^2 + z^2 + y*z + x*(y+z) + b1) # = 0
S^2 - E2 + b1 # = 0


### Sum =>
(x^3 + y^3 + z^3) + b*E2 + b1*S - 3*R # = 0
S^3 - 3*E2*S + 3*E3 + b*E2 + b1*S - 3*R

### Relations:
# E2 = S^2 + b1
# 3*E2 = 2*S^2 - b*S + 3*b1
# 3*E3 = - (2*S^3 - 5*E2*S + 2*b1*S - 3*R)
# 3*E3 = - (S^3 - 3*E2*S + b*E2 + b1*S - 3*R)

### =>
3*S^2 + 3*b1 - (2*S^2 - b*S + 3*b1) # = 0
S^2 + b*S # = 0

### Case 2:
# x == y, but != z;
x^3 + b*x*z + b1*x - R # = 0
z^3 + b*x^2 + b1*z - R # = 0

### [classic approach]
b^3*x^3*z^3 + b^4*x^5 + b^3*b1*z*x^3 - b^3*R*x^3 # = 0
(x^3 + b1*x - R)^3 - b^4*x^5 + b^2*b1*(x^3 + b1*x - R)*x^2 + b^3*R*x^3 # = 0
(x^3 + b*x^2 + b1*x - R)^3 - 3*b*x^2*(x^3 + b*x^2 + b1*x - R)^2 +
	+ 3*b^2*x^4*(x^3 + b*x^2 + b1*x - R) + b^2*b1*(x^3 + b*x^2 + b1*x - R)*x^2 +
	- b^3*(x^3 + b*x^2 + b1*x - R)*x^3 # = 0
(x^3 + b*x^2 + b1*x - R) *
	((x^3 + b*x^2 + b1*x - R)^2 - 3*b*x^2*(x^3 + b*x^2 + b1*x - R) +
	+ 3*b^2*x^4 + b^2*b1*x^2 - b^3*x^3) # = 0
# P[3] * P[6]
(x^3 + b*x^2 + b1*x - R) *
	(x^6 - b*x^5 + b^2*x^4 + 2*b1*x^4 - (b^3 + b*b1 + 2*R)*x^3 +
	+ (b^2*b1 + b1^2 + b*R)*x^2 - 2*b1*R*x + R^2) # = 0


### Solver:
solve.Htxy.S3P3 = function(R, b, bc, debug=TRUE) {
	S = c(-b[1]); # it seems this is the ONLY solution;
	if(debug) print(S);
	E2 = S^2 + bc[1];
	E3 = - (2*S^3 - 5*E2*S + 2*bc[1]*S - 3*R[1]) / 3;
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	if(any(round0(x) == 0)) print("Division by 0!")
	yz.s = S - x; yz = E3 / x;
	yz.d = sqrt(yz.s^2 - 4*yz + 0i);
	y = (yz.s + yz.d) / 2;
	z = (yz.s - yz.d) / 2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = rbind(sol, sol[,c(1,3,2)])
	### Case: x == y, but != z
	coeff = c(1, - b[1], b[1]^2 + 2*bc[1], - (b[1]^3 + b[1]*bc[1] + 2*R),
		(b[1]^2*bc[1] + bc[1]^2 + b[1]*R), - 2*bc[1]*R, R^2);
	x = roots(coeff);
	z = - (x^3 + bc[1]*x - R[1]) / (b*x);
	sol2 = cbind(x=x, y=x, z=z);
	sol = rbind(sol, sol2);
	return(sol);
}

### Examples:

R = -1
b = 2
bc = 3
sol = solve.Htxy.S3P3(R, b, bc)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b*y*z + bc*x # - R
y^3 + b*x*z + bc*y # - R
z^3 + b*x*y + bc*z # - R


### Derivation:
p = list(
	x = c(3,2,1,0),
	b = c(0,1,0,0),
	b1= c(0,0,1,0),
	R = c(0,0,0,1),
	coeff = c(1,1,1,-1)
	)


################################
################################

###########################
### Type: x^n + b*x*y*z ###
###########################

###############
### Order 3 ###
###############

# x^3 + b*x*y*z + b1*x = R
# y^3 + b*x*y*z + b1*y = R
# z^3 + b*x*y*z + b1*z = R

### Solution:

### Case 3:
# x != y != z

### Diff =>
(x-y)*(x^2 + y^2 + x*y + b1) # = 0
(x-z)*(x^2 + z^2 + x*z + b1) # = 0
(y-z)*(y^2 + z^2 + y*z + b1) # = 0
### Diff(Diff) =>
(y-z)*S # = 0
# S = 0!

### Sum(Diff) =>
2*(S^2 - 2*E2) + E2 + 3*b1 # = 0
# 3*E2 = 2*S^2 + 3*b1

### Sum =>
S^3 - 3*E2*S + 3*E3 + 3*b*E3 + b1*S - 3*R # = 0
# 3*(b+1)*E3 = - (S^3 - 3*E2*S + b1*S - 3*R)


###########
### Case 2:
# x = y, but != z
x^3 + b*x^2*z + b1*x - R # = 0
z^3 + b*x^2*z + b1*z - R # = 0

### [classic approach]
(x^3 + b1*x - R)^3 + b^3*x^9 + b^3*b1*x^7 + b^2*b1*x^7 +
	+ b^2*b1^2*x^5 - b^2*b1*R*x^4
(b^3+1)*x^9 + (b^3 + b^2 + 3)*b1*x^7 - 3*R*x^6 + (b^2+3)*b1^2*x^5 - (b^2+6)*b1*R*x^4 +
	+ b1^3*x^3 + 3*R^2*x^3 - 3*b1^2*R*x^2 + 3*b1*R^2*x - R^3
### ((b+1)*x^3 + b1*x - R) * P[6]
(1 - b + b^2)*x^6 +
	+ b1*(b^2 - b + 2)*x^4 + R*(b - 2)*x^3 +
	+ b1^2*x^2 - 2*R*b1*x + R^2


### Solver:
solve.Htxyz.S3P3 = function(R, b, bc, debug=TRUE) {
	S = c(0);
	if(debug) print(S);
	E2 = (2*S^2 + 3*bc[1]) / 3; # bc[1]
	E3 = - (S^3 - 3*E2*S + bc[1]*S - 3*R) / (3*(b[1]+1));
	len = length(S);
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])));
	S = rep(S, each=3); E3 = rep(E3, each=3);
	if(any(round0(x) == 0)) print("Division by 0!")
	yz.s = S - x; yz = E3 / x;
	yz.d = sqrt(yz.s^2 - 4*yz + 0i);
	y = (yz.s + yz.d) / 2;
	z = (yz.s - yz.d) / 2;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z));
	sol = rbind(sol, sol[,c(1,3,2)])
	
	### Case: x == y, but != z;
	bb = b[1]^2 - b[1] + 1;
	coeff = c(bb, 0, bc[1]*(bb + 1), R[1]*(b[1] - 2), bc[1]^2, - 2*R[1]*bc[1], R[1]^2);
	x = roots(coeff);
	z = - (x^3 + bc[1]*x - R[1]) / (b[1]*x^2);
	sol2 = cbind(x=x, y=x, z=z);
	sol = rbind(sol, sol2);
	return(sol);
}

### Examples:

R = -1
b = 2
bc = 3
sol = solve.Htxyz.S3P3(R, b, bc)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^3 + b*x*y*z + bc*x # - R
y^3 + b*x*y*z + bc*y # - R
z^3 + b*x*y*z + bc*z # - R

