########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
### with Additional Symmetry
###
### draft v.0.1b


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


####################
####################

##############
### Theory ###
##############

### Symmetric Eq 2
# n2 = n3 = n

### System:
# x1^n1 + x2^n1 + y1^n1 + y2^n1 = R1
# (x1*y1)^n + (x2*y2)^n = R2
# x1*x2*y1 + x1*x2*y2 + x1*y1*y2 + x2*y1*y2 = R3
# x1*x2*y1*y2 = R4

### Symmetries:
# - if (x1, x2, y1, y2) is a solution,
#   so is also the C2 permutation: (x2, x1, y2, y1);
# - there is additional symmetry:
#   the permutation (y1, y2, x1, x2) & its C2 permutation
#   are also solutions;

### Note:
# - for the non-symmetric types, see file:
#   Poly.System.S4.C2.R;


#############

### Solution:

### Step 1:
# (Eq 2) * (x1*y1)^n =>
(x1*y1)^(2*n) - R2*(x1*y1)^n + R4^n # = 0
# - solve for (x1*y1) & (x2*y2);

### Step 2:
### Eq 1 & Eq 3:

### Case 1: n1 = 1;
(x1 + y1) + (x2 + y2) - R1 # = 0
(x2*y2)*(x1 + y1) + (x1*y1)*(x2 + y2) - R3 # = 0
# - solve for: (x1 + y1) & (x2 + y2);

### Step 3:
### Solve 2x2:
# x1 + y1 = s1;
# x1 * y1 = p1;
# and:
# x2 + y2 = s2;
# x2 * y2 = p2;


### Step 2:
### Case 2: n1 > 1

### Case: n1 = 2
(x1^2 + y1^2) + (x2^2 + y2^2) - R1 # = 0
(x2*y2)*(x1 + y1) + (x1*y1)*(x2 + y2) - R3 # = 0

### Eq s:
s1^2 + s2^2 - 2*(x1*y1) - 2*(x2*y2) - R1 # = 0
(x2*y2)*s1 + (x1*y1)*s2 - R3 # = 0
# - solve for: s1 & s2;
# Note:
# - Eq s (2) is NOT symmetric;
### Substitution =>
(x1*y1)^2*s1^2 + ((x2*y2)*s1 - R3)^2 - 2*(x1*y1)^3 - R1*(x1*y1)^2 - 2*R4*(x1*y1) # = 0
((x1*y1 + x2*y2)^2 - 2*R4)*s1^2 - 2*R3*(x2*y2)*s1 +
	+ R3^2 - 2*(x1*y1)^3 - R1*(x1*y1)^2 - 2*R4*(x1*y1) # = 0


### Case: n1 > 2
# - similar, but more complicated;
### Ex: n1 = 3
(x1^3 + y1^3) + (x2^3 + y2^3) - R1 # = 0
(x2*y2)*(x1 + y1) + (x1*y1)*(x2 + y2) - R3 # = 0
# =>
s1^3 + s2^3 - 3*(x1*y1)*s1 - 3*(x2*y2)*s2 - R1 # = 0
(x2*y2)*s1 + (x1*y1)*s2 - R3 # = 0

# TODO: implement & check;


### Solver

solve.S4C2 = function(R, n, debug=TRUE) {
	if(length(n) > 2) stop("Ht: Not yet implemented!");
	if(n[1] > 3) stop("Only up to power 3 implemented!");
	len = 2*n[2];
	# Step 1:
	r = roots(c(1, - R[2], R[4]^n[2]));
	m = unity(n[2], all=TRUE);
	r = rootn(r, n=n[2]);
	p1 = lapply(r, function(r) r * m);
	p1 = unlist(p1);
	p2 = R[4] / p1;
	if(debug) { print(p1); print(p2); }
	# Step 2:
	if(n[1] == 1) {
		s1 = sapply(seq(len), function(id) {
			(p1[id]*R[1] - R[3]) / (p1[id] - p2[id]);
		})
		s2 = R[1] - s1;
	} else stop("TODO");
	# Step 3:
	solve2 = function(id, s, p) {
		r = roots(c(1, -s[id], p[id]));
	}
	xy1 = lapply(seq(len), solve2, s1, p1);
	xy1 = matrix(unlist(xy1), nrow=2);
	xy2 = lapply(seq(len), solve2, s2, p2);
	xy2 = matrix(unlist(xy2), nrow=2);
	sol = rbind(xy1, xy2);
	sol = t(sol); sol = sol[, c(1,3,2,4)];
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(sol);
}

### Examples:

R = c(1,2,-1,3)
n = c(1, 2)
sol = solve.S4C2(R, n=n)

test.S4C2(sol, n=n)


### Ex 2:
R = c(3,-1,2,3)
n = c(1, 2)
sol = solve.S4C2(R, n=n)

test.S4C2(sol, n=n)


########################
########################

################
### Variants ###
################

### Variants Eq 4:
# - break the decomposition into (x1*y1) & (x2*y2);
# - but the original symmetry is still preserved;


### System:
# x1^n1 + x2^n1 + y1^n1 + y2^n1 = R1
# (x1*y1)^n + (x2*y2)^n = R2
# x1*x2*y1 + x1*x2*y2 + x1*y1*y2 + x2*y1*y2 = R3
# x1*y2 + x2*y1 = R4


### Case: n = 1
x1*y1 + x2*y2 - R2 # = 0
x1*x2*(y1 + y2) + y1*y2*(x1 + x2) - R3 # = 0
x1*y2 + x2*y1 - R4 # = 0

### Solution:

# let:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;
S  = s1 + s2;
E4 = p1 * p2;

### Eq 2 + Eq 3 =>
s1*s2 - R2 - R4 # = 0

### Eq 2 * Eq 3 =>
p2*(s1^2 - 2*p1) + p1*(s2^2 - 2*p2) - R2*R4 # = 0
p2*s1^2 + p1*s2^2 - 4*p1*p2 - R2*R4 # = 0

### Transformed System:
# Eq 1: n1 = 1 vs n1 > 1;
s1 + s2 - R1 # = 0
s1*s2 - R2 - R4 # = 0
p2*s1 + p1*s2 - R3 # = 0
p2*s1^2 + p1*s2^2 - 4*p1*p2 - R2*R4 # = 0


### Solver:

solve.S4C2.SymVar1 = function(R, debug=TRUE, all=FALSE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	# Step 1:
	s1 = roots(c(1, -R1, R2 + R4));
	s2 = R1 - s1;
	if(debug) print(s1);
	len = length(s1);
	# Step 2:
	b1 = s2^2*s1 - s2*s1^2 - 4*R3;
	b0 = R3*s1^2 - s1*R4*R2;
	p1 = sapply(seq(len), function(id) {
		roots(c(4*s2[id], b1[id], b0[id]));
	});
	p1 = as.vector(p1);
	s1 = rep(s1, each = 2);
	s2 = rep(s2, each = 2);
	p2 = (R3 - p1*s2) / s1;
	# Step 3:
	solve2 = function(id, s, p) {
		roots(c(1, -s[id], p[id]));
	}
	len = length(s1);
	x12 = lapply(seq(len), solve2, s1, p1);
	x12 = matrix(unlist(x12), nrow=2);
	y12 = lapply(seq(len), solve2, s2, p2);
	y12 = matrix(unlist(y12), nrow=2);
	sol = rbind(x12, y12);
	sol = t(sol);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(sol);
}

### Examples

###
R = c(-4,3,-1,-2)
sol = solve.S4C2.SymVar1(R, all=T)

test.S4C2.Var(sol, n = c(1,1), type="x1y2")


###############
### Derivation:
4*s2*p1^2 + (s2^2*s1 - s2*s1^2 - 4*R3)*p1 + R3*s1^2 - s1*R4*R2 # = 0

