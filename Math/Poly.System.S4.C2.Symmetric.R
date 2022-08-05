########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
### with Additional Symmetry
###
### draft v.0.1d


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


####################
####################

####################
### Higher Power ###
####################

### Eq 2 & Eq 4: Equal Power
# - slightly easier;
# (x1*y1)^n + (x2*y2)^n = R2
# (x1*y2)^n + (x2*y1)^n = R4

### A + B =>
(x1^n + x2^n)*(y1^n + y2^n) - R2 - R4 # = 0

### A * B =>
p2^n*(x1^(2*n) + x2^(2*n)) + p1^n*(y1^(2*n) + y2^(2*n)) - R2*R4 # = 0


### Case: n = 2

### A + B =>
(x1^2 + x2^2)*(y1^2 + y2^2) - R2 - R4 # = 0
(s1^2 - 2*p1)*(s2^2 - 2*p2) - R2 - R4 # = 0
### A * B =>
p2^2*(x1^4 + x2^4) + p1^2*(y1^4 + y2^4) - R2*R4 # = 0
p2^2*(s1^4 - 4*p1*s1^2 + 2*p1^2) + p1^2*(s2^4 - 4*p2*s2^2 + 2*p2^2) - R2*R4 # = 0
p2^2*(s1^4 - 4*p1*s1^2) + p1^2*(s2^4 - 4*p2*s2^2) + 4*p1^2*p2^2 - R2*R4 # = 0


### Transformed System
# (Case: n1 = 1)
s1 + s2 - R1 # = 0
(s1^2 - 2*p1)*(s2^2 - 2*p2) - R2 - R4 # = 0
s1*p2 + s2*p1 - R3 # = 0
p2^2*(s1^4 - 4*p1*s1^2) + p1^2*(s2^4 - 4*p2*s2^2) + 4*p1^2*p2^2 - R2*R4 # = 0


### Solver:

solve.S4C2.Var_x1y1P2 = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S4C2.Var_x1y1P2(R);
	s1 = roots(coeff);
	if(debug) print(s1);
	s2 = R[1] - s1;
	p1 = solveP.S4C2.Var_x1y1P2.x0(R, s1);
	p2 = (R[3] - s2*p1) / s1;
	# Step 3
	solve2 = function(id, s, p) {
		r = roots(c(1, -s[id], p[id]));
	}
	len = length(s1);
	x12 = lapply(seq(len), solve2, s1, p1);
	x12 = matrix(unlist(x12), nrow=2);
	x12 = t(x12);
	#
	x12sq = x12^2; x12q = x12sq^2;
	y1sq = (R[2]*x12sq[,1] - R[4]*x12sq[,2]) / (x12q[,1] - x12q[,2]);
	y2sq = (R[2] - x12sq[,1]*y1sq) / x12sq[,2];
	yd = (y1sq - y2sq) / s2;
	y1 = (s2 + yd) / 2;
	y2 = (s2 - yd) / 2;
	#
	sol = cbind(x12, y1, y2);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(sol);
}
coeff.S4C2.Var_x1y1P2 = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(1, - 8*R1, 26*R1^2, - 42*R1^3,
		- 12*R3*R1 + 28*R1^4 - 8*R4 - 8*R2,
		72*R3*R1^2 + 14*R1^5 + 48*R1*R4 + 48*R1*R2,
		- 16*R3^2 - 180*R3*R1^3 - 42*R1^6 - 116*R1^2*R4 - 116*R1^2*R2,
		80*R3^2*R1 + 240*R3*R1^4 + 34*R1^7 + 140*R1^3*R4 + 140*R1^3*R2,
		- 160*R3^2*R1^2 - 180*R3*R1^5 - 13*R1^8 + 16*R3*R1*R4 - 80*R1^4*R4 + 14*R4^2 +
			+ 16*R3*R1*R2 - 80*R1^4*R2 + 36*R4*R2 + 14*R2^2,
		160*R3^2*R1^3 + 72*R3*R1^6 + 2*R1^9 - 64*R3*R1^2*R4 + 8*R1^5*R4 - 56*R1*R4^2 +
			- 64*R3*R1^2*R2 + 8*R1^5*R2 - 144*R1*R4*R2 - 56*R1*R2^2,
		- 80*R3^2*R1^4 - 12*R3*R1^7 + 96*R3*R1^3*R4 + 12*R1^6*R4 + 82*R1^2*R4^2 +
			+ 96*R3*R1^3*R2 + 12*R1^6*R2 + 220*R1^2*R4*R2 + 82*R1^2*R2^2,
		16*R3^2*R1^5 - 64*R3*R1^4*R4 - 4*R1^7*R4 - 50*R1^3*R4^2 - 64*R3*R1^4*R2 +
			- 4*R1^7*R2 - 156*R1^3*R4*R2 - 50*R1^3*R2^2,
		16*R3*R1^5*R4 - 4*R3*R1*R4^2 + 8*R1^4*R4^2 - 8*R4^3 + 16*R3*R1^5*R2 + 8*R3*R1*R4*R2 +
			+ 48*R1^4*R4*R2 + 8*R4^2*R2 - 4*R3*R1*R2^2 + 8*R1^4*R2^2 + 8*R4*R2^2 - 8*R2^3,
		8*R3*R1^2*R4^2 + 2*R1^5*R4^2 + 16*R1*R4^3 - 16*R3*R1^2*R4*R2 - 4*R1^5*R4*R2 +
			- 16*R1*R4^2*R2 + 8*R3*R1^2*R2^2 + 2*R1^5*R2^2 - 16*R1*R4*R2^2 + 16*R1*R2^3,
		- 4*R3*R1^3*R4^2 - 8*R1^2*R4^3 + 8*R3*R1^3*R4*R2 + 8*R1^2*R4^2*R2 - 4*R3*R1^3*R2^2 +
			+ 8*R1^2*R4*R2^2 - 8*R1^2*R2^3,
		0,
		R4^4 - 4*R4^3*R2 + 6*R4^2*R2^2 - 4*R4*R2^3 + R2^4
	);
	return(coeff);
}
solveP.S4C2.Var_x1y1P2.x0 = function(R, s1) {
	p1 = coeff.S4C2.Var_x1y1P2.x0(R);
	p1 = sapply(s1, function(s1) {
		sum(p1 * s1^seq(11, 0, by=-1));
	})
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	div = 2*s1^9 - 13*s1^8*R1 + 36*s1^7*R1^2 - 55*s1^6*R1^3 + 50*s1^5*R1^4 - 27*s1^4*R1^5 +
		+ 8*s1^3*R1^6 - s1^2*R1^7;
	p1 = p1 / (64 * div);
	return(p1);
}
coeff.S4C2.Var_x1y1P2.x0 = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(48, - 336*R1, 1008*R1^2, - 64*R3 - 1680*R1^3,
		320*R3*R1 + 1680*R1^4 - 64*R4 - 64*R2,
		- 640*R3*R1^2 - 1008*R1^5 + 320*R1*R4 + 320*R1*R2,
		640*R3*R1^3 + 336*R1^6 - 640*R1^2*R4 - 640*R1^2*R2,
		- 320*R3*R1^4 - 48*R1^7 + 640*R1^3*R4 + 640*R1^3*R2,
		64*R3*R1^5 - 320*R1^4*R4 + 16*R4^2 - 320*R1^4*R2 - 32*R4*R2 + 16*R2^2,
		64*R1^5*R4 - 48*R1*R4^2 + 64*R1^5*R2 + 96*R1*R4*R2 - 48*R1*R2^2,
		48*R1^2*R4^2 - 96*R1^2*R4*R2 + 48*R1^2*R2^2,
		- 16*R1^3*R4^2 + 32*R1^3*R4*R2 - 16*R1^3*R2^2
	);
	return(coeff);
}

# TODO:
# Case: s2 = 0; s1 = R1;

### Examples:

R = c(2,3,-1,5)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")


### Ex 2:
R = c(2,3,-1,-1)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")


##########

### Debug:
x1 =  1.1725398503 + -1.4433080165i;
x2 = -0.7162047206 + -0.2205315378i;
y1 =  0.7158321338 + 0.7331960503i;
y2 =  0.8278327365 + 0.9306435040i;
sol = cbind(x1,x2,y1,y2);


R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;


### Derivation
p1 = toPoly.pm("s1 + s2 - R1")
p2 = toPoly.pm("s1*p2 + s2*p1 - R3")
p3 = toPoly.pm("(s1^2 - 2*p1)*(s2^2 - 2*p2) - R2 - R4")
p4 = toPoly.pm("p2^2*(s1^4 - 4*p1*s1^2) + p1^2*(s2^4 - 4*p2*s2^2) + 4*p1^2*p2^2 - R2*R4")

pR = solve.lpm(p1, p2, p3, p4, xn=c("s2", "p2", "p1"))
sort(unique(pR[[3]]$Rez$s1))

pD = div.pm(pR[[3]]$Rez, toPoly.pm("(s1 - R1)^7"), "s1")
p = toCoeff(pD$Rez, xn="s1", print=TRUE)


