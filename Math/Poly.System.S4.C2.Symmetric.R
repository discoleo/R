########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
### with Additional Symmetry
###
### draft v.0.1k-edit


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

### Power: only Eq 2

### System: Eq 2 & Eq 4
# (x1*y1)^n + (x2*y2)^n = R2
# x1*y2 + x2*y1 = R4


#############
### n = 2 ###
#############

### System:
# x1 + x2 + y1 + y2 = R1
# (x1*y1)^2 + (x2*y2)^2 = R2
# x1*x2*(y1 + y2) + y1*y2*(x1 + x2) = R3
# x1*y2 + x2*y1 = R4

### Solution:

### Transform 1:
s1 + s2 - R1 # = 0
s1^2*s2^2 - 2*R4*s1*s2 - 2*p1*p2 - R2 + R4^2 # = 0
s1*p2 + s2*p1 - R3 # = 0
p1*s2^2+ p2*s1^2 - 4*p1*p2 - R4*s1*s2 + R4^2 # = 0

### Transform 2:
S - R1 # = 0
ps^2 - 2*R4*ps - 2*E4 - R2 + R4^2 # = 0
ps*sp^2 + E4*S^2 - 4*E4*ps - R3*sp*S + R3^2 # = 0
(4*E4 + R4*ps - R4^2)^2 - (4*E4 + R4*ps - R4^2)*sp*(S^2 - 2*ps) +
	+ E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0

### Solver:

solve.S4C2.SymVar1.P2P1 = function(R, debug = TRUE, all = FALSE) {
	coeff = coeff.S4C2.SymVar1.P2P1(R);
	ps = roots(coeff);
	if(debug) print(ps);
	### Step 2:
	s1 = sapply(ps, function(ps) {
		roots(c(1, - R[1], ps));
	});
	s1 = as.vector(s1);
	if(debug) print(s1);
	s2 = R[1] - s1;
	ps = rep(ps, each=2);
	# p:
	sp = solve.S4C2.SymVar1.P2P1sp(ps, R);
	p1 = (s1*sp - R[3]) / (s1 - s2);
	p2 = sp - p1;
	### Step 3:
	len = length(s1);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, - s1[id], p1[id]));
	});
	x12 = t(x12);
	# y:
	y1 = (x12[,1] * s2 - R[4]) / (x12[,1] - x12[,2]);
	y2 = s2 - y1;
	#
	sol = cbind(x12, y1, y2);
	if(all) sol = rbind(sol, sol[, c(2,1,4,3)]);
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(sol)
}
coeff.S4C2.SymVar1.P2P1 = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	#
	return(c(4, R1^2 - 16*R4, - 2*R1^2*R4 - 4*R1*R3 + 22*R4^2 - 12*R2,
		2*R3^2 - 12*R4^3 + 24*R2*R4 + R1^2*R4^2 - R1^2*R2 + 6*R1*R3*R4,
		2*(R4^2 - 2*R2)*(R4^2 - 2*R2 - R1*R3)));
}
# [old]
coeff.S4C2.SymVar1.P2P1.old = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	# R = c(4,1,1,3)
	if(all( (R - c(4,1,1,3)) == 0)) {
		return(c(2, -16, 37, -25, 21));
	}
	# R = c(2,1,-1,3)
	if(all( (R - c(2,1,-1,3)) == 0)) {
		return(c(2, -22, 85, -127, 63));
	}
	# R = c(2,1,1,3)
	if(all( (R - c(2,1,1,3)) == 0)) {
		return(c(2, -22, 77, -91, 35));
	}
	# R = c(2,1,1,-3)
	if(all( (R - c(2,1,1,-3)) == 0)) {
		return(c(2, 26, 101, 125, 35));
	}
	# R = c(-2,1,1,-3)
	if(all( (R - c(-2,1,1,-3)) == 0)) {
		return(c(2, 26, 109, 161, 63));
	}
	# R = c(2,-1,1,-3)
	if(all( (R - c(2,-1,1,-3)) == 0)) {
		return(c(2, 26, 113, 201, 99));
	}
	# R = c(0,1,-2,-3)
	if(all( (R - c(0,1,-2,-3)) == 0)) {
		return(c(2, 24, 93, 130, 49));
	}
	# R = c(0,-1,-2,-3)
	if(all( (R - c(0,-1,-2,-3)) == 0)) {
		return(c(2, 24, 105, 202, 121));
	}
	# R = c(0,0,1,4)
	if(all( (R - c(0,0,1,4)) == 0)) {
		return(c(2, -32, 176, -383, 256));
	}
	stop("TODO");
}
solve.S4C2.SymVar1.P2P1sp = function(ps, R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	S  = R1;
	# x0
	x0_coeff = c(- 192, 176*S^2 + 512*R4,
		- 68*S^4 - 416*S^2*R4 - 480*R4^2 + 320*R2,
		32*R3^2 + 13*S^6 + 144*S^4*R4 + 320*S^2*R4^2 + 192*R4^3 - 240*S^2*R2 - 384*R4*R2,
		- 16*R3^2*S^2 - S^8 - 26*S^6*R4 - 86*S^4*R4^2 - 96*S^2*R4^3 - 32*R4^4 + 76*S^4*R2 +
			+ 192*S^2*R4*R2 + 128*R4^2*R2 - 128*R2^2,
		2*R3^2*S^4 + 2*S^8*R4 + 13*S^6*R4^2 + 12*S^4*R4^3 + 16*S^2*R4^4 - 13*S^6*R2 +
			- 24*S^4*R4*R2 - 64*S^2*R4^2*R2 + 64*S^2*R2^2,
		- S^8*(R4^2 - R2) - 2*S^4*(R4^2 - 2*R2)^2);
	x0 = sapply(ps, function(ps) {
		sum(x0_coeff * ps^seq(6, 0, by=-1));
	});
	# div
	div_coeff = c(128, - 128*S^2 - 192*R4,
		32*R3*S + 40*S^4 + 192*S^2*R4 + 64*R4^2 - 128*R2,
		- 16*R3*S^3 - 4*S^6 - 60*S^4*R4 - 64*S^2*R4^2 + 128*S^2*R2,
		2*R3*S^5 + 6*S^6*R4 + 20*S^4*R4^2 - 40*S^4*R2,
		- 2*S^6*R4^2 + 4*S^6*R2);
	div = sapply(ps, function(ps) {
		sum(div_coeff * ps^seq(5, 0, by=-1));
	});
	ps = x0 / div;
	return(ps);
}

### Examples:

###
R = c(4,1,1,3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 2:
R = c(2,1,-1,3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 3:
R = c(2,1,1,3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 4:
R = c(2,1,1,-3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 5:
R = c(-2,1,1,-3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 6:
R = c(2,-1,1,-3)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")


### Ex 7:
R = c(0,1,-2,-3)
sol = solve.S4C2.SymVar1.P2P1(R)

round0(test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2"))


### Ex 8:
R = c(0,-1,-2,-3)
sol = solve.S4C2.SymVar1.P2P1(R)

round0(test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2"))


### Ex 9:
R = c(4,-1,-3,5)
sol = solve.S4C2.SymVar1.P2P1(R)

test.S4C2.Var(sol, n=c(1,2,2,1,1), type="x1y2")



##############
### Derivation

### Transform 1:

### (Eq 4)^2 =>
(x1*y2)^2 + (x2*y1)^2 + 2*E4 - R4^2 # = 0
# Eq 2 + (Eq 4)^2 =>
(x1^2 + x2^2)*(y1^2 + y2^2) + 2*E4 - R2 - R4^2 # = 0
(s1^2 - 2*p1)*(s2^2 - 2*p2) + 2*E4 - R2 - R4^2 # = 0
### Eq 2-bis:
s1^2*s2^2 - 2*(p1*s2^2 + p2*s1^2) + 6*E4 - R2 - R4^2 # = 0

###
A = x1*y1 + x2*y2;
B = x1*y2 + x2*y1; # = R[4]
#
A*B - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
A + B - s1*s2 # = 0
# =>
R4*(s1*s2 - R4) - p1*(s2^2 - 2*p2) - p2*(s1^2 - 2*p1) # = 0
p1*(s2^2 - 2*p2) + p2*(s1^2 - 2*p1) - R4*s1*s2 + R4^2 # = 0
### Eq 4-bis:
p1*s2^2+ p2*s1^2 - 4*E4 - R4*s1*s2 + R4^2 # = 0

### Eq 2-bis & Eq 4-bis =>
s1^2*s2^2 - 2*R4*s1*s2 - 2*E4 - R2 + R4^2 # = 0

### Transformed System [1]:
s1 + s2 - R1 # = 0
s1^2*s2^2 - 2*R4*s1*s2 - 2*p1*p2 - R2 + R4^2 # = 0
s1*p2 + s2*p1 - R3 # = 0
p1*s2^2+ p2*s1^2 - 4*p1*p2 - R4*s1*s2 + R4^2 # = 0

### Transform 2:

### Eq 3:
A2 = s1*p2 + s2*p1; # = R3
B2 = s1*p1 + s2*p2;
# =>
A2 + B2 - (p1 + p2)*(s1 + s2) # = 0
A2 * B2 - s1*s2*(p1^2+p2^2) - p1*p2*(s1^2+s2^2) # = 0
# =>
A2 + B2 - sp*S # = 0
A2 * B2 - ps*(sp^2 - 2*E4) - E4*(S^2 - 2*ps) # = 0
# =>
R3*(sp*S - R3) - ps*(sp^2 - 2*E4) - E4*(S^2 - 2*ps) # = 0
ps*(sp^2 - 2*E4) + E4*(S^2 - 2*ps) - R3*sp*S + R3^2 # = 0
ps*sp^2 + E4*S^2 - 4*E4*ps - R3*sp*S + R3^2 # = 0

### [for Eq 4]
A2 = p2*s2^2+ p1*s1^2;
B2 = p1*s2^2+ p2*s1^2;
# =>
A2 + B2 - (p1+p2)*(s1^2 + s2^2) # = 0
A2 * B2 - p1*p2*(s1^4 + s2^4) - (s1*s2)^2*(p1^2 + p2^2) # = 0
# =>
A2 + B2 - sp*(S^2 - 2*ps) # = 0
A2 * B2 - E4*(S^4 - 4*ps*S^2 + 2*ps^2) - ps^2*(sp^2 - 2*E4) # = 0
# =>
B2*(B2 - sp*(S^2 - 2*ps)) + E4*(S^4 - 4*ps*S^2 + 2*ps^2) + ps^2*(sp^2 - 2*E4) # = 0
B2^2 - B2*sp*(S^2 - 2*ps) + E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0
# B2 = 4*E4 + R4*ps - R4^2 # = 0
# =>
(4*E4 + R4*ps - R4^2)^2 - (4*E4 + R4*ps - R4^2)*sp*(S^2 - 2*ps) +
	+ E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0


### Transformed System [2]:
S - R1 # = 0
ps^2 - 2*R4*ps - 2*E4 - R2 + R4^2 # = 0
ps*sp^2 + E4*S^2 - 4*E4*ps - R3*sp*S + R3^2 # = 0
(4*E4 + R4*ps - R4^2)^2 - (4*E4 + R4*ps - R4^2)*sp*(S^2 - 2*ps) +
	+ E4*(S^4 - 4*ps*S^2) + ps^2*sp^2 # = 0

# TODO: solve;

### Debug
# id = 1; x1 = sol[id,1]; x2 = sol[id,2]; y1 = sol[id,3]; y2 = sol[id,4];
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;
sp = p1 + p2; ps = s1 * s2;
S  = s1 + s2; # = R1;
E4 = p1 * p2;

###
library(nleqslv)

# NO complex roots!
psys = function(x) {
	y = numeric(4);
	y[1] = sum(x) - R[1];
	y[2] = (x[1]*x[3])^2 + (x[2]*x[4])^2 - R[2];
	y[3] = x[1]*x[2]*(x[3] + x[4]) + x[3]*x[4]*(x[1] + x[2]) - R[3];
	y[4] = x[1]*x[4] + x[2]*x[3] - R[4];
	return(y);
}
jac = function(x) {
	m = matrix(1, 4, 4);
	f3 = function(x) x[1]*x[2] + x[1]*x[3] + x[2]*x[3];
	m[2,] = 2*c(x[1]*x[3]^2, x[1]^2*x[2], x[3]*x[4]^2, x[3]^2*x[4]);
	m[3,] = c(f3(x[-1]), f3(x[-2]), f3(x[-3]), f3(x[-4]));
	m[4,] = x[c(4, 3, 2, 1)];
	return(m)
}

R  = c(4,1,1,3)
x0 = c(1,1/2,-1/3,2);
sol = nleqslv(x0, psys, jac=jac, method="Newton")
sol = matrix(sol$x, nrow=1);

###
p1 = toPoly.pm("ps^2 - 2*R4*ps - 2*E4 - R2 + R4^2")
p2 = toPoly.pm("ps*sp^2 + E4*S^2 - 4*E4*ps - R3*sp*S + R3^2")
p3 = toPoly.pm("(4*E4 + R4*ps - R4^2)^2 - (4*E4 + R4*ps - R4^2)*sp*(S^2 - 2*ps) +
	+ E4*(S^4 - 4*ps*S^2) + ps^2*sp^2")
#
pR = solve.lpm(p1, p2, p3, xn = c("E4", "sp"))
str(pR) # 392 monomes
max(pR[[2]]$Rez$ps)

p0 = pR[[2]]$Rez[pR[[2]]$Rez$ps == 0, ]
p0$coeff = p0$coeff / 2;
p0 = drop.pm(p0)
p0 = div.pm(p0, toPoly.pm("(R4^2 - 2*R2)*(R4^2 - 2*R2 - R3*S)"), "R4")$Rez
p0 = sort.pm(p0, "R4")

pR = pR[[2]]$Rez;

R = c(2,1,-1,3)
cp = tapply(seq(nrow(pR)), pR$ps, function(id) eval.pm(pR[id,], list(S=R[1], R2=R[2], R3=R[3], R4=R[4], ps=1)))
paste0(rev(cp) / 1024, ", ", collapse="")


### R3 = 0
g1 = toPoly.pm("ps^2 - 2*R4*ps - 2*E4 - R2 + R4^2")
g2 = toPoly.pm("(4*E4 + R4*ps - R4^2)^2 + ps*E4*(S^2 - 4*ps)")

gR = solve.pm(g1, g2, xn=c("E4"))
str(gR)
gR = gR$Rez
gR = sort.pm(gR, "ps")
print.pm(gR, lead="ps")


####################
####################

### Power: Both Eqs

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
	# Special cases
	sol0 = solve.S4C2.Var_x1y1P2Special(R, debug=debug, all=all);
	if(sol0$isSpecial) {
		return(sol0$sol);
	}
	### [old] based on P[16]
	# coeff = coeff.S4C2.Var_x1y1P2(R);
	# s1 = roots(coeff);
	### P[8]
	R24 = R[2] + R[4]; R24dsq = (R[2] - R[4])^2;
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(1, 2*R1^2, - 4*(3*R1*R3 + 2*R24),
		- 4*(R1^2*R24 - 4*R3^2), 2*(8*R1*R3*R24 + 7*R24^2 + 4*R2*R4),
		2*R1^2*R24dsq,
		- 4*(R1*R3*R24dsq + 2*R24^3 - 8*R2*R4*R24), 0, R24dsq^2);
	ps  = roots(coeff);
	len = length(ps);
	s1  = sapply(seq(len), function(id) {
		roots(c(1, -R[1], ps[id]));
	})
	s1 = as.vector(s1);
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
### [old] P[16]
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
solve.S4C2.Var_x1y1P2Special = function(R, debug=TRUE, all=FALSE) {
	if(R[2] != R[4]) return(list(isSpecial = FALSE));
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	solve2 = function() {
		len = length(x1);
		y12 = sapply(seq(len), function(id) {
			roots(c(1, -s2[id], p2[id]));
		});
		y12 = t(y12);
		sol = cbind(x1, x2, y12);
		# NO rbind( , y2, y1): overlaps with all = TRUE;
		sol = rbind(sol, sol[ , c(3,4,1,2)]);
		return(sol)
	}
	### Case: x1 == x2
	# 2*x1^4 - 3*R1*x1^3 + R1^2*x1^2 - R3*x1 - R2
	coeff = c(2, - 3*R1, R1^2, - R3, - R2);
	x1 = roots(coeff);
	x2 = x1;
	if(debug) print(x1);
	s2 = R1 - 2*x1;
	p2 = (s2^2 - R2/x1^2) / 2;
	sol = solve2();
	# sol = rbind(sol, sol[ , c(3,4,1,2)]);
	### Case: x1 == -x2
	# s1 = 0;
	x1 = rootn( - R3 / R1, 2);
	x1 = c(x1, -x1);
	if(debug) print(x1);
	x2 = - x1;
	p2 = (R1^2 - R2/x1^2) / 2;
	s2 = rep(R1, length(x1));
	sol2 = solve2();
	#
	sol  = rbind(sol, sol2);
	#
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(list(isSpecial=TRUE, sol=sol));
}

# TODO:
# Special Cases:
# R2 = R4;
# s2 = 0; s1 = R1;
# s1 = 0; s2 = R1;

### Examples:

R = c(2,3,-1,5)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")

# E4 actually satisfies a P[8]
# - implemented using ps = s1*s2: satisfies also a P[8];
E4 = apply(sol, 1, prod);
x  = 2 * E4;
26269 - 52740*x + 37652*x^2 - 12296*x^3 + 2670*x^4 - 884*x^5 + 260*x^6 - 32*x^7 + x^8


### Ex 2:
R = c(2,3,-1,-1)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")


### Ex 3: Special
R = c(2,3,-1,3)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")


### Ex 4: Special
R = c(-1,3,2,3)
sol = solve.S4C2.Var_x1y1P2(R)

test.S4C2.Var(sol, n=c(1,2,2,2,2), type="x1y2")


##########

### Debug:
R = c(2,3,-1,5)
x1 =  1.1725398503 + -1.4433080165i;
x2 = -0.7162047206 + -0.2205315378i;
y1 =  0.7158321338 + 0.7331960503i;
y2 =  0.8278327365 + 0.9306435040i;
sol = cbind(x1,x2,y1,y2);

# x1 = sol[1,1]; x2 = sol[1,2]; y1 = sol[1,3]; y2 = sol[1,4];
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;
sp = p1 + p2; ps = s1 * s2;
S  = s1 + s2; # = R1;
E4 = p1 * p2;


### Derivation
p1 = toPoly.pm("s1 + s2 - R1")
p2 = toPoly.pm("s1*p2 + s2*p1 - R3")
p3 = toPoly.pm("(s1^2 - 2*p1)*(s2^2 - 2*p2) - R2 - R4")
p4 = toPoly.pm("p2^2*(s1^4 - 4*p1*s1^2) + p1^2*(s2^4 - 4*p2*s2^2) + 4*p1^2*p2^2 - R2*R4")

pR = solve.lpm(p1, p2, p3, p4, xn=c("s2", "p2", "p1"))
sort(unique(pR[[3]]$Rez$s1))

pD = div.pm(pR[[3]]$Rez, toPoly.pm("(s1 - R1)^7"), "s1")
p = toCoeff(pD$Rez, xn="s1", print=TRUE)

###
source("Polynomials.Helper.BigNumbers.R")
pD = pD$Rez
pD$coeff = pD$coeff * 2^16
pD = replace.pm(pD, data.frame(s1=c(1,0), R1=c(0,1), coeff=c(1,1/2)), xn="s1")
eval.pm(pD, list(R1=R1, R2=R2, R3=R3, R4=R4, s1 = s1 - R1/2))
#
pD$coeff = as.bigz(pD$coeff)
pD = square.pm(pD, xn="s1")
# eval.pm(pD, list(R1=R1, R2=R2, R3=R3, R4=R4, s1 = (s1 - R1/2)^2))
pX = data.frame(s1=c(1,0), R1=c(0,2), coeff=c(4,1))
pX$coeff = as.bigz(pX$coeff) / 4;
# TODO: fix bug with bigq;
pD = replace.pm(pD, pX, xn="s1")

### P[8] in ps
pSq = toPoly.pm("ps^8 + 2*R1^2*ps^7 - 4*(3*R1*R3 + 2*R2 + 2*R4)*ps^6 +
	- 4*(R1^2*(R2 + R4) - 4*R3^2)*ps^5 + 2*(8*R1*R3*(R2+R4) + 7*R2^2 + 7*R4^2 + 18*R2*R4)*ps^4 +
	+ 2*R1^2*(R2 - R4)^2*ps^3 +
	- 4*(R1*R3*(R2 - R4)^2 + 2*(R2^3 + R4^3) - 2*R2*R4*(R2 + R4))*ps^2 + (R2-R4)^4")

eval.pm(pSq, list(ps=ps[1], R1=R1, R2=R2, R3=R3, R4=R4))

# Test
pSq$coeff = as.bigz(pSq$coeff)
mult.pm(pSq, pSq)



### Derivation E4: P[8]
# - overflows!
# p2 = replace.fr.pm(p2, data.frame(E4=1, coeff=1), data.frame(p1=1, coeff=1), xn="p2")
# ...

### Eq 1:
E4*(S^2 - 2*ps) + ps*(sp^2 - 2*E4) - R3*S*sp + R3^2 # = 0
### Eq 2:
ps^2 - 2*R3*S + 2*ps*sp + 4*E4 - R2 - R4 # = 0
### Eq 3:
# based on:
p2^2*(s1^4 - 4*p1*s1^2) + p1^2*(s2^4 - 4*p2*s2^2) + 4*p1^2*p2^2 - R2*R4 # = 0
p2^2*s1^4 + p1^2*s2^4 - 4*E4*(p2*s1^2 + p1*s2^2) + 4*E4^2 - R2*R4 # = 0
# p2*s1^2 + p1*s2^2 = R3*S - ps*sp =>
(R3*S - ps*sp)^2 - 2*E4*ps^2 - 4*E4*(R3*S - ps*sp) + 4*E4^2 - R2*R4 # = 0


p1 = toPoly.pm("E4*(S^2 - 2*ps) + ps*(sp^2 - 2*E4) - R3*S*sp + R3^2")
p2 = toPoly.pm("ps^2 - 2*R3*S + 2*ps*sp + 4*E4 - R2 - R4")
p3 = toPoly.pm("(R3*S - ps*sp)^2 - 2*E4*ps^2 - 4*E4*(R3*S - ps*sp) + 4*E4^2 - R2*R4")

# TODO: still huge!
pR = solve.lpm(p1, p2, p3, xn = c("ps", "sp"))
str(pR)

