########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
###
### draft v.0.1i


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


# - moved to file:
#   Poly.System.S4..C2.Helper.R;


####################
####################

##############
### Theory ###
##############

### System:
# x1^n1 + x2^n1 + y1^n1 + y2^n1 = R1
# x1^n2*y1^n3 + x2^n2*y2^n3 = R2
# x1*x2*y1 + x1*x2*y2 + x1*y1*y2 + x2*y1*y2 = R3
# x1*x2*y1*y2 = R4

### Symmetries:
# - if (x1, x2, y1, y2) is a solution,
#   so is also the C2 permutation: (x2, x1, y2, y1);

##################
### Basic Type ###
##################

### Symmetric Eq 2
# n2 = n3 = n
# (x1*y1)^n + (x2*y2)^n = R2
# Note:
# - has additional symmetry;
# - if (x1, x2, y1, y2) is a solution,
#   so is also the permutation: (y1, y2, x1, x2) & its C2 permutation;

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


###############
###############

###############
### Ht Type ###
###############

### Non-Symmetric Eq 2
# x1^n2*y1^n3 + x2^n2*y2^n3 = R2

### Case: n2 = 2
# x1^2*y1 + x2^2*y2 = R2

### Solution:

# Translate system to new variables:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;

# let:
A = x1^2*y1 + x2^2*y2; # = R2
B = x1^2*y2 + x2^2*y1;
# =>
A + B - s2*(s1^2 - 2*p1) # = 0
A * B - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0
# A = R2 =>
R2*(s2*(s1^2 - 2*p1) - R2) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0

### Transformed System:
s1 + s2 - R1 # = 0
R2*s2*(s1^2 - 2*p1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 - R2^2 # = 0
p2*s1 + p1*s2 - R3 # = 0
p1*p2 - R4 # = 0


### Solver:
solve.S4C2.HtP21 = function(R, debug=TRUE, all=FALSE) {
	### Special cases:
	# check/solve when s1 = 0 or s2 = 0;
	tol = 1E-12;
	sol0 = solve.S4C2.HtP21Cases(R, debug=debug, all=all, tol=tol);
	#
	coeff = coeff.S4C2.HtP21(R);
	s1 = roots(coeff);
	if(sol0$isSpecial) {
		if(sol0$type == 1) {
			isZ = abs(s1) < tol;
			s1  = s1[ ! isZ];
		}
	}
	s2 = R[1] - s1;
	if(debug) print(s1);
	p1 = solveP.S4C2.HtP21(s1, R);
	p2 = (R[3] - p1*s2) / s1;
	# Special Case: s1 == 0 is handled above;
	# Step 3:
	solve2 = function(id, s, p) {
		r = roots(c(1, -s[id], p[id]));
	}
	len = length(s1);
	x12 = lapply(seq(len), solve2, s1, p1);
	x12 = matrix(unlist(x12), nrow=2);
	x12 = t(x12);
	# y12 = lapply(seq(len), solve2, s2, p2);
	# y12 = matrix(unlist(y12), nrow=2);
	# robust:
	y1 = (x12[,2]^2*s2 - R[2]) / (s1 * (x12[,2] - x12[,1]));
	y2 = s2 - y1;
	sol = cbind(x12, y1, y2);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	# add back the correct root;
	if(sol0$isSpecial) sol = rbind(sol, sol0$sol);
	colnames(sol) = c("x1", "x2", "y1", "y2");
	return(sol);
}
solve.S4C2.HtP21Cases = function(R, debug=TRUE, all=FALSE, tol=1E-12) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	### Case: s1 == 0
	# R2 + R3 = 0
	if(abs(R2 + R3) < tol) {
		x1 = rootn(R2 / R1, n=2);
		x2 = - x1;
		if(debug) print(paste("Special: s1 = 0; x1 = ", x1));
		# Case: R2 = R3 = 0, R4 != 0;
		if(R3 == 0) {
			if(R4 != 0) {
				# NO solution
				sol = cbind(x1=NA, x2=NA, y1=NA, y2=NA);
				warning("No solutions!");
				return(list(isSpecial=TRUE, sol=sol, type=1));
			} else {
				sol = cbind(x1=0, x2=0, y1=R1, y2=0);
				warning("Indeterminate solution!");
				return(list(isSpecial=TRUE, sol=sol, type=1));
			}
		}
		p2 = R1 * R4 / R3;
		y1 = roots(c(1, -R1, p2));
		y2 = R1 - y1;
		sol = cbind(x1, x2, y1, y2);
		if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
		return(list(isSpecial = TRUE, sol=sol, type=1));
	}
	### Case: s2 == 0
	if(abs(R1^2*(4*R4 - R1*R3) - R2^2) < tol) {
		y1 = rootn( - R3/R1, n=2);
		y2 = - y1;
		if(debug) print(paste("Special: s2 = 0; y1 = ", y1));
		dx = R2 / (R1 * y1);
		x1 = (R1 + dx) / 2;
		x2 = R1 - x1;
		sol = cbind(x1, x2, y1, y2);
		if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
		return(list(isSpecial = TRUE, sol=sol, type=2));
	}
	return(list(isSpecial=FALSE));
}
coeff.S4C2.HtP21 = function(R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(R4, - R1*R4, - (R2^2 + R2*R3),
		(R3*R1*R2 + 2*R1*R2^2 + R3*R4 + 2*R4*R2),
		(- R1^2*R2^2 + 3*R3*R1*R4 - 9*R4^2),
		(- 2*R1^2*R4*R2 - R3^3 - 6*R1*R4^2 - 3*R3^2*R2 - 3*R3*R2^2 - 2*R2^3),
		(2*R1*R2^3 - R1^2*R4^2 + R3^2*R1*R2 + 2*R3*R1*R2^2 + 4*R3^2*R4 + 10*R3*R4*R2 + 10*R4*R2^2),
		(- 2*R1*R4*R2^2 - 2*R3*R1*R4*R2),
		- R2^2 * (R2^2 + 2*R3*R2 + R3^2));
	return(coeff);
}
solveP.S4C2.HtP21 = function(s1, R) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	p1 = - R4*s1^6 + 2*R1*R4*s1^5 - R1^2*R4*s1^4 + 2*R2*R4*s1^3 + R4*R3*s1^3 +
		- 4*R2*R1*R4*s1^2 - 2*R1*R4*R3*s1^2 + 2*R2*R1^2*R4*s1 + R1^2*R4*R3*s1;
	div = R2*s1^5 - 3*R2*R1*s1^4 - 3*R4*s1^4 + 3*R2*R1^2*s1^3 + 5*R1*R4*s1^3 + R2^2*s1^2 +
		- R2*R1^3*s1^2 - R1^2*R4*s1^2 + 2*R2*R3*s1^2 + R3^2*s1^2 - 2*R2^2*R1*s1 - R1^3*R4*s1 +
		- 4*R2*R1*R3*s1 - 2*R1*R3^2*s1 + R2^2*R1^2 + 2*R2*R1^2*R3 + R1^2*R3^2;
	p1 = p1 / div;
	return(p1);
}

### Examples

### Ex 1:
R = c(2,3,-2,5)
sol = solve.S4C2.HtP21(R)
x1 = sol[,1]; x2 = sol[,2]; y1 = sol[,3]; y2 = sol[,4];

test.S4C2(sol, n=c(1,2,1))


### Ex 2:
R = c(-1,3,2,5)
sol = solve.S4C2.HtP21(R)
x1 = sol[,1]; x2 = sol[,2]; y1 = sol[,3]; y2 = sol[,4];

test.S4C2(sol, n=c(1,2,1))


### Ex 3: Special
R = c(3,1,-1,3)
sol = solve.S4C2.HtP21(R)

test.S4C2(sol, n=c(1,2,1))


### Ex 4: Special
# TODO: remove false root;
R = c(-1,3,5,1)
sol = solve.S4C2.HtP21(R)

test.S4C2(sol, n=c(1,2,1))


###################
### Derivation: A*B
A*B - p2*(x1^4 + x2^4) - p1^2*(s2^2 - 2*p2) # = 0
A*B - p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) - p1^2*(s2^2 - 2*p2) # = 0
A*B - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 # = 0

### Solution:

### Eq 3:
p2*s1 + p1*s2 - R3 # = 0
p1*p2*s1 + p1^2*s2 - R3*p1 # = 0
R4*s1 + p1^2*(R1 - s1) - R3*p1 # = 0

### Eq 2:
R2*s2*(s1^2 - 2*p1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 - R2^2 # = 0
R2*s2*p1*(s1^2 - 2*p1) - R4*(s1^4 - 4*p1*s1^2) - p1^3*s2^2 - R2^2*p1 # = 0

### Debug
R = c(2,3,-2,5)
x1 =  1.0279315562 + 0.0000000000i;
x2 = -1.1340431142 + 0.0000000000i;
y1 = -1.2703309529 + -0.0000000000i;
y2 =  3.3764425109 + -0.0000000000i;

###
p1 = toPoly.pm("R4*s1 + p1^2*(R1 - s1) - R3*p1")
p2 = toPoly.pm("R2*(R1 - s1)*p1*(s1^2 - 2*p1) - R4*(s1^4 - 4*p1*s1^2) - p1^3*(R1 - s1)^2 - R2^2*p1")
pR = solve.pm(p1, p2, "p1")
pMAX = max(pR$Rez$s1)
pR$Rez$coeff = sign(pR$Rez$coeff[pR$Rez$s1 == pMAX]) * pR$Rez$coeff;
#
pR$Rez = div.pm(pR$Rez, toPoly.pm("(s1-R1)^4"), "s1")$Rez
pR$Rez = sort.pm(pR$Rez, "s1")
str(pR$Rez)

p = sapply(unique(pR$Rez$s1), function(id) {
	p = pR$Rez[pR$Rez$s1 == id, , drop=FALSE];
	p$s1 = NULL;
	# paste0("(", as.character.pm(p), ")*", "s1", "^", id, " +");
	paste0("(", as.character.pm(p), "),");
})
cat(p, sep = rep("\n", length(p)))


######################

######################
### Eq 4: Variants ###
######################

#################
### Variant Eq 4:
x1*x2 + y1*y2 - R4 # = 0

### Transformed System:
s1 + s2 - R1 # = 0
R2*s2*(s1^2 - 2*p1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 - R2^2 # = 0
p2*s1 + p1*s2 - R3 # = 0
p1 + p2 - R4 # = 0

# TODO: solve;


#################
### Variant Eq 4:
x1*y1 + x2*y2 - R4 # = 0

### Transformed System:
s1 + s2 - R1 # = 0
R2*s2*(s1^2 - 2*p1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 - R2^2 # = 0
p2*s1 + p1*s2 - R3 # = 0
### Eq 4:
# - simple version: works only with Ht21;
R4*s1 - p1*s2 - R2 # = 0
# - complex version:
R4*s1*s2 - p2*s1^2 - p1*s2^2 + 4*p1*p2 - R4^2 # = 0

### Derivation:
A = R4;
B = x1*y2 + x2*y1;
A + B - s1*s2 # = 0;
A * B - p2*(s1^2 - 2*p1) - p1*(s2^2 - 2*p2) # = 0

# TODO: solve;



###########################

###########################
### Eq 1: Higher Powers ###
###########################


### Eq 1: Power 2 (part)
# - simplified version: only x^2;

### System:
# x1^2 + x2^2 + y1 + y2 = R1
# x1^n2*y1^n3 + x2^n2*y2^n3 = R2
# x1*x2*y1 + x1*x2*y2 + x1*y1*y2 + x2*y1*y2 = R3
# x1*x2*y1*y2 = R4


### Solution:

# Translate system to new variables:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;

### Transformed System:
s1^2 - 2*p1 + s2 - R1 # = 0
R2*s2*(s1^2 - 2*p1) - p2*(s1^4 - 4*p1*s1^2) - p1^2*s2^2 - R2^2 # = 0
p2*s1 + p1*s2 - R3 # = 0
p1*p2 - R4 # = 0

# TODO: solve;

