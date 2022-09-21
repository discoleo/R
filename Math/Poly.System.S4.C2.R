########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
###
### draft v.0.2e


####################

### Helper Functions


source("Poly.System.S4.C2.Helper.R")


# - moved to file:
#   Poly.System.S4.C2.Helper.R;


test.S4C2.E2aP1 = function(sol, R = NULL, n=c(2,1), type=0) {
	x1 = sol[,1]; x2 = sol[,2];
	y1 = sol[,3]; y2 = sol[,4];
	err1 = x1 + x2; err2 = y1 + y2;
	err3 = x1^n[1]*y1^n[2] + x2^n[1]*y2^n[2];
	err4 = if(type == 0) { x1*x2*y1*y2; }
		else if(type == 1) { x1*x2 + y1*y2; }
		else NA;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

# Generate Classic Polynomial
classicPoly.E2a = function(nx, ny, print=TRUE) {
	pX = toPoly.pm("s1 - x1");
	pY = toPoly.pm("s2 - y1");
	pE = toPoly.pm("x1^nx*y1^ny + x2^nx*y2^ny - R3");
	pE = replace.pm(pE, pX, "x2");
	pE = replace.pm(pE, pY, "y2");
	pE4 = toPoly.pm("x1*x2*y1*y2 - E4")
	pE4 = replace.pm(pE4, pX, "x2");
	pE4 = replace.pm(pE4, pY, "y2");
	
	pR = solve.pm(pE4, pE, "y1");
	pR = pR$Rez;
	pR = sort.pm(pR, "x1", c("s1", "s2", "E4"));
	pR$coeff = - pR$coeff;
	if(print) print.pm(pR, lead="x1");
	invisible(pR);
}
classicPoly.E2ax2 = function(nx, ny, print=FALSE) {
	pX = toPoly.pm("s1 - x1");
	pY = toPoly.pm("s2 - y1");
	pE1 = toPoly.pm("x1^nx*y1^ny + x2^nx*y2^ny - R3");
	pE1 = replace.pm(pE1, pX, "x2");
	pE1 = replace.pm(pE1, pY, "y2");
	pE2 = toPoly.pm("x1^ny*y1^nx + x2^ny*y2^nx - R4");
	pE2 = replace.pm(pE2, pX, "x2");
	pE2 = replace.pm(pE2, pY, "y2");
	
	pR = solve.pm(pE1, pE2, "y1");
	pR = pR$Rez;
	pR = orderVars.pm(pR, c("R3", "R4"), last=FALSE);
	pR = sort.pm(pR, "x1", c("s1", "s2"));
	# pR$coeff = - pR$coeff;
	if(print) print.pm(pR, lead="x1");
	invisible(pR);
}


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


### Note:
# - the systems with additional symmetry moved to file:
#   Poly.System.S4.C2.Symmetric.R;


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


############################
############################

#####################
### Separate Vars ###
#####################

# - Eqs. 1 & 2: Separate variables;
# - Eqs. 3 & 4: Entangled variables;

### System:
# x1^n1 + x2^n1 = R1
# y1^n2 + y2^n2 = R2
# x1^n3 * y1^n4 + x2^n3 * y2^n4 = R3
# x1*x2*y1*y2 = R4

### Symmetries:
# - if (x1, x2, y1, y2) is a solution,
#   so is also the C2 permutation: (x2, x1, y2, y1);

### Case:
# n1 = n2 = 1;
# n3 = 2; n4 = 1;

### Solution:

### Transform:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;

A = x1^2*y1 + x2^2*y2;
B = x1^2*y2 + x2^2*y1;

### A + B =>
A + B - (y1 + y2)*(x1^2 + x2^2) # = 0
A + B - s2*(s1^2 - 2*p1) # = 0
A + B + 2*s2*p1 - s1^2*s2 # = 0

### A * B =>
A * B - p1^2*(y1^2 + y2^2) - p2*(x1^4 + x2^4) # = 0
A * B - p1^2*(s2^2 - 2*p2) - p2*(s1^4 - 4*p1*s1^2 + 2*p1^2) # = 0
A * B - p1^2*s2^2 - p2*(s1^4 - 4*p1*s1^2) # = 0
A * B - p1^2*s2^2 - p2*s1^4 + 4*s1^2*E4 # = 0

# =>
A^2 + (2*s2*p1 - s1^2*s2)*A + p1^2*s2^2 + p2*s1^4 - 4*s1^2*E4 # = 0

### Transformed System:
s1 - R1 # = 0
s2 - R2 # = 0
A - R3 # = 0
A^2 + (2*s2*p1 - s1^2*s2)*A + p1^2*s2^2 + p2*s1^4 - 4*s1^2*E4 # = 0
p1*p2 - R4 # = 0

### Solution Tr. System:
s2^2*p1^3 + 2*s2*A*p1^2 + (A^2 - 4*s1^2*E4 - s1^2*s2*A)*p1 + E4*s1^4 # = 0

### Debug:
x1 = x[1]; x2 = x[3]; y1 = x[2]; y2 = x[4];
R1 = x1 + x2; R2 = y1 + y2;
R3 = x1^2*y1 + x2^2*y2;
R4 = x1*x2*y1*y2;

### Solver:

solve.S4C2.E21a = function(R, debug=TRUE, all=FALSE) {
	s1 = R[1]; s2 = R[2]; E21a = R[3]; E4 = R[4];
	# Special Cases:
	isSpecial = FALSE;
	if(round0(s1^2*s2 - 4*E21a) == 0) {
		isSpecial = TRUE;
		warning("Special case!");
		x = s1/2;
		y = roots(c(1, -s2, E4 / x^2));
		sol0 = cbind(x1=x, x2=x, y1=y[1], y2=y[2]);
		# 2nd Set:
		if(round0(s1^2 - 16*E4) == 0) {
			# another set: x1 == x3;
			coeff = c(1, s1^2);
		} else {
			coeff = c(4, 3*s1^2, -16*s1^2*E4);
		}
		# Note: includes also the permutation;
		# - but useful for generating the Classic polynomial;
	} else {
		coeff = c(s2^2, 2*s2*E21a, (E21a^2 - s1^2*s2*E21a - 4*s1^2*E4), E4*s1^4);
	}
	#
	p1 = roots(coeff);
	if(debug) print(p1);
	# Step 2:
	x1 = sapply(p1, function(p) roots(c(1, -s1, p)));
	x1 = as.vector(x1);
	x2 = s1 - x1;
	# Step 3: robust
	# p1, p2 = NOT needed;
	y1 = (x2^2*s2 - E21a) / (x2^2 - x1^2);
	y2 = s2 - y1;
	#
	sol = cbind(x1, x2, y1, y2);
	if(isSpecial) sol = rbind(sol, sol0);
	if(all && ! isSpecial) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
test.S4C2.E21a = function(sol, R=NULL) {
	test.S4C2.E2aP1(sol=sol, R=R, n=c(2,1));
}

### Examples:

### Ex 1:
R = c(-2,-3,5,-1)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)

poly.calc(sol[,1]) * 9


### Ex 2:
R = c(1,-2,3,-1)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)


### Ex 3: b4 = 0
s1 = 2;
R = c(s1,1,3*s1^2/2, -2)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)

poly.calc(sol[,1])

x = sol[,1];
32 - 88*x - 4*x^2 + 40*x^3 - 6*x^5 + x^6


### Ex 4: b4 = 0
s1 = 2;
R = c(s1,1,3*s1^2/2, 1/4)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)

x = sol[,1];
-4 - 16*x - 40*x^2 + 40*x^3 - 6*x^5 + x^6


### Ex 5: Special Case
R = c(-2,1,1, 1/4)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)


### Ex 6: Special Case
R = c(-2,1,1, 3/4)
sol = solve.S4C2.E21a(R)

test.S4C2.E21a(sol)


#################
### Classic Poly:

s1 = R[1]; s2 = R[2]; R3 = R[3]; E4 = R[4];
x1 = sol[,1];

s2^2*x1^6 - 3*s1*s2^2*x1^5 + (3*s1^2*s2^2 - 2*R3*s2)*x1^4 +
	- s1*s2*(s1^2*s2 - 4*R3)*x1^3 - (3*R3*s1^2*s2 + 4*s1^2*E4 - R3^2)*x1^2 +
	+ (R3*s1^3*s2 + 4*s1^3*E4 - R3^2*s1)*x1 - s1^4*E4 # = 0

### Derivation:

nx = 2; ny = 1;
pR = classicPoly.E2a(nx, ny);
pR = sort.pm(pR, "x1", xn2 = c("s1", "s2", "E4"))
str(pR);


############
### Variant:

### System:
# x1 + x2 = R1
# y1 + y2 = R2
# x1^2*y1 + x2^2*y2 = R3
# x1*x2 + y1*y2 = R4

### Transformed System:
s1 - R1 # = 0
s2 - R2 # = 0
A - R3 # = 0
A^2 + (2*s2*p1 - s1^2*s2)*A + p1^2*s2^2 + p2*s1^4 - 4*s1^2*p1*p2 # = 0
p1 + p2 - R4 # = 0

### Solution Tr. System:
(4*s1^2 + s2^2)*p1^2 + (2*s2*A - s1^4 - 4*R4*s1^2)*p1 + R4*s1^4 + A^2 - s1^2*s2*A # = 0


### Solver:

solve.S4C2.E21aVar = function(R, debug=TRUE, all=FALSE) {
	s1 = R[1]; s2 = R[2]; E21a = R[3]; R4 = R[4];
	# Special Cases:
	isSpecial = FALSE;
	if(round0(s1^2*s2 - 4*E21a) == 0) {
		isSpecial = TRUE;
		warning("Special case!");
		x = s1/2;
		y = roots(c(1, -s2, R4 - x^2));
		sol0 = cbind(x1=x, x2=x, y1=y[1], y2=y[2]);
		# 2nd Set:
		if(round0(s1^2 + s2^2 - 4*R4) == 0) {
			# another set: x1 == x3;
			# coeff = c();
			# Note: redundant; y = s2/2;
			if(all) sol0 = rbind(sol0, sol0[ , c(2,1,4,3)]);
			return(sol0);
		} else {
			coeff = c(16*s1^2 + 4*s2^2, s1^2*(3*s2^2 - 16*R4));
		}
	} else {
		coeff = c((4*s1^2 + s2^2), (2*s2*E21a - s1^4 - 4*R4*s1^2),
			R4*s1^4 + E21a^2 - s1^2*s2*E21a);
	}
	#
	p1 = roots(coeff);
	if(debug) print(p1);
	# Step 2:
	x1 = sapply(p1, function(p) roots(c(1, -s1, p)));
	x1 = as.vector(x1);
	x2 = s1 - x1;
	# Step 3: robust
	# p1, p2 = NOT needed;
	y1 = (x2^2*s2 - E21a) / (x2^2 - x1^2);
	y2 = s2 - y1;
	#
	sol = cbind(x1, x2, y1, y2);
	if(isSpecial) sol = rbind(sol, sol0);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
test.S4C2.E21aVar = function(sol, R=NULL) {
	test.S4C2.E2aP1(sol=sol, R=R, n=c(2,1), type=1);
}

### Examples:

### Ex 1:
R = c(-2,-3,5,-1)
sol = solve.S4C2.E21aVar(R)

test.S4C2.E21aVar(sol)


### Ex 2: Special Case
s1 = 2; s2 = -3;
R = c(s1, s2, s2*s1^2/4, -1)
sol = solve.S4C2.E21aVar(R)

print(R)
test.S4C2.E21aVar(sol)


### Ex 3: Special Case
s1 = 2; s2 = -3;
R = c(s1, s2, s1^2*s2/4, (s1^2 + s2^2) / 4);
sol = solve.S4C2.E21aVar(R)

print(R)
test.S4C2.E21aVar(sol)


############################
############################

### Case: x^3*y
# n3 = 3; n4 = 1;

### System:
# x1 + x2 = R1
# y1 + y2 = R2
# x1^3*y1 + x2^3*y2 = R3
# x1*x2*y1*y2 = R4


### Solution:

### Transform:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;

A = x1^3*y1 + x2^3*y2;
B = x1^3*y2 + x2^3*y1;

### A + B =>
A + B - (y1 + y2)*(x1^3 + x2^3) # = 0
A + B - s2*(s1^3 - 3*p1*s1) # = 0
A + B + 3*s1*s2*p1 - s1^3*s2 # = 0

### A * B =>
A * B - p1^3*(y1^2 + y2^2) - p2*(x1^6 + x2^6) # = 0
A * B - p1^3*(s2^2 - 2*p2) - p2*(s1^6 - 6*p1*s1^4 + 9*p1^2*s1^2 - 2*p1^3) # = 0
A * B - p1^3*s2^2 - p2*(s1^6 - 6*p1*s1^4 + 9*p1^2*s1^2) + 4*p1^3*p2 # = 0
A * B - p1^3*s2^2 - p2*s1^6 + 6*s1^4*E4 - 9*p1*s1^2*E4 + 4*p1^2*E4 # = 0

# =>
A^2 + (3*s1*s2*p1 - s1^3*s2)*A + A*B # = 0
# where (A*B) is the previous eq;


### Transformed System:
s1 - R1 # = 0
s2 - R2 # = 0
A - R3 # = 0
s2^2*p1^3 - 4*E4*p1^2 + 3*s1*(3*s1*E4 + s2*A)*p1 + s1^6*p2 - 6*s1^4*E4 + A^2 - s1^3*s2*A # = 0
p1*p2 - R4 # = 0

### Solution Tr. System:
s2^2*p1^4 - 4*E4*p1^3 + 3*s1*(3*s1*E4 + s2*A)*p1^2 +
	- (6*s1^4*E4 - A^2 + s1^3*s2*A)*p1 + s1^6*E4 # = 0


### Solver:

solve.S4C2.E31a = function(R, debug=TRUE, all=FALSE) {
	s1 = R[1]; s2 = R[2]; E31a = R[3]; E4 = R[4];
	# Special Cases:
	isSpecial = FALSE;
	if(round0(s1^3*s2 - 8*E31a) == 0) {
		isSpecial = TRUE;
		warning("Special case!");
		x = s1/2;
		y = roots(c(1, -s2, E4 / x^2));
		sol0 = cbind(x1=x, x2=x, y1=y[1], y2=y[2]);
		if(all) sol0 = rbind(sol0, sol0[, c(2,1,4,3)]);
		# 2nd Set:
		if(round0((s1*s2)^2 - 16*E4) == 0) {
			# another set: x1 == x3;
			coeff = c(4, s1^2, 4*s1^4);
		} else {
			coeff = c(16*s2^2, 4*s1^2*s2^2 - 64*E4, 7*s1^4*s2^2 + 128*s1^2*E4, - 64*s1^4*E4);
		}
		# Note: includes also the permutation;
		# - but useful for generating the Classic polynomial;
	} else {
		coeff = c(s2^2, - 4*E4, 3*s1*(3*s1*E4 + s2*E31a),
			- (6*s1^4*E4 - E31a^2 + s1^3*s2*E31a), s1^6*E4);
	}
	#
	p1 = roots(coeff);
	if(debug) print(p1);
	# Step 2:
	x1 = sapply(p1, function(p) roots(c(1, -s1, p)));
	x1 = as.vector(x1);
	x2 = s1 - x1;
	# Step 3: robust
	# p1, p2 = NOT needed;
	y1 = (x2^3*s2 - E31a) / (x2^3 - x1^3);
	y2 = s2 - y1;
	#
	sol = cbind(x1, x2, y1, y2);
	if(isSpecial) sol = rbind(sol, sol0);
	if(all && ! isSpecial) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
test.S4C2.E31a = function(sol, R=NULL) {
	test.S4C2.E2aP1(sol=sol, R=R, n=c(3,1));
}

### Examples:

### Ex 1:
R = c(-2,-3,5,-1)
sol = solve.S4C2.E31a(R)

test.S4C2.E31a(sol)


### Ex 2: Special Case
s1 = 2
R = c(s1, -3, -3*s1^3/8, -1)
sol = solve.S4C2.E31a(R)

print(R)
test.S4C2.E31a(sol)


### Ex 3: Special Case
s1 = 2; s2 = -3;
R = c(s1, s2, s1^3*s2/8, (s1*s2)^2 / 16);
sol = solve.S4C2.E31a(R)

print(R)
test.S4C2.E31a(sol)


### Ex 4:
R = c(1,-2,5, -6)
sol = solve.S4C2.E31a(R)

test.S4C2.E31a(sol)

poly.calc(sol[,1]) * 4


### Classic Poly:

x1 = sol[,1];
s1 = R[1]; s2 = R[2]; R3 = R[3]; E4 = R[4];
s2^2*x1^8 - 4*s1*s2^2*x1^7 + (6*s1^2*s2^2 + 4*E4)*x1^6 - 4*s1*(s1^2*s2^2 + 3*E4)*x1^5 +
	+ (s1^4*s2^2 + 21*s1^2*E4 + 3*R3*s1*s2)*x1^4 - s1^2*(22*s1*E4 + 6*R3*s2)*x1^3 +
	+ (15*s1^4*E4 + 4*R3*s1^3*s2 - R3^2)*x1^2 +
	- (6*s1^5*E4 + R3*s1^4*s2 - R3^2*s1)*x1 + s1^6*E4 # = 0

############
### Variant:

### System:
# x1 + x2 = R1
# y1 + y2 = R2
# x1^3*y1 + x2^3*y2 = R3
# x1*x2 + y1*y2 = R4

### Transformed System:
s1 - R1 # = 0
s2 - R2 # = 0
A - R3 # = 0
s2^2*p1^3 - 4*p1^3*p2 + 3*s1*(3*s1*p1*p2 + s2*A)*p1 + s1^6*p2 - 6*s1^4*p1*p2 + A^2 - s1^3*s2*A # = 0
p1 + p2 - R4 # = 0

### Solution Tr. System:
4*p1^4 - (9*s1^2 - s2^2 + 4*R4)*p1^3 + (6*s1^4 + 9*s1^2*R4)*p1^2 +
	- (s1^6 + 6*s1^4*R4 - 3*s1*s2*E31a)*p1 +
	+ s1^6*R4 - s2*s1^3*E31a + E31a^2 # = 0


### Solver:

solve.S4C2.E31aVar = function(R, debug=TRUE, all=FALSE) {
	s1 = R[1]; s2 = R[2]; E31a = R[3]; R4 = R[4];
	# Special Cases:
	isSpecial = FALSE;
	if(round0(s1^3*s2 - 8*E31a) == 0) {
		isSpecial = TRUE;
		warning("Special case!");
		x = s1/2;
		y = roots(c(1, -s2, R4 - x^2));
		sol0 = cbind(x1=x, x2=x, y1=y[1], y2=y[2]);
		if(all) sol0 = rbind(sol0, sol0[, c(2,1,4,3)]);
		# 2nd Set:
		if(round0(s1^2 + s2^2 - 4*R4) == 0) {
			# another set: x1 == x3;
			s1sq = s1^2;
			coeff = c(16, -32*s1sq, s1sq*(16*s1sq + 9*s2^2));
		} else {
			coeff = c(64, - 128*s1^2 + 16*s2^2 - 64*R4,
				4*s1^2*(16*s1^2 + s2^2 + 32*R4), s1^4*(7*s2^2 - 64*R4) );
		}
		# Note: includes also the permutation;
		# - but useful for generating the Classic polynomial;
	} else {
		coeff = c(4, - (9*s1^2 - s2^2 + 4*R4), (6*s1^4 + 9*s1^2*R4),
			- (s1^6 + 6*s1^4*R4 - 3*s1*s2*E31a),
			s1^6*R4 - s2*s1^3*E31a + E31a^2);
	}
	#
	p1 = roots(coeff);
	if(debug) print(p1);
	# Step 2:
	x1 = sapply(p1, function(p) roots(c(1, -s1, p)));
	x1 = as.vector(x1);
	x2 = s1 - x1;
	# Step 3: robust
	# p1, p2 = NOT needed;
	y1 = (x2^3*s2 - E31a) / (x2^3 - x1^3);
	y2 = s2 - y1;
	#
	sol = cbind(x1, x2, y1, y2);
	if(isSpecial) sol = rbind(sol, sol0);
	if(all && ! isSpecial) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
test.S4C2.E31aVar = function(sol, R=NULL) {
	test.S4C2.E2aP1(sol=sol, R=R, n=c(3,1), type=1);
}

### Examples:

### Ex 1:
R = c(-2,-3,5,-1)
sol = solve.S4C2.E31aVar(R)

test.S4C2.E31aVar(sol)


### Ex 2: Special Case
s1 = 2
R = c(s1, -3, -3*s1^3/8, -1)
sol = solve.S4C2.E31aVar(R)

print(R)
test.S4C2.E31aVar(sol)


### Ex 3: Special Case
s1 = 2; s2 = -3;
R = c(s1, s2, s1^3*s2/8, (s1^2 + s2^2) / 4);
sol = solve.S4C2.E31aVar(R)

print(R)
test.S4C2.E31aVar(sol)


### Ex 4:
R = c(1,-3,1,-2)
sol = solve.S4C2.E31aVar(R)

test.S4C2.E31aVar(sol)


############################
############################

##################
### 2 Eqs E21a ###
##################

### System:
# x1^n1 + x2^n1 = R1
# y1^n2 + y2^n2 = R2
# x1^2*y1 + x2^2*y2 = R3
# x1*y1^2 + x2*y2^2 = R4

### Symmetries:
# - if (x1, x2, y1, y2) is a solution,
#   so is also the C2 permutation: (x2, x1, y2, y1);

### Case:
# n1 = n2 = 1;

### Solution:

### Transform:
s1 = x1 + x2; s2 = y1 + y2;
p1 = x1 * x2; p2 = y1 * y2;

A1 = x1^2*y1 + x2^2*y2;
A2 = x1*y1^2 + x2*y2^2;

### C2-Transform:
### Transformed System:
s1 - R1 # = 0
s2 - R2 # = 0
# A1 = R3; A2 = R4;
A1^2 + (2*s2*p1 - s1^2*s2)*A1 + p1^2*s2^2 + p2*s1^4 - 4*s1^2*p1*p2 # = 0
A2^2 + (2*s1*p2 - s2^2*s1)*A2 + p2^2*s1^2 + p1*s2^4 - 4*s2^2*p1*p2 # = 0

### Solution Tr. System:
3*s2^2*p1^2 - (s1^2*s2^2 - 2*s2*A1 + 4*s1*A2)*p1 - A1^2 + s1^3*A2 # = 0

### Debug:
x = sqrt(c(2,3,5,7));
x[3] = - x[3];
x1 = x[1]; x2 = x[3]; y1 = x[2]; y2 = x[4];
R1 = x1 + x2; R2 = y1 + y2;
R3 = x1^2*y1 + x2^2*y2;
R4 = x1*y1^2 + x2*y2^2;

### Solver:

solve.S4C2.E21x2 = function(R, debug=TRUE, all=TRUE) {
	coeff = coeff.S4C2.E21x2(R);
	p1 = roots(coeff);
	if(debug) print(p1);
	p2 = solve.S4C2.E21x2.p2(p1, R);
	# Step 2:
	s1 = R[1]; s2 = R[2]; A1 = R[3]; A2 = R[4];
	x13 = sapply(p1, function(p) roots(c(1, -s1, p)));
	x13 = t(x13);
	# Step 3: robust
	y1 = (x13[,2]^2*s2 - A1) / (x13[,2]^2 - x13[,1]^2);
	y2 = s2 - y1;
	#
	sol = cbind(x13, y1, y2);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
coeff.S4C2.E21x2 = function(R) {
	s1 = R[1]; s2 = R[2]; A1 = R[3]; A2 = R[4];
	coeff = c(3*s2^2, - (s1^2*s2^2 - 2*s2*A1 + 4*s1*A2), - A1^2 + s1^3*A2);
	return(coeff);
}
solve.S4C2.E21x2.p2 = function(p1, R) {
	s1 = R[1]; s2 = R[2]; A1 = R[3]; A2 = R[4];
	s1sq = s1^2;
	p2 = (p1*s2 + A1)^2 - s1sq*s2*A1;
	p2 = p2 / (s1sq*(4*p1 - s1sq));
	return(p2);
}
test.S4C2.E21x2 = function(sol, R=NULL) {
	x1 = sol[,1]; x2 = sol[,2];
	y1 = sol[,3]; y2 = sol[,4];
	err1 = x1 + x2; err2 = y1 + y2;
	err3 = x1^2*y1 + x2^2*y2;
	err4 = x1*y1^2 + x2*y2^2;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		err = err - R;
	}
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = c(-2,-1,1,3)
sol = solve.S4C2.E21x2(R)

test.S4C2.E21x2(sol)

poly.calc(sol[,1]) * 3


### Ex 2:
R = c(1,-2,1,3)
sol = solve.S4C2.E21x2(R)

test.S4C2.E21x2(sol)


### Ex 3:
R = c(-1,-2,5,3)
sol = solve.S4C2.E21x2(R)

test.S4C2.E21x2(sol)


### Derivation:
pA1 = toPoly.pm("A1^2 + (2*s2*p1 - s1^2*s2)*A1 + p1^2*s2^2 + p2*s1^4 - 4*s1^2*p1*p2");
pA2 = toPoly.pm("A2^2 + (2*s1*p2 - s2^2*s1)*A2 + p2^2*s1^2 + p1*s2^4 - 4*s2^2*p1*p2");
pR = solve.pm(pA1, pA2, "p2")
pR$Rez = sort.pm(pR$Rez, c("s1", "s2"))
pR$Rez$coeff = - pR$Rez$coeff;
str(pR)

pR = pR$Rez;
pDiv = toPoly.pm("3*s2^2*p1^2 - (s1^2*s2^2 - 2*s2*A1 + 4*s1*A2)*p1 - A1^2 + s1^3*A2")
div.pm(pR, pDiv, "s1")
toCoeff(pDiv, "p1")

### Factorization (workout)
# Coefficient b1:
# pR = full P[4]; # before div;
tmp = pR[pR$p1 == 1, ]
tmp$coeff = tmp$coeff * 3;
tmp$p1 = 0; tmp = drop.pm(tmp);
pB01 = toPoly.pm("- A1^2 + s1^3*A2")
pB02 = toPoly.pm("s1^4*s2^2 - 2*s1^2*s2*A1 - s1^3*A2 + A1^2")
pB3 = toPoly.pm("- 20*s1^2*s2^2 + 28*s2*A1 - 8*s1*A2")
tmp = diff.pm(tmp, mult.pm(pB3, pB01))
div.pm(tmp, diff.pm(mult.pm(pB02, 3), mult.pm(pB01, 5)), "s1")

