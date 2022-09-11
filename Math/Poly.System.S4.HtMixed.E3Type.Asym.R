########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### E3-Type: Asymmetric
###
### draft v.0.1e


### E3-Type:
# SUM(x1^n1 * x2^n2 * x3^n3) = R3;
# where: n1 != n3 (Asymmetric);


####################
####################

### Helper Functions

source("Poly.System.S4.HtMixed.Basic.Helper.R")

test.S4HtMixed.E211a = function(sol, R = NULL, n=1) {
	err = test.S4HtMixed.En3(sol, R=R, n=1, nE=c(2,1,1), type="E2");
	# Correct E2 order;
	# TODO: move correction to test.S4HtMixed.En3 ???
	err = err[ c(1,3,2,4), ];
	return(err);
}
test.S4HtMixed.E422a = function(sol, R = NULL, n=1) {
	err = test.S4HtMixed.En3(sol, R=R, n=1, nE=c(4,2,2), type="E2");
	# Correct E2 order;
	# TODO: move correction to test.S4HtMixed.En3 ???
	err = err[ c(1,3,2,4), ];
	return(err);
}


####################

####################
### Type: E211a  ###
####################

E211a = x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2;

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - R2 # = 0
x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0

### Solution:

### Eq 3:
# - see formula for E211a in file:
#   Poly.System.S4.C2.Formulas.R;


### Transformed System:
S - R1 # = 0
sp + ps - R2 # = 0
E211a - R3 # = 0
E4 - R4 # = 0

### Solver:

solve.S4Ht.E211a = function(R, debug=TRUE, all=FALSE) {
	coeff = coeff.S4Ht.E211a(R);
	ps = roots(coeff);
	sp = R[2] - ps;
	if(debug) print(ps);
	# Step 2:
	s1 = sapply(ps, function(ps) {
		roots(c(1, -R[1], ps));
	})
	s1 = as.vector(s1);
	s2 = R[1] - s1;
	sp = rep(sp, each = 2);
	ps = rep(ps, each = 2);
	# Step 3: robust;
	p1 = solve.S4Ht.E211a.p1(list(ps=ps, sp=sp, s1=s1, s2=s2), R=R);
	p2 = sp - p1;
	# Step 4:
	len = length(s1);
	x13 = sapply(seq(len), function(id) {
		roots(c(1, -s1[id], p1[id]));
	})
	x13 = t(x13);
	x1 = x13[,1]; x3 = x13[,2];
	# x2*(p1*x1 + p2*x3) + x4*(p1*x3 + p2*x1) = R3;
	x2  = (p1*x3 + p2*x1)*s2 - R[3];
	div = (x1 - x3)*(p1 - p2);
	x2  = - x2 / div;
	x4  = s2 - x2;
	#
	sol = cbind(x1, x2, x3, x4);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
coeff.S4Ht.E211a = function(R) {
	S = R[1]; R2 = R[2]; E211a = R[3]; E4 = R[4];
	coeff = c(1, - 6*R2, 2*E211a + 15*R2^2 - 14*E4,
		- S^2*E211a - E4*S^2 - 8*E211a*R2 - 20*R2^3 + 66*E4*R2,
		E4*S^4 + 4*S^2*E211a*R2 + 3*E211a^2 - E4*S^2*R2 + 12*E211a*R2^2 - 14*E4*E211a +
			+ 15*R2^4 - 122*E4*R2^2 + 81*E4^2,
		- 4*E4*S^4*R2 - S^2*E211a^2 - 6*S^2*E211a*R2^2 + 4*E4*S^2*E211a - 8*E211a^2*R2 +
			+ 13*E4*S^2*R2^2 + 4*E4^2*S^2 - 8*E211a*R2^3 + 38*E4*E211a*R2 - 6*R2^5 +
			+ 110*E4*R2^3 - 264*E4^2*R2,
		6*E4*S^4*R2^2 - 8*E4^2*S^4 + 3*S^2*E211a^2*R2 + 2*E211a^3 + 4*S^2*E211a*R2^3 +
			- 8*E4*S^2*E211a*R2 + 7*E211a^2*R2^2 - 14*E4*E211a^2 - 23*E4*S^2*R2^3 +
			+ 28*E4^2*S^2*R2 + 2*E211a*R2^4 - 32*E4*E211a*R2^2 + 32*E4^2*E211a + R2^6 +
			- 48*E4*R2^4 + 296*E4^2*R2^2 - 224*E4^3,
		- 4*E4*S^4*R2^3 + 16*E4^2*S^4*R2 - 3*S^2*E211a^2*R2^2 + 4*E4*S^2*E211a^2 +
			- 2*E211a^3*R2 - S^2*E211a*R2^4 + 4*E4*S^2*E211a*R2^2 - 2*E211a^2*R2^3 +
			+ 24*E4*E211a^2*R2 + 16*E4*S^2*R2^4 - 64*E4^2*S^2*R2^2 + 8*E4*E211a*R2^3 +
			- 32*E4^2*E211a*R2 + 8*E4*R2^5 - 128*E4^2*R2^3 + 384*E4^3*R2,
		E211a^4 + E4*S^4*R2^4 - 8*E4^2*S^4*R2^2 + 16*E4^3*S^4 + S^2*E211a^2*R2^3 +
			- 4*E4*S^2*E211a^2*R2 - 8*E4*E211a^2*R2^2 + 32*E4^2*E211a^2 - 4*E4*S^2*R2^5 +
			+ 32*E4^2*S^2*R2^3 - 64*E4^3*S^2*R2 + 16*E4^2*R2^4 - 128*E4^3*R2^2 + 256*E4^4);
	return(coeff);
}
solve.S4Ht.E211a.p1 = function(s, R) {
	ps = s$ps; sp = s$sp; s1 = s$s1; s2 = s$s2;
	E211a = R[3]; E4 = R[4];
	p1 = - ps^2*E4^2*sp + E211a*ps*E4*sp^2 - E4*s2^2*sp^4 - E4^2*s1^2*sp^2 + 5*E4^2*s2^2*sp^2 +
		+ 4*E4^3*s1^2 - 4*E4^3*s2^2 + 4*E4^2*sp^3 - 16*E4^3*sp - E211a^2*E4*sp;
	div = ps^2*E4*sp^2 - ps^2*E4^2 - E211a*ps*sp^3 + E211a*ps*E4*sp + s2^2*sp^5 +
		+ E4*s1^2*sp^3 - 6*E4*s2^2*sp^3 - 4*E4*sp^4 - 4*E4^2*s1^2*sp + 8*E4^2*s2^2*sp +
		+ E211a^2*sp^2 + 20*E4^2*sp^2 - 16*E4^3 - E211a^2*E4;
	p1 = - p1 / div;
	return(p1);
}

### Examples:

### Ex 1:
R = c(-1,2,3,-2)
sol = solve.S4Ht.E211a(R);

test.S4HtMixed.E211a(sol)


### Ex 2:
R = c(-5,2,3,-2)
sol = solve.S4Ht.E211a(R);

test.S4HtMixed.E211a(sol)


### Ex 3:
# - small precision error in one of the roots;
R = c(-5,2,-3,1)
sol = solve.S4Ht.E211a(R);

test.S4HtMixed.E211a(sol)


###############
### Derivation:
p1 = toPoly.pm(...) # copy polynomial for E211a;
pSp = toPoly.pm("R2 - ps")

pR = replace.pm(p1, pSp, "sp")
pR = sort.pm(pR, c("ps", "E211a", "S"))

toCoeff(pR, "ps", print=TRUE)


#####################
#####################

###################
### Type: E422a ###
###################

E422a = x1^4*x2^2*x3^2 + x2^4*x3^2*x4^2 + x3^4*x4^2*x1^2 + x4^4*x1^2*x2^2;


coeffFactory = function() {
	pE = toPoly.pm("E422a^2 - SA2B2*E422a + A2B2")
	pSA2B2 = toPoly.pm("- 2*(sp^2*S - 2*E4*S)*p1s1 +
		- (sp^2 - 2*E4)*(- 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4)"); # actually: - (A2 + B2);
	pSA2B2$coeff = - pSA2B2$coeff;
	pA2B2 = toPoly.pm("- E4^2*((2*S*p1s1 - 2*sp*S^2 + 2*sp*ps + ps^2 + 4*E4)^2 - 2*A1B1) +
		- A1B1*(sp^4 - 4*E4*sp^2 + 2*E4^2)"); # actually: - (A2 * B2);
	pA2B2$coeff = - pA2B2$coeff;
	pA1B1 = toPoly.pm("(4*E4*S + 2*sp*ps*S - sp*S^3)*p1s1 +
		- E4*(S^4 + 4*sp*S^2 - 4*ps*S^2 + 2*ps^2 - 4*sp*ps) + 4*E4^2 +
		+ sp^2*S^4 + sp^2*ps^2 - 3*ps*sp^2*S^2");
	pR = replace.pm(pA2B2, pA1B1, "A1B1")
	#
	pE = replace.pm(pE, pSA2B2, "SA2B2");
	pE = replace.pm(pE, pR, "A2B2");
	#
	pP1S1 = toPoly.pm("p1s1^2 - sp*S*p1s1 + ps*sp^2 + E4*S^2 - 4*ps*E4")
	pR = solve.pm(pE, pP1S1, "p1s1")
	# 213 Monomials;
	pR = pR$Rez
	#
	pE2 = toPoly.pm("ps + sp - E2");
	pR = solve.pm(pE2, pR, "sp");
	pR = pR$Rez
	pR = sort.pm(pR, "ps", xn2="E422a");
	#
	coeff = as.coeff.pm(pR, "ps");
	return(invisible(coeff));
}
eval.lpm = function(p, vals) {
	len = length(p);
	rez = sapply(seq(len), function(id) {
		p1 = p[[id]];
		if(is.numeric(p1)) return(p1);
		rez = eval.pm(p1, vals);
		return(rez);
	});
	return(rez);
}
solverFactory = function() {
	coeff.p = coeffFactory();
	FUN = function(R, debug=TRUE, all=FALSE) {
		coeff = eval.lpm(coeff.p, list(S=R[1], E2=R[2], E422a=R[3], E4=R[4]));
		ps = roots(coeff);
		sp = R[2] - ps;
		if(debug) { print(coeff); print(sp); }
		# Step 2:
		S = R[1];
		s1 = sapply(ps, function(ps) {
			roots(c(1, - S, ps));
		})
		s1 = as.vector(s1);
		s2 = S - s1;
		sp = rep(sp, each=2);
		ps = rep(ps, each=2);
		# Step 3:
		# TODO: robust
		p1 = sapply(sp, function(sp) {
			roots(c(1, - sp, R[4]));
		})
		p1 = as.vector(p1);
		sp = rep(sp, each=2);
		p2 = sp - p1;
		s1 = rep(s1, each=2); s2 = rep(s2, each=2);
		# Step 4:
		len = length(s1);
		x13 = sapply(seq(len), function(id) {
			roots(c(1, - s1[id], p1[id]));
		});
		x13 = t(x13);
		x1 = x13[,1]; x3 = x13[,2];
		#
		R24 = R[3] + p2*(x1^2 + x3^2)*(p1^2 + p2^2);
		c1  = s2*((x1*p1)^2 + (x3*p2)^2);
		c2  = s2*((x3*p1)^2 + (x1*p2)^2);
		div = c2 - c1;
		x2 = (c2*s2 - R24) / div;
		x4 = s2 - x2;
		#
		sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
		if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
		return(sol);
	}
	return(FUN);
}
solver.S4Ht.E422a = solverFactory();

### Examples:

# TODO:
# robust: only 1/2 are TRUE roots;

###
R = c(-1, 2, -3, 3)
sol = solver.S4Ht.E422a(R)

test.S4HtMixed.E422a(sol)


#####################
#####################

#####################
### Type: 2 Eqs   ###
### E211a + E112a ###
#####################

E211a = x1^2*x2*x3 + x2^2*x3*x4 + x3^2*x4*x1 + x4^2*x1*x2;
E112a = x1*x2*x3^2 + x2*x3*x4^2 + x3*x4*x1^2 + x4*x1*x2^2;

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2 + x1*x3 + x1*x4 + x2*x3 + x2*x4 + x3*x4 - R2 # = 0
E211a - R3 # = 0
E112a - R4 # = 0


### Solution:

### Eqs 3 & 4:
# - see formula for E211a &
#   for the sum (E211a + E112a) in file:
#   Poly.System.S4.C2.Formulas.R;


### Transformed System:
S - R1 # = 0
sp + ps - R2 # = 0
E211a - R3 # = 0
sp*ps - R3 - R4 # = 0

### Solver:

# TODO: robust
solve.S4Ht.E211aE112a = function(R, debug=TRUE, all=FALSE) {
	sp = roots(c(1, -R[2], R[3] + R[4]));
	ps = R[2] - sp;
	### Step 2:
	# TODO: robust
	len = length(sp);
	E4 = sapply(seq(len), function(id) {
		coeff = coeff.S4Ht.E211aE112a(list(sp=sp[id], ps=ps[id]), R=R);
		roots(coeff);
	})
	E4 = as.vector(E4);
	sp = rep(sp, each=2);
	ps = rep(ps, each=2);
	if(debug) print(E4);
	### Step 3:
	# TODO: robust
	len = length(ps);
	s1 = sapply(seq(len), function(id) {
		roots(c(1, - R[1], ps[id]));
	})
	s1 = as.vector(s1);
	s2 = R[1] - s1;
	sp = rep(sp, each=2);
	ps = rep(ps, each=2);
	E4 = rep(E4, each=2);
	### Step 4:
	p1 = solve.S4Ht.E211aE112a.p1(list(s1=s1, s2=s2, sp=sp, ps=ps, E4=E4), R);
	p2 = sp - p1;
	### Step 5:
	len = length(s1);
	x13 = sapply(seq(len), function(id) {
		roots(c(1, - s1[id], p1[id]));
	})
	x13 = t(x13);
	x1 = x13[,1]; x3 = x13[,2];
	# x2*(p1*x1 + p2*x3) + x4*(p1*x3 + p2*x1) = R3;
	x2  = (p1*x3 + p2*x1)*s2 - R[3];
	div = (x1 - x3)*(p1 - p2);
	x2  = - x2 / div;
	x4  = s2 - x2;
	#
	sol = cbind(x1, x2, x3, x4);
	if(all) sol = rbind(sol, sol[ , c(2,1,4,3)]);
	return(sol);
}
coeff.S4Ht.E211aE112a = function(s, R) {
	S = R[1]; E211a = R[3];
	sp = s$sp; ps = s$ps;
	coeff = c(256,
		- 128*sp^2 + 128*sp*ps + 32*ps^2 - 64*sp*S^2 - 64*ps*S^2 + 16*S^4,
		16*sp^4 - 64*sp^3*ps + 8*sp^2*ps^2 + 8*sp*ps^3 + ps^4 + 32*sp^3*S^2 +
			+ 32*sp^2*ps*S^2 - 4*sp*ps^2*S^2 - 8*sp^2*S^4 - 32*sp*ps*E211a + 32*E211a^2,
		8*sp^5*ps - 8*sp^4*ps^2 - 2*sp^3*ps^3 - 4*sp^5*S^2 - 4*sp^4*ps*S^2 + sp^3*ps^2*S^2 +
			+ sp^4*S^4 + 8*sp^3*ps*E211a - 8*sp^2*ps^2*E211a - 2*sp*ps^3*E211a + 4*sp^2*ps*S^2*E211a +
			- 8*sp^2*E211a^2 + 8*sp*ps*E211a^2 + 2*ps^2*E211a^2 - 4*sp*S^2*E211a^2,
		sp^6*ps^2 + 2*sp^4*ps^2*E211a - sp^4*ps*S^2*E211a - 2*sp^3*ps*E211a^2 + sp^2*ps^2*E211a^2 +
			+ sp^3*S^2*E211a^2 - 2*sp*ps*E211a^3 + E211a^4);
	return(coeff);
}
solve.S4Ht.E211aE112a.p1 = function(s, R) {
	ps = s$ps; sp = s$sp; s1 = s$s1; s2 = s$s2;
	E4 = s$E4;;
	E211a = R[3];
	# same as for E211a System;
	p1 = - ps^2*E4^2*sp + E211a*ps*E4*sp^2 - E4*s2^2*sp^4 - E4^2*s1^2*sp^2 + 5*E4^2*s2^2*sp^2 +
		+ 4*E4^3*s1^2 - 4*E4^3*s2^2 + 4*E4^2*sp^3 - 16*E4^3*sp - E211a^2*E4*sp;
	div = ps^2*E4*sp^2 - ps^2*E4^2 - E211a*ps*sp^3 + E211a*ps*E4*sp + s2^2*sp^5 +
		+ E4*s1^2*sp^3 - 6*E4*s2^2*sp^3 - 4*E4*sp^4 - 4*E4^2*s1^2*sp + 8*E4^2*s2^2*sp +
		+ E211a^2*sp^2 + 20*E4^2*sp^2 - 16*E4^3 - E211a^2*E4;
	p1 = - p1 / div;
	return(p1);
}
test.S4HtMixed.E211aE112a = function(sol, R=NULL) {
	err = test.S4HtMixed.En3(sol, R=R, n=1, nE=c(2,1,1), type="E2");
	err = err[ c(1,3,2,4), ];
	#
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	E112a = x1*x2*x3^2 + x2*x3*x4^2 + x3*x4*x1^2 + x4*x1*x2^2;
	err[4, ] = E112a;
	return(err);
}

### Examples

### Ex 1:
R = c(-1,2,-3,-2)
sol = solve.S4Ht.E211aE112a(R)

test.S4HtMixed.E211aE112a(sol)


### Ex 2:
R = c(2,-3,-1,4)
sol = solve.S4Ht.E211aE112a(R)

test.S4HtMixed.E211aE112a(sol)

