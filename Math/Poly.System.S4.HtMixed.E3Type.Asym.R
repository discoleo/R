########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### E3-Type: Asymmetric
###
### draft v.0.1a


### E3-Type:
# SUM(x1^n1 * x2^n2 * x3^n3) = R3;
# where: n1 != n3 (Asymmetric);


####################
####################

### Helper Functions

source("Poly.System.S4.HtMixed.Basic.Helper.R")


####################

####################
### Type: E211a  ###
####################

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
# Eq E211a - R3 = 0
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
	# Step 3:
	# TODO: robust;
	len = length(sp);
	p1 = sapply(seq(len), function(id) {
		roots(c(1, -sp[id], R[4]));
	})
	s1 = rep(s1, each=2); s2 = rep(s2, each=2);
	sp = rep(sp, each=2);
	p1 = as.vector(p1);
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

### Examples:

R = c(-1,2,3,-2)
sol = solve.S4Ht.E211a(R);

# TODO: robust;
test.S4HtMixed.En3(sol, n=1, nE=c(2,1,1), type="E2")


### Derivation:
p1 = toPoly.pm(...) # copy polynomial for E211a;
pSp = toPoly.pm("R2 - ps")

pR = replace.pm(p1, pSp, "sp")
pR = sort.pm(pR, c("ps", "E211a", "S"))

toCoeff(pR, "ps", print=TRUE)

