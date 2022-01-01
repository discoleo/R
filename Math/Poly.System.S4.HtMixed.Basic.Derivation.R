########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Hetero-Symmetric S4: Mixed
### Basic Types: Derivation
###
### draft v.0.1b


### Derivation of the Formulas


########################

### Helper Functions

source("Poly.System.S4.HtMixed.Basic.Helper.R")


########################
########################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2 + x2*x3 + x3*x4 + x4*x1 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0


### Derivation:

# - classic approach: P[2] o P[2];

### E2a:
(x1+x3)*(S - x1 - x3) - R2 # = 0
# x13s = x1 + x3;
x13s^2 - S*x13s + R2 # = 0

### E3 =>
x1*x3*(S - x1 - x3) + x2*x4*(x1+x3) - R3 # = 0
(x1*x3)^2*(S - x1 - x3) - R3*x1*x3 + R4*(x1+x3) # = 0


### Solution: based on "classic" approach
solve.S4Ht.P1old = function(R, debug=FALSE) {
	xs  = roots(c(1, -R[1], R[2]));
	x13 = sapply(seq(length(xs)), function(id) roots(c(R[1] - xs[id], -R[3], R[4]*xs[id])));
	xs = rep(xs, each=2); x13 = as.vector(x13);
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	# x2, x4:
	xs = R[1] - xs; x24 = R[4] / x13;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1, x2, x3, x4)
	return(sol)
}

### Eq for E2:
R = c(1,-1,2,3)
sol = solve.S4Ht.P1old(R)
round0(poly.calc(apply(sol, 1, e2.f)[1:2]) * R[2])


### Eq:
R[2]*E2^2 - (R[1]*R[3] + 2*R[2]^2)*E2 +
	+ R[1]^2*R[4] + R[1]*R[2]*R[3] + R[2]^3 - 4*R[2]*R[4] + R[3]^2

test.S4HtMixed(sol, n=1)


####################
### E2a: Order 2 ###
####################

###############
### Order 1 ###
###############

x1 + x2 + x3 + x4 - R1 # = 0
(x1*x2)^2 + (x2*x3)^2 + (x3*x4)^2 + (x4*x1)^2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0


###############
### Derivation:

# - classic approach: P[2] o P[2];

### E2a:
# let: xs = x1 + x3 => x2 + x4 = S - xs;
(x1^2 + x3^2)*((S - xs)^2 - 2*x2*x4) - R2 # = 0
(xs^2 - 2*x1*x3)*(x1*x3*(S - xs)^2 - 2*R4) - R2*x1*x3 # = 0


### E3 =>
x1*x3*(S - x1 - x3) + x2*x4*(x1+x3) - R3 # = 0
(x1*x3)^2*(S - xs) - R3*x1*x3 + R4*xs # = 0


### Solution: based on "classic" approach
coeff.S4HtM.Ord2.P1 = function(R) {
	S = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c(1, - 4*S, 6*S^2, - 4*S^3, S^4 - 2*S*R3 - 8*R4 - 2*R2,
		4*S^2*R3 + 16*R4*S + 4*R2*S,
		- 2*S^3*R3 - 12*R4*S^2 - 2*R2*S^2 - 4*R3^2,
		4*R4*S^3 + 4*S*R3^2,
		- 8*R4*S*R3 + 2*R2*S*R3 + 16*R4^2 - 8*R2*R4 + R2^2);
	return(coeff);
}
solve.S4HtM.Ord2.P1old = function(R, debug=TRUE) {
	coeff = coeff.S4HtM.Ord2.P1(R);
	xs  = roots(coeff);
	if(debug) print(xs);
	#
	S = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	x13 = 4*R4*xs^3 - 6*R4*xs^2*S + 2*R4*xs*S^2;
	div = - R2*xs + 4*R4*xs + xs^5 + R2*S - 4*R4*S - 3*xs^4*S + 3*xs^3*S^2 - xs^2*S^3 +
		+ 2*xs^2*R3 - 4*xs*S*R3 + 2*S^2*R3;
	x13 = x13 / div;
	#
	xd = sqrt(xs^2 - 4*x13 + 0i);
	x1 = (xs + xd)/2; x3 = (xs - xd)/2;
	# x2, x4:
	xs = R[1] - xs; x24 = R[4] / x13;
	xd = sqrt(xs^2 - 4*x24 + 0i);
	x2 = (xs + xd)/2; x4 = (xs - xd)/2;
	sol = cbind(x1, x2, x3, x4)
	return(sol)
}

###
R = c(1,-1,2,3)
sol = solve.S4HtM.Ord2.P1old(R)

test.S4HtMixed(sol, n=1, nE2=2)


round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2])


R = c(1,  1,  1,  2,   7)
R = c(1,  1,  1,  3,  11)
R = c(1,  1,  1,  4,  15)
R = c(1,  2,  1,  2,   6)
R = c(1,  2,  1,  3,  10)
R = c(1,  2,  1,  4,  14)
R = c(1, -7,  1,  9,  43)
R = c(2, -3,  3, -5,  17)
R = c(2, -3,  3, sqrt(2)) # - 24
R = c(2, -3,  3, 2^(1/3)) # - 24
R = c(2, -3,  3, sqrt(3)) # - 24
R = c(2, -3,  3, 3^(1/3)) # - 24
R = c(2, -3,  3, sqrt(5)) # - 24
#
R = c(1, sqrt(2), 1, sqrt(2)) # -4
R = c(1, sqrt(2), 1, -sqrt(2)) # 12
R = c(1, -sqrt(2), 1, sqrt(2)) # -12
R = c(1, -sqrt(2), 3, sqrt(2)) # -36
R = c(1, -sqrt(3), 3, sqrt(3)) # -36
R = c(1, -sqrt(3), 5, sqrt(3)) # -60
#
R = c(1, 2*sqrt(3), 5, 1) # 24
which.coeff(R, sq=3)
#
which.coeff(R <- c(1, 2*sqrt(3), 5, 1), sq=3)
which.coeff(R <- c(1, 2*2^(1/3), 5, 1), sq=2^2, pow=3, DIFF=DIFF)

sol = solve.S4HtM.Ord2.P1old(R)
round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2])


### x = E2:
(4*R[4] - R[2])*x^4 - 2*(R[3]^2 + R[1]^2*R[4])*x^3 +
	+ (4*(R[1]*R[3] - 2*R[4])*R[2] - 8*R[1]*R[3]*R[4] + R[1]^2*R[3]^2 + 2*R[2]^2)*x^2 +
	+ (4*R[1]^3*R[3]*R[4] + 4*R[1]*R[3]^3 + 2*(R[1]^2*R[4] + R[3]^2)*R[2])*x +
	- R[1]^4*R[4]^2 - 2*R[1]^3*R[3]^3 - R[2]^3 - R[3]^4 + 4*R[2]^2*R[4] + 2*R[1]^2*R[3]^2*R[4] +
		+ 8*R[1]*R[2]*R[3]*R[4] - 4*R[1]*R[2]^2*R[3] - 5*R[1]^2*R[2]*R[3]^2;

whichHasPower(4, id=2, type=2)
whichHasPower(R <- c(1,1,-2,-3), id=2, type=2)
polyR(R)
which.sq(DIFF(R), sq=2)


### Hack Coefficients
polyR = function(R) {
	sol = solve.S4HtM.Ord2.P1old(R, debug=FALSE)
	p = round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2]);
	return(p);
}
which.coeff = function(R, sq=2, id=3, pow=2, DIFF=NULL, print=TRUE, digits=6, iter=1000) {
	sol = solve.S4HtM.Ord2.P1old(R, debug=FALSE)
	p = round0(poly.calc(e2.f(sol)[c(1,3,5,7)])) * (4*R[4] - R[2]);
	if(print) print(p);
	x = p[id];
	if( ! is.null(DIFF)) {
		x = x - DIFF(R);
	}
	return(which.sq(x, sq=sq, pow=pow, digits=digits, iter=iter))
}

# S^2:
DIFF = function(R) 4*(R[1]*R[3] - 2*R[4])*R[2] - 8*R[1]*R[3]*R[4];
# S^1:
DIFF = function(R) 4*R[1]^3*R[3]*R[4] + 4*R[1]*R[3]^3 + 2*(R[1]^2*R[4] + R[3]^2)*R[2];
# S^0:
DIFF = function(R) - R[2]^3 + 4*R[2]^2*R[4] - 2*R[1]^3*R[3]^3 +
	+ 2*R[1]^2*R[3]^2*R[4] - R[1]^4*R[4]^2 - R[3]^4 + 8*R[1]*R[2]*R[3]*R[4] +
	- 4*R[1]*R[2]^2*R[3] - 5*R[1]^2*R[2]*R[3]^2;

###
p1 = toPoly.pm("(xs^2 - 2*x13)*(x13*(S - xs)^2 - 2*R4) - R2*x13")
p2 = toPoly.pm("x13^2*(S - xs) - R3*x13 + R4*xs")
#
pR = solve.pm(p2, p1, "x13")
str(pR)
pR = div.pm(pR$Rez, toPoly.pm("xs^2 - 2*S*xs + S^2"), "xs")
pR$Rez = sort.pm(pR$Rez, "xs", xn2=c("S", "R4", "R3"))
print.pm(pR$Rez, lead="xs")
print.coeff(pR$Rez, "xs")


########################
########################

####################
### Type: E121a  ###
### Order 1      ###
####################

x1 + x2 + x3 + x4 - R1 # = 0
x1*x2^2*x3 + x2*x3^2*x4 + x3*x4^2*x1 + x4*x1^2*x2 - R2 # = 0
x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4 - R3 # = 0
x1*x2*x3*x4 - R4 # = 0


### [old]

### Simple version:
# - computes only E2 (via x13);
# - used to derive the proper Formula;
coeff.X13.E121aP1 = function(R) {
	S = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	coeff = c((4*R4 + R2), - (R4*S^2 + R3^2), 4*R3*R4*S - 8*R4^2 - 2*R4*R2,
		- R4*(R4*S^2 + R3^2), R4^2*(4*R4 + R2) );
	return(coeff);
}
e2.old = function(R, debug=TRUE) {
	coeff = coeff.X13.E121aP1(R);
	x13 = roots(coeff);
	if(debug) print(x13);
	x24 = R[4] / x13;
	S = R[1]; R3 = R[3];
	xs = - x13*R3 + x13^2*S;
	xs = xs / (x13^2 - R[4]);
	E2a = xs * (S - xs);
	E2  = E2a + x13 + x24;
	return(E2);
}

### Solver [old]
coeff.S4Ht.E121aP1.old = function(R) {
	p = replace.pm(pR$Rez, c(S = R[1], E121a = R[2], E3 = R[3], E4 = R[4]));
	coeff = coef.pm(p, "E2");
	return(coeff);
}
E22a.E121aP1.old = function(R, E2) {
	len = length(E2);
	# c(E3, S, E4, E2, E121a)
	e22a = sapply(seq(len),
		function(id) eval.pm(pR$x0, list(E3=R[3], S=R[1], E4=R[4], E2=E2[id], E121a=R[2])));
	eDiv = sapply(seq(len),
		function(id) eval.pm(pR$div, list(E3=R[3], S=R[1], E4=R[4], E2=E2[id], E121a=R[2])));
	return( e22a / eDiv);
}
solve.S4HtM.E121aP1.old = function(R, sort=TRUE, all.sol=FALSE, debug=TRUE) {
	coeff = coeff.S4Ht.E121aP1.old(R);
	E2 = roots(coeff);
	if(debug) print(E2);
	E22a = E22a.E121aP1(R, E2);
	len  = length(E2);
	#
	sol = lapply(seq(len), function(id) {
		RS = R; RS[2] = E22a[id];
		solve.S4HtM.Ord2Base(RS, E2[id], sort=sort, all.sol=all.sol)
	})
	sol = do.call(rbind, sol);
	if(sort) sol = sort.sol(sol, ncol=1, useRe=TRUE, mod.first=FALSE);
	return(sol);
}
### [other older solvers]
solve.x3.S4HtM.E121P1.old = function(R, E2, x1) {
	# pR needs to be computed!
	# [see below]
	S = R[1]; E121 = R[2]; E3 = R[3]; R4 = R[4]; # E4
	E2d = E2 - x1*(S-x1);
	x0 = eval.pm(pR[[2]]$x0, list(x1=x1, S=S, E121=E121, E2d=E2d, E3=E3, R4=R4));
	div = eval.pm(pR[[2]]$div, list(x1=x1, S=S, E121=E121, E2d=E2d, E3=E3, R4=R4));
	return(x0 / div);
}
solve.S4Ht.E2aE2 = function(id, R, E2) {
	# [old] [not used anymore]
	S = R[1]; E3 = R[3]; E4 = R[4]; E2 = E2[id];
	# E2a*E2^2 - (S*E3 + 2*E2a^2)*E2 +
	#	+ S^2*E4 + S*E2a*E3 + E2a^3 - 4*E2a*E4 + E3^2
	coeff = c(1, - 2*E2, E2^2 + E3*S - 4*E4,
		- S*E3*E2 + S^2*E4 + E3^2);
	return(roots(coeff));
}


### Derivation:

# non-efficient!
pE121a = toPoly.pm("E2a^2 - E22a - 2*E121a - 4*E4")
pE21 = polyE2Ord1(); # E2, E2a
pE22 = polyE2Ord2(); # E2, E22a

pR1 = solve.pm(pE121a, pE21, "E2a");
pR = solve.pm(pR1$Rez, pE22, "E22a");

div = toPoly.pm("E3^2*S^2 - 8*E3*E4*S - 6*E3*E121a*S + 16*E4^2 + 24*E4*E121a + 9*E121a^2");
pR$Rez = div.pm(pR$Rez, div, "S")$Rez
# pR$Rez$coeff = - pR$Rez$coeff;
pR$Rez = sort.pm(pR$Rez, xn="E2", xn2=c("E4", "E3", "S"))

# TODO:
# Note:
# - div: from 1200 Monomials => 598 Monomials;
# - E2^12, S^14;
# - correct roots: should be only E2^2;
# print.pm(pR$Rez, lead="E2")
# print.coeff(pR$Rez, "E2")
pLead = pR$Rez[pR$Rez$E2 == 12, ]
pLead$E2 = NULL
xgcd = gcd.pm(pLead); pLead$coeff = pLead$coeff / xgcd;
pLead = sort.pm(pLead, xn="S", xn2=c("E4", "E3"))
pLead$E4 = pLead$E4 - min(pLead$E4)
print.pm(pLead, lead="S")


### Lead: E2^12
S = R[1]; E121a = R[2]; E3 = R[3]; E4 = R[4];
E4^2 * (4*E4 + E121a) * (E3^2*S^2 - 8*E3*E4*S - 6*E3*E121a*S + 16*E4^2 + 24*E4*E121a + 9*E121a^2)


###
R = c(7,-1,2,1)
sol = solve.S4HtM.E121aP1(R)
test.S4HtMixed.En3(sol, n=1, nE=c(1,2,1))

apply(sol[c(17:20, 36,37, 38,39), ], 1, e2.f)
poly.calc(apply(sol[c(36,39), ], 1, e2.f)) * 4 * (4*R[4] + R[2]) * R[4]


### [old]
R = c(7,-1,2,1)
e2 = roots(coeff.S4Ht.E121aP1(R))
pr = poly.calc(e2[c(1,2)]);
(pr[2] - round(pr[2])) * 4 * (4*R[4] + R[2]) * R[4]

###
R = c(1,1,1,1) # (1,1,1,1) is different;
poly.calc(sort(e2.old(R, debug=F))[c(1,4)]) * 4 * (4*R[4] + R[2]) * R[4]

which.coeff(R <- c(1, sqrt(2),1,1), id=1, sq=2)

### E2:
4*(4*R[4] + R[2])*R[4]*E2^2 - (R[2] + 8*R[4])*(R[1]^2*R[4] + R[3]^2)*E2 +
	+ R[1]^4*R[4]^2 + R[1]*R[2]^2*R[3] - R[2]^3 + R[3]^4 - 4*R[2]^2*R[4] + 2*R[1]^2*R[3]^2*R[4];


solveBase = function(R, debug=FALSE) {
	mult = 4 * (4*R[4] + R[2]) * R[4];
	r = sort(e2.old(R, debug=debug));
	if(debug) print(r);
	round0(poly.calc(r[c(1,4)])) * mult;
}
which.coeff = which.coeff.gen(FUN=solveBase);


eval.pm(pR$Rez[pR$Rez$E2 == 0, ], c(S=R[1], E121a=R[2], E3=R[3], E4=R[4], E2=1))

####
p1 = toPoly.pm("x13^2*(S-xs) + R4*xs - R3*x13")
p2 = toPoly.pm("x13^2*(S-xs)^2 - 2*R4*x13 + R4*(xs^2 - 2*x13) - R2*x13")
#
pR2 = solve.pm(p1, p2, "xs")
str(pR2)
pR2$Rez = sort.pm(pR2$Rez, xn="x13", xn2 = c("S", "R4", "R3"))
print.pm(pR2$Rez, lead="x13")
print.coeff(pR2$Rez, "x13")


### Robust:
# - pR needed for the robust computation of the solutions;
pE2d = toPoly.pm("E2 - x1*(S - x1)");
pE2  = toPoly.pm("x1*x3^2*(S-xs) + R4 - E2d*x1*x3");
pE3  = toPoly.pm("(x1*x3)^2*(S - xs) + R4*xs - E3*x1*x3");
pxs  = toPoly.pm("x1 + x3");
pE2x = replace.pm(pE2, pxs, "xs");
pE3x = replace.pm(pE3, pxs, "xs"); # redundant;
# How to derive short/compact expression for x3?
# [DONE: differently]
# [old & complicated] =>
pE121old = toPoly.pm("(x1*x3)^2 * (S^2 + xs^2 - 2*S*xs) - 2*R4*x1*x3 + R4*xs^2 - 2*R4*x1*x3 - x1*x3*E121");
pE121 = toPoly.pm("(E3*x1*x3 - R4*xs) * (S-xs) - 2*R4*x1*x3 + R4*xs^2 - 2*R4*x1*x3 - x1*x3*E121");
pE121x = replace.pm(pE121, pxs, "xs");
pR = solve.pm(pE2x, pE121x, "x3")

### [old] complicated
# pR = solve.lpm(pE2, pE3, pE121, xn=c("xs", "x3"))
# pR[[2]]$x0$coeff = - pR[[2]]$x0$coeff; pR[[2]]$div$coeff = - pR[[2]]$div$coeff;
# print.pm(pR[[2]]$x0, lead="S") # 134 Monomials;

