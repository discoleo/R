########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S4:
### Mixed Type with Resonances
###
### draft v.0.1d


### Heterogeneous Symmetric
### Polynomial Systems: 4 Variables
### Mixed: Hetero + Symmetric

### Example:
# x1^p*x2^n*x3^n + x2^p*x3^n*x4^n + x3^p*x4^n*x1^n + x4^p*x1^n*x2^n = R1
# x1^k + x2^k + x^3^k + x^4^k = R2
# ...


###############

###############
### History ###
###############


### draft v.0.1d:
# - S4MHt311 + P5:
#   robust computation of x2;
### draft v.0.1c - v.0.1c-cases:
# - another example with Resonances:
#   x1^3*x2*x3 + b*x4^5 = Ru;
# - special cases: all equal vs other cases;
### draft v.0.1b:
# - [started] classic roots;
# - more entangled roots; [v.0.1b-more]
### draft v.0.1a:
# - Proof of concept:
#   Rotations Entangled with Roots of Unity;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions

### Rotate roots & Entangle with Roots of Unity
rotate.roots = function(x, n=5, r.m = NULL) {
	x0 = x;
	# Roots of unity
	m = unity(n, all=TRUE);
	if(is.matrix(x0) && nrow(x0) > 1) {
		x = lapply(seq(nrow(x0)), function(id1)
			t(sapply(seq_along(m), function(id2) x[id1,]*m[id2])));
		x = do.call(rbind, x);
	} else {
		x = sapply(seq_along(m), function(id) x*m[id]);
		x = t(x);
	}
	if(is.null(r.m)) {
		r.m = matrix(
			c(1,0,2,3,  1,2,0,4,  3,2,4,0,  1,4,3,0,  1,3,4,2), ncol=5
		);
	}
	if(is.matrix(x0) && nrow(x0) > 1) {
		sol = lapply(seq(nrow(x0)), function(id)
			t(apply(r.m, 2, function(pow) x0[id,] * m[2]^pow)));
		sol = do.call(rbind, sol);
	} else {
		sol = apply(r.m, 2, function(pow) x0 * m[2]^pow);
		sol = t(sol);
	}
	x = rbind(x, sol);
	# Rotate / Cycle
	x = rbind(x, x[,c(2,3,4,1)], x[,c(3,4,1,2)], x[,c(4,1,2,3)]);
	rownames(x) = NULL
	return(x)
}
### E
calc.E = function(x) {
	if(is.matrix(x) && nrow(x) > 1) {
		sol = t(apply(x, 1, calc.E));
		colnames(sol) = c("S", "E2", "E3", "E4");
		return(sol);
	}
	S = sum(x); E4 = prod(x);
	E2 = x[1]*(S - x[1]) + x[2]*(x[3]+x[4]) + x[3]*x[4];
	E3 = E4 * sum(1/x);
	E = cbind(S=S, E2=E2, E3=E3, E4=E4);
	return(E);
}
reduce.E = function(x1, E) {
	if( ! is.matrix(E)) E = matrix(E, nrow=1);
	S = E[,1]; E2 = E[,2]; E4 = E[,4];
	S  = S - x1;
	E2 = E2 - x1*S;
	E3 = E4 / x1;
	return(cbind(S, E2, E3));
}


#####################
#####################

##################
### Rotations: ###
### x1^n*x2*x3 ###
##################

### Order: 3+1+1
# x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 = R1
# x1^5 + x2^5 + x3^5 + x4^5 = R2
# x1^10 + x2^10 + x3^10 + x4^10 = R3
# (x1*x2)^5 + (x2*x3)^5 + (x3*x4)^5 + (x4*x1)^5 = R4


### Classic Solution:

# a toy model for computing robust roots;
solve.byE.S4M5 = function(E, R, mpow=0, debug=TRUE) {
	if( ! is.matrix(E)) E = matrix(E, nrow=1);
	x1 = sapply(seq(nrow(E)), function(id) roots(c(1, -E[id,1], E[id,2], -E[id,3], E[id,4])));
	if(debug) print(x1);
	S  = rep(E[,1], each=4);
	E2 = rep(E[,2], each=4);
	E4 = rep(E[,4], each=4);
	S  = S - x1;
	E2 = E2 - x1*S;
	E3 = E4 / x1;
	### TODO:
	# - robust computation of x2;
	# - use: solve.x2.byE.S4MHt311() in file:
	#   Poly.System.S4.HtMixed.Resonances.Helper.R;
	### non-robust:
	A = roots(c(1, -R[2], R[4]));
	A = rep(A, each=length(x1));
	B = R[2] - A;
	x1 = rep(x1, 2); # NOT each;
	x3p5 = A - x1^5;
	x3 = rootn(x3p5, n=5);
	if(mpow != 0) {m = unity(5, FALSE); x3 = x3 * m^mpow;}
	#
	S  = rep(S, 2);
	E2 = rep(E2, 2);
	E3 = rep(E3, 2);
	s24 = S - x3; x24 = E3 / x3;
	x2 = x1*s24*x3^3 - x24^2*(x1+x3) + x24*x1*s24^2 - R[1];
	x2 = - x2 / ((x1 - x3)*(x1*x3*(x1 + x3) - x24*s24));
	x4 = s24 - x2;
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	return(sol);
}
solve.byx1.S4M5.classic = function(x1, R) {
	A = roots(c(1, -R[2], R[4]));
	A = rep(A, each=length(x1));
	B = R[2] - A;
	x1 = rep(x1, each=2);
	x3p5  = A - x1^5;
	x24p5 = (R[2]^2 - R[3] - 2*R[4]) / 2 - x1^5*x3p5;
	# non-robust
	# x2^10 - B*x2^5 + x24p5 = 0
	len = length(x3p5);
	x2p5 = sapply(seq(len), function(id) roots(c(1, - B[id], x24p5[id])));
	x2p5 = as.vector(x2p5);
	# repeat x1
	x1 = rep(x1, each=2);
	B = rep(B, each=2);
	x3p5 = rep(x3p5, each=2);
	x4p5 = B - x2p5;
	# TODO: ???
	m = unity(5, all=TRUE);
	x = lapply(seq_along(x1), function(id) {
		x2 = rootn(x2p5[id], 5) * m;
		x3 = rootn(x3p5[id], 5) * m;
		x4 = rootn(x4p5[id], 5) * m;
		expand.grid(x1[id], x2, x3, x4);
	})
	sol = do.call(rbind, x);
	names(sol) = paste0("x", seq(4));
	return(sol);
}
test.R1 = function(x) {
	sum(x[1]^3*x[2]*x[3], x[2]^3*x[3]*x[4], x[3]^3*x[4]*x[1], x[4]^3*x[1]*x[2]);
}

### Test:

### by E:
R = c(1, -2, -1, 3);
# x = see below;
S = sum(x); E4 = prod(x);
E2 = x[1]*(S - x[1]) + x[2]*(x[3]+x[4]) + x[3]*x[4];
E3 = E4 * sum(1/x);
E = cbind(S=S, E2=E2, E3=E3, E4=E4);
sol = solve.byE.S4M5(E, R, mpow=1);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];


### Classic by x1
# x1 & R: see below;
x = solve.byx1.S4M5.classic(x1, R);
R1r = apply(x, 1, test.R1);
R1r = round(R1r, 2);
# only 2*5 = 10 root-tuples are real!
x[R1r == R[1], ]


### Test
x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 # - R1
x1^5 + x2^5 + x3^5 + x4^5 # - R2
x1^10 + x2^10 + x3^10 + x4^10 # - R3
(x1*x2)^5 + (x2*x3)^5 + (x3*x4)^5 + (x4*x1)^5 # - R4


### Derivation:

### let:
A = x1^5 + x3^5;
B = x2^5 + x4^5;
# Eq 2 =>
A + B - R2 # = 0
# Eq 4 =>
A*B - R4 # = 0
# =>
A^2 - R2*A + R4 # = 0

### (Eq 2)^2 - Eq 3 - 2*Eq 4 =>
2*(x1*x3)^5 + 2*(x2*x4)^5 - R2^2 + R3 + 2*R4 # = 0
# (x1*x3)^5 + (x2*x4)^5 = (R2^2 - R3 - 2*R4) / 2;
# if (x1*x3)^5 is known =>
x2^10 - B*x2^5 + (R2^2 - R3 - 2*R4) / 2 - (x1*x3)^5 # = 0

### Eq 2 =>
# x1^5 + x3^5 = R2 - (x2^5 + x4^5) # sq =>
x1^10 + x3^10 + 2*(x1*x3)^5 - (x2^10 + x4^10 + 2*(x2*x4)^5 - 2*R2*(x2^5 + x4^5) + R2^2) # = 0
x2^10 + x4^10 + 2*(x2*x4)^5 - (x1^10 + x3^10 + 2*(x1*x3)^5 - 2*R2*(x1^5 + x3^5) + R2^2) # = 0
### Sum => [redundant]


### Debug:
R = c(1, -2, -1, 3);
R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
x1 =  0.0451568676 + 0.1440038706i;
x2 =  0.9916018740 + 0.5052387172i;
x3 = -0.1370521599 - 1.1076829770i;
x4 =  0.7725898004 + 0.1223347600i;
x = c(x1, x2, x3, x4);


x = rotate.roots(x, n=5);
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];


p1 = toPoly.pm("x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 - R1")
p2 = toPoly.pm("x2^5 - x25");
p4 = toPoly.pm("x4^5 - x45");
p45 = data.frame(x45=1, coeff=1);
pR = solve.pm(p1, p2, "x2")
pR$Rez = replace.pm(pR$Rez, p45, "x4", 5);
# Memory overflow!
pR4 = solve.pm(pR$Rez, p4, "x4")


p1 = toPoly.pm("x1^3*x2*x3 + x2^2*x3*x24 + x3^3*x4*x1 + x4^2*x1*x24 - R1")
p1 = replace.pm(p1, toPoly.pm("s24 - x2"), "x4", 1);
p2 = toPoly.pm("x2^2 - s24*x2 + x24");
pR = solve.pm(p1, p2, "x2")


p1 = toPoly.pm("x1^3*x2*x3 + E3*x2^2 + x3^3*x4*x1 + x4^3*x1*x2 - R1")
p1 = replace.pm(p1, toPoly.pm("S - x2 - x3"), "x4", 1);
pS = toPoly.pm("x2^2*x3 + x2*x3^2 + E3 - x2*x3*S");
p2 = toPoly.pm("x2^2*S - x2^3 + E3 - E2*x2");
#
pR1 = solve.pm(p1, pS, "x3");
pR  = solve.pm(pR1$Rez, p2, "x2");
str(pR)
# - polynomial with 7039 monomials;
# - robust formula for x2: ~300 / 366 monomials;

### Robust computation of x2:
# - see solve.x2.byE.S4MHt311() in file:
#   Poly.System.S4.HtMixed.Resonances.Helper.R;


#############################
#############################

######################
### Mixed: 3+1+1   ###
### & with Unknown ###
######################

# x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 = R1
# x1^3*x2*x3 + b*x4^5 = Ru
# x2^3*x3*x4 + b*x1^5 = Ru
# x3^3*x4*x1 + b*x2^5 = Ru
# x4^3*x1*x2 + b*x3^5 = Ru
# where Ru = unknown


### Solution:

### Case 1:
# x1 = x2 = x3 = x4
# - trivial case;
x1^5 - R1 / 4 # = 0

### Case 2:
# x1 = x3, x2 = x4;
# x1^4*x2 + x2^4*x1 = R1/2
# x1^4*x2 + b*x2^5 = Ru
# x2^4*x1 + b*x1^5 = Ru

### Case 3:
# - all distinct;


### Sol Case 2:

### Diff =>
b*(x1^5 - x2^5) - x1*x2*(x1^3 - x2^3) # = 0
### x1 ! = x2 =>
b*(x1^4 + x2^4 + x1*x2*(x1^2+x2^2) + (x1*x2)^2) - x1*x2*(x1^2 + x2^2 + x1*x2) # = 0
# S = (x1 + x2);
b*(S^4 - 4*x1*x2*S^2 + x1*x2*(S^2 - 2*x1*x2) + 3*(x1*x2)^2) - x1*x2*(S^2 - x1*x2) # = 0
b*(S^4 - 3*x1*x2*S^2 + (x1*x2)^2) - x1*x2*(S^2 - x1*x2) # = 0
b*S^4 - (3*b+1)*x1*x2*S^2 + (b+1)*(x1*x2)^2 # = 0

### Eq 1 (the "Sum") =>
x1*x2*(S^3 - 3*x1*x2*S) - R1/2 # = 0
x1*x2*S^3 - 3*x1^2*x2^2*S - R1/2 # = 0

### Derivation:
p1 = toPoly.pm("b*S^4 - 3*b*x12*S^2 - x12*S^2 + b*x12^2 + x12^2")
p2 = toPoly.pm("x12*S^3 - 3*x12^2*S - R1/2")
pR = solve.pm(p1, p2, "x12")
pR$Rez$coeff = -4 * pR$Rez$coeff;
print.p(pR$Rez, "S")

### Eq S:
12*(b^2 - 2*b)*S^10 + 12*R1*(1 + 4*b + 9*b^2)*S^5 + 3*R1^2*(b+1)^2

### Solver:
solve.S4M311.Case2 = function(R, b, debug=TRUE, all.rotated=FALSE) {
	coeff = c(12*(b^2 - 2*b), 12*R[1]*(1 + 4*b + 9*b^2), 3*R[1]^2*(b+1)^2);
	S = roots(coeff);
	S = rootn(S, 5);
	m = unity(5, all=TRUE);
	if( ! all.rotated) S = sapply(S, function(S) S*m);
	S = as.vector(S);
	if(debug) print(S);
	xy = - 3*b*S^5 + R[1]*(b+1)/2;
	xy = xy / -(2*(4*b+1)*S^3);
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = (S - xy.d) / 2;
	sol = cbind(x1=x, x2=y, x3=x, x4=y);
	if(all.rotated) {
		sol = rotate.roots(sol, n=5);
	} else sol = rbind(sol, sol[,c(2,1,2,1)]);
	return(sol);
}

### Example:
R = 2
b = 3
sol = solve.S4M311.Case2(R, b, all.rotated=TRUE);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];


### Test
x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 # - R1
tmp1 = x1^3*x2*x3 + b*x4^5
tmp2 = x2^3*x3*x4 + b*x1^5
tmp3 = x3^3*x4*x1 + b*x2^5
tmp4 = x4^3*x1*x2 + b*x3^5
cbind(tmp1, tmp2, tmp3, tmp4)


### Debug:
R = 2; b = 3;
x1 = -0.9511358311 - 0.3938106618i;
x2 =  1.0126064890 - 0.2503814432i;
x3 =  0.4904683198 + 0.9794972522i;
x4 =  0.9652251325 - 0.1464569610i;
x = c(x1, x2, x3, x4);


# 40 roots derived from a base-root
x = rotate.roots(x, n=5);
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];

### another set of roots
x1 = x2 = x3 = x4 = rootn(R[1]/4, 5);
x = c(x1, x2, x3, x4);
x = rotate.roots(x, n=5);
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];


##################
##################
##################

##################
### Resonances ###
##################

### Order: 3+1+1
### 4 Variables
# i = 4; p = c(3,1,1)
# p = Divisors(3^3*5^2);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(15);
# - the computed ones do NOT work;
# - but there are quasi-non-trivial solutions for p = 5;
### Ex: p = 5 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^0, m^4);
# (1,2,0), (2,0,4), (0,4,1), (4,1,2)
# (1,3,4), (3,4,2), (4,2,1), (2,1,3)
# (1,4,3), (4,3,0), (3,0,1), (0,1,4)
# (1,0,2), (0,2,3), (2,3,1), (3,1,0)
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)
