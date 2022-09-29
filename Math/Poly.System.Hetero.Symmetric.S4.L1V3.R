########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
### Leading 1, NL Mixed: V3
###
### draft v.0.1f-Special


### Type L1 NLm V3
# - Leading: 1 var;
# - Non-Leading / Side-Chain: Mixed, 3 vars;

### Subtypes:
# V3a: x1^n + b*x1*x2*x3 = R
# V3b: x1^n + b*x2*x3*x4 = R
# V3c: x1^n + b*x1*x2*x4 = R


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


test.S4Ht.V3b = function(sol, b, R=NULL, n=2) {
	x1 = sol[,1]; x2 = sol[,2];
	x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + b*x2*x3*x4;
	err2 = x2^n + b*x1*x3*x4;
	err3 = x3^n + b*x1*x2*x4;
	err4 = x4^n + b*x1*x2*x3;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		err = err - R;
	}
	notNAN = apply(err, 2, function(x) ! any(is.nan(x)));
	err[ , notNAN] = round0(err[ , notNAN]);
	return(err);
}

####################

### History

### v.0.1a
# - moved specific code from file:
#   Poly.System.Hetero.Symmetric.S4.R


#############################
#############################

### Theory

#############
### Type V3a:
# TODO


#############
### Type V3b:

# - System is easily transformed to a system in {S, E2, E3, E4};
#  -- the full Diff simplifies tremendously the system;
#     (but works only for higher powers and if NO additional terms)
# - alternative: C2-Decomposition;

### Eq 1: Sum =>
# S[n] + b*E3 - 4*R = 0

### Eq 2: Sum(x1*...) =>
# S[n+1] + 4*b*E4 - R*S = 0

### Eq 3: Sum(x1^2*...) =>
# S[n+2] + b*E4*S - R*(S^2 - 2*E2) = 0

### Eq 4: Sum(x1^3*...) =>
# S[n+3] + b*E4*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) = 0

### Alternative: Diff Eqs
# - simplify tremendously the solution;
# - but do NOT work with small powers / certain powers (?);

### Eq 3: Diff(Eq 1 - Eq 3) =>
(x1^n - x3^n) - b*x2*x4*(x1 - x3) # = 0
(x2^n - x4^n) - b*x1*x3*(x2 - x4) # = 0

### Eq 4: full Diff =>
(x1^n - x2^n) - b*x3*x4*(x1 - x2) # = 0
# ...
# TODO: more work;


#############
### Type V3c:
# - is trivially decomposable:
(x1^n - x3^n) + b*x2*x4*(x1 - x3) # = 0
(x2^n - x4^n) + b*x1*x3*(x2 - x4) # = 0
(x1^n + x3^n) + b*x2*x4*(x1 + x3) - 2*R # = 0
(x2^n + x4^n) + b*x1*x3*(x2 + x4) - 2*R # = 0
# TODO: Theory & C2-Decomposition;


###############################
###############################

##################
### Type: V3b  ###
##################

###############
### Order 2 ###
###############

### V3b:
### x1^2 + b*x2*x3*x4 = R

### Solution:

### Case 1:
# x1 = x2 = x3 = x4
x^2 + b*x^3 - R # = 0

### Case 2:
# x1 = x2 = x3, but != x4;
x1^2 + b*x1^2*x4 - R # = 0
x4^2 + b*x1^3 - R # = 0

### Case 3:
# x1 = x2 and x3 = x4
# - same as: x1 = x3 and x2 = x4,
#   only the variable index varies slightly;
x1^2 + b*x1*x3^2 - R # = 0
x3^2 + b*x3*x1^2 - R # = 0

### Case 4:
# x1 = x2 and x3 != x4
x1^2 + b*x1*x3*x4 - R # = 0
x3^2 + b*x1^2*x4 - R # = 0
x4^2 + b*x1^2*x3 - R # = 0

### Case 5: x[1:4] all distinct;
# - NO such solutions found!


### Equations:

### Eq 1: Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Eq 2: Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### Eq 3: Sum(x1^2*...) =>
S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S + b*E4*S - 4*E4 # = 0
# Reduction =>
b*E2*E3 - 2*R*E2 + E3*S - 3*b*E4*S - 4*E4 # = 0

### [old]
### Diff Eq[i] - Eq[i+1] =>
# - non-robust for powers < ???;
(x1 - x2)*(x1 + x2 - b*x3*x4) # = 0
(x2 - x3)*(x2 + x3 - b*x1*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum =>
3*S - b*E2 # = 0

### =>
#   b*E2 = 3*S; [does NOT work]
# - b^2*E3 = b*S^2 - 6*S - 4*b*R
# - 4*b^3*E4 = b^2*S^3 - 12*b*S^2 - b^2*R*S + 18*S + 12*b*R

#########

### Eq S:
(b*S + 1) * (b*S^3 + 4*S^2 - 64*R) * (b*S^2 - 2*S - 4*b*R) *
(b^2*S^4 - b*S^3 + 7*S^2 - 2*b^2*R*S^2 - 24*b*R*S + b^2*R^2 - 28*R) # = 0


###########
### Solver:

solve.S4Ht.V3bP2 = function(R, b, debug=TRUE, all=FALSE) {
	### Special Cases: x = 0
	solSp = NULL;
	if(round0(b[1]^2*R - 1) == 0) {
		warning("Special Case: x1 = 0;");
		solSp = solve.S4Ht.V3bP2.Sp0(R, b, all=all);
	}
	### Case: x1 == x2 == x3, but != x4;
	sol3 = solve.S4Ht.V3bP2.CaseX3(R, b);
	
	### Case: x1 == x2, x3 == x4, but x1 != x3;
	sol22 = solve.S4Ht.V3bP2.CaseX22(R, b);
	
	### Case: x1 = x2, x3 != x4; [only 1 real case]
	# contains also the variants: x1 == x2 == x3 != x4;
	sol2x = solve.S4Ht.V3bP2.CaseX2(R, b, debug=debug);
	sol22 = rbind(sol22, sol2x);
	
	### Case: x[1] != x[j]
	# actually: x2 = x3 = x4;
	b1 = b[1];
	coeff = c(b^2, - b, 7 - 2*b^2*R, - 24*b*R, b^2*R^2 - 28*R);
	S = roots(coeff);
	if(debug) print(S);
	#
	E2 = 18*b1^3*S^5 - 114*b1^2*S^4 + 204*b1*S^3 - 18*R*b1^3*S^3 + 336*S^2 +
		+ 156*R*b1^2*S^2 - 1152*R*b1*S - 1344*R + 384*R^2*b1^2;
	div = 62*b1^3*S^3 - 296*b1^2*S^2 + 296*b1*S - 32*R*b1^3*S + 672 - 416*R*b1^2;
	E2 = E2 / div;
	E3 = - (S^2 - 2*E2 - 4*R) / b1;
	E4 = (E2*S + b*E3*S - 3*E3 - 3*R*S) / (4*b1);
	### Robust: x1 is actually linear!
	# - formula valid only for Case: x2 = x3 = x4!
	# - without this simplifying assumption, the formula will become
	#   probably humongous; (e.g. to work also for the (b*S + 1) case)
	# - but there are NO distinct solutions anyway (for Power = 2);
	# - the (b*S + 1) case does NOT work with this formula;
	x1 = - 3*b*S^5 - b*R*S^3 + 81*R*S^2 + 6*b^2*E4*S^2 + 162*b*E4*S +
		+ 27*R^2 + 729*E4 - b^2*R*E4;
	div = - 6*b*S^4 + 27*S^3 + 3*b*R*S^2 + 81*R*S + 3*b^2*E4*S - b*R^2 + 27*b*E4;
	x1 = x1 / div;
	### x2: robust;
	# - Note: x2 = x3 = x4 =>
	#   shortcut possible: x2 = (S - x1) / 3;
	s3 = S - x1;
	p3 = (R - x1^2) / b1;
	e2 = (E3 - p3) / x1;
	x2 = (- R*b*E4 - b*e2*E4 + b*s3^2*E4 - R*p3 - p3*e2);
	x2 = x2 / (b*s3*E4 - R^2 - 2*R*e2 - e2^2 + p3*s3 + R*s3^2);
	# x3:
	s2 = s3 - x2;
	p2 = p3 / x2;
	x3 = (b*E4 - p2*s2) / (R + p2 - s2^2);
	x4 = s2 - x3;
	#
	sol = cbind(x1, x2, x3, x4);
	if( ! is.null(solSp)) sol = rbind(sol, solSp);
	return(list(sol=sol, sol3=sol3, sol22=sol22))
}
solve.S4Ht.V3bP2.CaseX3 = function(R, b) {
	### Case: x1 == x2 == x3, but != x4;
	# - solved using: x & y (not S);
	coeff.3eq = c(b^2, - b, 1, 0, - R);
	x3 = roots(coeff.3eq);
	y = (R - x3^2) / b / x3^2;
	S.3eq = 3*x3 + y;
	sol = list(sol=cbind(x123=x3, x4=y), S=S.3eq);
	return(sol);
}
solve.S4Ht.V3bP2.CaseX22 = function(R, b) {
	### Case: x1 == x3, x2 == x4, x1 != x2;
	# - solved using: S2 = x + y;
	b1 = b[1];
	S2 = roots(c(b1, -1, - b1*R))
	xy2 = S2 / b1;
	x2 = sapply(seq_along(S2), function(id) roots(c(1, -S2[id], xy2[id])))
	y2 = x2[2:1, ];
	x2 = as.vector(x2); y2 = as.vector(y2);
	sol = cbind(x1=x2, x2=x2, x3=y2, x4=y2); # + many permutations;
	return(sol);
}
solve.S4Ht.V3bP2.CaseX2 = function(R, b, debug=TRUE) {
	### Case: x1 == x2, x3 != x4; [only 1 real case]
	# - solved using: S2 = x + y;
	b1 = b[1];
	coeff = c(b1^2, b1, 1 - 2*b1^2*R, -2*b1*R, b1^2*R^2);
	S34 = roots(coeff);
	S34 = c(1/b1, S34); # the real case: (b*S - 1) * P[4]
	if(debug) print(S34);
	p = S34^2 - R;
	x34.d = sqrt(S34^2 - 4*p + 0i); # TODO: +/-;
	x3 = (S34 + x34.d) / 2; x4 = (S34 - x34.d) / 2;
	# robust
	x1 = ((b1*x4-1)*R + x3^2) / (b1^2*p*x4);
	#
	sol = cbind(x1=x1, x2=x1, x3=x3, x4=x4)
	return(sol);
}
solve.S4Ht.V3bP2.Sp0 = function(R, b, all=FALSE) {
	# Check: round0(b[1]^2*R - 1) == 0;
	x2 = rootn(R, 2);
	if(sign(b) != sign(R)) x2 = - x2;
	# Case: x2 = x3 = x4;
	sol = cbind(0, x2, x2, x2);
	# Case: 2 of x[i] distinct;
	sol = rbind(sol, cbind(0, x2, - x2, - x2));
	sol = rbind(sol, cbind(0, - x2, x2, - x2));
	sol = rbind(sol, cbind(0, - x2, - x2, x2));
	# C2-permutations work as well;
	# - all permutations = actually + all cyclic permutations;
	if(all) sol = rbind(sol, sol[ , c(3,1,4,2)]);
	return(sol);
}
test.S4Ht.V3bP2 = function(sol, b, R=NULL) {
	test.S4Ht.V3b(sol, b=b, R=R, n=2);
}

### Examples:

R = -5;
b = 3
sol.all = solve.S4Ht.V3bP2(R=R, b=b)
sol = sol.all$sol;
x1 = sol[,1]; x2 = sol[,2];
x3 = sol[,3]; x4 = sol[,4];

### Test
test.S4Ht.V3bP2(sol, b=b)
test.S4Ht.V3bP2(sol.all$sol22, b=b)


### Ex 2:
R = 2; b = - sqrt(1/R);
sol.all = solve.S4Ht.V3bP2(R=R, b=b)
sol = sol.all$sol;

test.S4Ht.V3bP2(sol, b=b)


### Test: Case 3 equal
x1=x2=x3 = sol.all$sol3$sol[,1];
x4 = sol.all$sol3$sol[,2];
x1^2 + b*x2*x3*x4 # - R
x4^2 + b*x1*x2*x3 # - R

poly.calc(sol$sol3$S) * b^2


### Cases:

### Case: x1 == x2 == x3, but != x4:
# - degenerates to a non-symmetric S2 system;
x^2 + b*x^2*x4 # - R
x4^2 + b*x^3 # - R

b^3*x^7 - b^2*R*x^4 + x^4 - 2*R*x^2 + R^2 # = 0
(x^2 + b*x^3 - R) * (b^2*x^4 - b*x^3 + x^2 - R) #= 0

### Case: x1 == x2
# degenerates to a non-symmetric S3 system:
x^2  + b*x*x3*x4 # - R
x3^2 + b*x^2*x4 # - R
x4^2 + b*x^2*x3 # - R

### SubCase: x1 == x2, x3 == x4
# x1 != x3;
x^2 + b*x*y^2 # - R
y^2 + b*x^2*y # - R
# - exactly solvable;
### Diff
# b*x*y = S2;
### Sum =>
S2^2 - 2*x*y + b*x*y*S2 - 2*R # = 0
b*S2^2 - S2 - b*R # = 0


###################
###################

###############
### Order 3 ###
###############

### V3b:
### x1^3 + b*x2*x3*x4 = R

### Solution:

### Eq 1: Sum =>
S^3 - 3*E2*S + (b+3)*E3 - 4*R # = 0


### Eq 2: Sum(x1*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + 4*b*E4 - R*S # = 0
S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4 + 4*b*E4 - R*S # = 0
# Reduction =>
E2*S^2 - 2*E2^2 + (b-1)*E3*S - 4*(b-1)*E4 - 3*R*S # = 0


### Eq 3: Sum(x1^2*...) =>
(x1^5 + x2^5 + x3^5 + x4^5) + b*E4*S - R*(S^2 - 2*E2) # = 0
S^5 - 5*E2*S^3 + 5*E3*S^2 - R*S^2 + 5*E2^2*S + b*E4*S - 5*E4*S - 5*E2*E3 + 2*R*E2 # = 0
# Reduction =>
2*E2*S^3 - 3*R*S^2 - 2*E3*S^2 + b*E3*S^2 - 5*E2^2*S + 5*E4*S - E4*b*S - 2*E2*R + 5*E2*E3 # = 0
b*E3*S^2 - 3*R*S^2 + E2^2*S - (7*b - 3)*E4*S + 2*R*E2 - 5*E2*E3 # = 0


### Eq 4: Sum(x1^3*...) =>
(x1^6 + x2^6 + x3^6 + x4^6) + b*E4*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
(S^6 - 6*E2*S^4 + 6*E3*S^3 + 9*E2^2*S^2 - 6*E4*S^2 - 12*E2*E3*S - 2*E2^3 + 3*E3^2 + 6*E2*E4) +
	+ b*E4*(S^2 - 2*E2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
# Reduction =>
3*b*E4*S^2 - (b + 1)*E2*E3*S + 3*E3^2 +
	+ (2*b + 2)*E2*E4 + 2*R*E2*S - 3*R*E3 # = 0


### Solver:

# Standard solver:
solve.S4Ht.L1V3bP3 = function(R, b, debug=TRUE) {
	isSpecial = round0(b[1] - c(-3, 1)) == 0;
	if(any(isSpecial)) {
		if(isSpecial[2]) {
			# b = 1;
			return(solve.S4Ht.L1V3bP3.B1(R, b=1));
		}
		stop("Special case: not yet implemented!");
	}
	coeff = coeff.S4Ht.L1V3bP3(R, b=b);
	S = roots(coeff);
	S = rootn(S, 3);
	m = unity(3, all=TRUE);
	S = sapply(S, function(S) S * m);
	if(debug) print(S);
	S = as.vector(S);
	#
	cc.all = coeff0.S4Ht.L1V3bP3(R, b=b);
	E2 = with(cc.all, {
		S3 = S^3;
		E2 = (d0*c1 - d7*c2 + (d4*c1 - d7*c4)*S3 + S3^2*d6*c1)*S^2;
		E2div = - d5*c1*S3^2 + (d7*c3 - d2*c1)*S3 + d7*c0;
		E2 = E2 / E2div;
	});
	E4 = (b - 1)*S^4 - 4*b*E2*S^2 - b*R*S + 13*R*S + 2*b*E2^2 + 6*E2^2;
	E4 = - E4 / (4*(b + 3)*(b - 1));
	E3 = S^3 - 3*E2*S - 4*R;
	E3 = - E3 / (b + 3);
	# TODO: robust;
	len = length(S);
	x1 = sapply(seq(len), function(id) {
		coeff = c(1, -S[id], E2[id], -E3[id], E4[id]);
		roots(coeff);
	})
	x1 = as.vector(x1);
	S  = rep(S, each=4);
	E3 = rep(E3, each=4);
	E4 = rep(E4, each=4);
	#
	s3 = S - x1;
	p3 = (R - x1^3) / b;
	e2 = (E3 - p3) / x1;
	#
	div = p3^2 + e2^3 - 2*p3*e2*s3 + e2*E4*b - s3^2*E4*b - 2*p3*R + 3*e2*s3*R - s3^3*R + R^2;
	x2 = - p3*e2^2 + p3^2*s3 + p3*E4*b - 2*e2*s3*E4*b + s3^3*E4*b - p3*s3*R - E4*b*R;
	x2 = - x2 / div;
	#
	s2 = s3 - x2;
	p2 = p3 / x2;
	x3 = (E4*b + p2^2 - p2*s2^2) / (R + 2*p2*s2 - s2^3);
	x4 = s2 - x3;
	#
	sol = cbind(x1, x2, x3, x4);
	return(sol);
}
coeff.S4Ht.L1V3bP3 = function(R, b) {
	# does NOT include solutions: (x1 = x3, x2 = x4);
	coeff = c(- 1 + 2*b - 2*b^2 + 2*b^3 - 2*b^5 + 2*b^6 - 2*b^7 + b^8,
		(21 + 6*b - 58*b^2 + 14*b^3 + 12*b^4 - 38*b^5 + 90*b^6 - 46*b^7 - b^8)*R,
		(- 363 + 558*b - 166*b^2 - 306*b^3 + 972*b^4 - 1030*b^5 +
			+ 390*b^6 - 54*b^7 - b^8)*R^2,
		(343 - 1862*b + 4186*b^2 - 5006*b^3 + 3400*b^4 - 1298*b^5 +
			+ 262*b^6 - 26*b^7 + b^8)*R^3
		);
	return(coeff);
}
# needed for E2:
coeff0.S4Ht.L1V3bP3 = function(R, b) {
	c0 = 8*(b + 3)*(b - 1)*(b - 7)*R; # (168 - 136*b - 40*b^2 + 8*b^3)*R;
	c1 = (90 - 18*b + 38*b^2 + 18*b^3);
	c2 = - (9 - 159*b - 45*b^2 + 3*b^3)*R;
	c3 = - (60 - 40*b + 28*b^2 + 16*b^3);
	c4 = 9 - 15*b + 3*b^2 + 3*b^3;
	#
	d0 = 144*R^2 - 240*b*R^2 + 48*b^2*R^2 + 48*b^3*R^2;
	d1 = 108 + 216*b + 144*b^2 + 40*b^3 + 4*b^4;
	d2 = (846 - 192*b + 28*b^2 + 80*b^3 + 6*b^4)*R;
	d3 = 216 - 198*b + 30*b^2 + 70*b^3 + 10*b^4;
	d4 = (- 180 + 507*b + 243*b^2 + 9*b^3 - 3*b^4)*R;
	d5 = - 198 + 156*b - 52*b^2 - 84*b^3 - 14*b^4;
	d6 = 36 - 51*b - 3*b^2 + 15*b^3 + 3*b^4;
	### Reduction:
	d6 = c1*d6; d5 = c1*d5 - c4*d1;
	d4 = c1*d4; d3 = c1*d3 - c3*d1;
	d2 = c1*d2 - c2*d1;
	d0 = c1*d0; d7 = - c0*d1;
	d1 = c1*d1; # Last!
	#
	d7 = c1*d7; d0 = c1*d0;
	d5 = c1*d5 - d3*c3; d2 = c1*d2 - d3*c0;
	d6 = c1*d6 - d3*c4; d4 = c1*d4 - d3*c2;
	c.l = list(c0=c0, c1=c1, c2=c2, c3=c3, c4=c4);
	d.l = list(d0=d0, d2=d2, d4=d4, d5=d5, d6=d6, d7=d7);
	return(c(c.l, d.l));
}

### Special Cases:
solve.S4Ht.L1V3bP3.Y3 = function(R, b, debug=TRUE, all=FALSE) {
	if(any(round0(b[1] - c(-3)) == 0)) {
		stop("Special case for b: Not yet implemented!");
	}
	coeff = coeff.S4Ht.L1V3bP3.Y3(R, b=b);
	S = roots(coeff);
	S = rootn(S, 3);
	m = unity(3, all=TRUE);
	S = sapply(S, function(S) S * m);
	if(debug) print(S);
	S = as.vector(S);
	#
	y = (- 9 + 54*b - 108*b^2 + 82*b^3 - 3*b^4)*S^4 +
		+ (- 477 + 2124*b - 1998*b^2 - 244*b^3 + 3*b^4)*R*S;
	ydiv = (- 54 + 270*b - 406*b^2 + 258*b^3 - 36*b^4)*S^3 +
		+ (- 756 + 3888*b - 5048*b^2 + 560*b^3 + 12*b^4)*R;
	y = y / ydiv;
	x = S - 3*y;
	#
	sol = cbind(x1=x, x2=y, x3=y, x4=y);
	# just one permutation;
	if(all) sol = rbind(sol, sol(, c(4,1,2,3)));
	return(sol);
}
coeff.S4Ht.L1V3bP3.Y3 = function(R, b) {
	# Excluded: solutions of type (x1 = x3, x2 = x4);
	# (i.e. all 4 equal)
	coeff = c(- (b - 1)*(b^2 + 1),
		2*(b^3 + 21*b^2 - 16*b - 10)*R,
		- (b - 7)^3*R^2
	);
	return(coeff);
}
### b == 1
solve.S4Ht.L1V3bP3.B1 = function(R, b=1, debug=TRUE, all=FALSE) {
	if(round0(b - 1) != 0) stop("Only Special case: b == 1!");
	### Case: x2=x3=x4 = y;
	x1 = 0;
	y = rootn(R, n=3);
	m = unity(3, all=TRUE);
	y = y*m;
	if(debug) print(y);
	sol = cbind(x1=x1, x2=y, x3=y, x4=y);
	### Case: S2Ht
	# TODO;
	### All:
	# - only 1 cyclic permutation added;
	if(all) sol = rbind(sol, sol[ , c(4,1,2,3)]);
	return(sol);
}
solve.S4Ht.L1V3bP3.B1.FALSE = function(R, b=1, debug=TRUE, all=FALSE) {
	# ONLY FALSE solutions!
	warning("False solutions!");
	if(round0(b - 1) != 0) stop("Only Special case: b == 1!");
	# 3*s1^6 - s1^5 + 2*R*s1^3 - 4*R^2;
	s1 = roots(c(3, -1, 0, 2*R,0,0, -4*R^2));
	if(debug) print(s1);
	p1 = (s1^3 + s1^2 - 2*R) / (4*s1);
	#
	y = R / s1^2;
	#
	len = length(s1);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -s1[id], p1[id]));
	})
	x1 = as.vector(x1);
	x2 = s1 - x1;
	y  = rep(y, each=2);
	#
	sol = cbind(x1=x1, x2=x2, x3=y, x4=y);
	# only the C2-permutation is added;
	if(all) sol = rbind(sol, sol[, c(2,1,4,3)]);
	return(sol);
}

### Test:
test.S4Ht.L1V3bP3 = function(sol, b, R=NULL) {
	err = test.S4Ht.V3b(sol, b=b, R=R, n=3);
	return(err);
}

### Examples:

###
R = -5
b = 4
sol = solve.S4Ht.L1V3bP3(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### [many FALSE solutions]
err = test.S4Ht.L1V3bP3(sol, b=b)

err = x4^3 + b*x1*x2*x3;
isR = round(err - R, 3) == 0
err[isR]
sol = sol[isR, ]; sol;


### [Case 3*y] Same set:
sol3 = solve.S4Ht.L1V3bP3.Y3(R, b)
test.S4Ht.L1V3bP3(sol3, b=b)


### Ex 2: Special Case
R = -5
b = -3
# TODO
sol3 = solve.S4Ht.L1V3bP3.Y3(R, b)
test.S4Ht.L1V3bP3(sol3, b=b)


### Ex 3: Special Case
R = 4
b = 7
sol3 = solve.S4Ht.L1V3bP3.Y3(R, b)
test.S4Ht.L1V3bP3(sol3, b=b)


### Ex 4: Special Case
R = 3
b = 1
sol = solve.S4Ht.L1V3bP3(R, b)
test.S4Ht.L1V3bP3(sol, b=b)


#################
### Alternatives:

# - NOT possible for power = 3;
# - create only FALSE solutions;

### Eq 3: Diff(Eq 1 - Eq 3) =>
x1^3 - x3^3 - b*x2*x4*(x1 - x3) # = 0
x2^3 - x4^3 - b*x1*x3*(x2 - x4) # = 0
# Case: x1 != x3; x2 ! = x4;
s1^2 - p1 - b*p2 # = 0
s2^2 - p2 - b*p1 # = 0
# Sum =>
S^2 - 2*ps - (b+1)*sp # = 0
S^2 - 2*E2 - (b-1)*sp # = 0
S^2 - 2*E2 - bd*sp # = 0
# =>
S^6 - (bd + 6)*E2*S^4 + bd^2*E3*S^3 + 4*(bd + 3)*E2^2*S^2 - bd^2*(bd + 4)*E4*S^2 +
	- 2*bd^2*E2*E3*S - 4*(bd + 2)*E2^3 + 4*(2*bd^2 + bd^3)*E2*E4 - bd^3*E3^2 # = 0
# Reduction =>
bd^2*(bd + 4)*E4*S^2 +
	- bd*(bd + 1)*E2*E3*S + (bd - 9)*R*E2*S +
	+ 2*(bd + 1)*E2^3 - 4*bd*(bd^2 + 3*bd + 3)*E2*E4 +
	+ (2*bd^3 + 3*bd^2 - 8*bd - 16)*E3^2 - 4*(bd^2 - 2*bd - 8)*R*E3 - 16*R^2 # = 0


### Eq 4: Diff all
3*S^2 - 5*E2 - b*E2 # = 0


### Debug:
R = 4; b = -3; bd = b - 1;
# Note: x2 == x4, which breaks Eq 3 above (temporarily);
x1 =  1.5329574364 - 0.2097804173i;
x2 =  0.2924017738 + 0.5064547285i;
x3 = -0.9481538887 + 1.2226898741i;
x4 =  0.2924017738 + 0.5064547285i;

#
R = -5; b = 4; bd = b - 1;
x1 = 0.9722777968 - 1.4936414930i;
x2 = -0.3889111187 + 0i;
x3 = 0.9722777968 + 1.4936414930i;
x4 = -0.3889111187 + 0i;

x = c(x1,x2,x3,x4)
s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
sp = p1 + p2; ps = s1 * s2;
S = s1 + s2; E4 = p1 * p2;
E2 = sp + ps;
E3 = p1*s2 + p2*s1;

# Eq for (sp, ps):
E3^2 - sp*E3*S + E4*S^2 + ps*sp^2 - 4*ps*E4 # = 0
# =>
sp^3 - E2*sp^2 + sp*E3*S - 4*sp*E4 - E3^2 - E4*S^2 + 4*E2*E4 # = 0


### Special Cases:

### Case: b = 1
x1^3 + x2*x3*x4 - R # = 0

### Sub-Case: x2=x3=x4 = y;
x^3 + y^3 - R # = 0
y^3 + x*y^2 - R # = 0

### Diff =>
x*(x^2 - y^2) # = 0

### Solutions:
# Sol 1: x == y;
# - all 4 roots equal;

# Sol 2: x = 0;
# y^3 = R;

# Sol 3: y = -x; (NO)
# - NO solution (except when R = 0);


### Sub-Case: x3=x4 = y;
# - but distinct from x1 & x2;
x1^3 + x2*y^2 - R # = 0
x2^3 + x1*y^2 - R # = 0
y^3 + x1*x2*y - R # = 0

### Diff: Eq 1 - Eq 2 =>
(x1 - x2)*(x1^2 + x2^2 + x1*x2 - y^2) # = 0

### Sol 1: x1 == x2;
# - S2Ht system, same as: x1 = x3, x2 = x4;
x^3 + x*y^2 - R # = 0
y^3 + y*x^2 - R # = 0

### Sol 2: x1 != x2;
# - all are FALSE solutions!
# s1 = x1 + x2;
s1^2 - x1*x2 - y^2 # = 0
s1^2 - 3*x1*x2*s1 + s1*y^2 - 2*R # = 0
y^3 + x1*x2*y - R # = 0
# Sum: y*Eq 1 + Eq 3 =>
s1^2*y - R # = 0
# =>
s1^6 - x1*x2*s1^4 - R^2 # = 0
s1^3 + s1^2 - 4*x1*x2*s1 - 2*R # = 0
# =>
3*s1^6 - s1^5 + 2*R*s1^3 - 4*R^2 # = 0
# - but FALSE solutions;

