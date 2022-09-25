########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
### Leading 1, NL Mixed: V3
###
### draft v.0.1c-ref


### Type L1 NLm V3
# - Leading: 1 var;
# - Non-Leading: Mixed, 3 vars;

### Subtypes:
# V3a: x1^n + b*x1*x2*x3 = R
# V3b: x1^n + b*x2*x3*x4 = R
# V3c: x1^n + b*x1*x2*x4 = R


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


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
#     (but works only if NO additional terms)
# - alternative: C2-Decomposition;

### Eq 1: Sum =>
# S[n] + b*E3 - 4*R = 0

### Eq 2: Sum(x1*...) =>
# S[n+1] + 4*b*E4 - R*S = 0

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
x1^2 + b*x1*x3^2 - R # = 0
x3^2 + b*x3*x1^2 - R # = 0

### Case 4:
# x1 = x2 and x3 != x4
x1^2 + b*x1*x3*x4 - R # = 0
x3^2 + b*x1^2*x4 - R # = 0
x4^2 + b*x1^2*x3 - R # = 0

### Case 5: x[1:4] all distinct;


### Eq 1: Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Eq 2: Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### [old]
### Diff Eq[i] - Eq[i+1] =>
# - non-robust for powers < ???;
(x1 - x2)*(x1 + x2 - b*x3*x4) # = 0
(x2 - x3)*(x2 + x3 - b*x1*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum =>
3*S - b*E2 # = 0

### =>
#   b*E2 = 3*S
# - b^2*E3 = b*S^2 - 6*S - 4*b*R
# - 4*b^3*E4 = b^2*S^3 - 12*b*S^2 - b^2*R*S + 18*S + 12*b*R

### Sum(x1^2*...) =>

### Eq S:
(b*S + 1) *
(b^2*S^4 - b*S^3 + 7*S^2 - 2*b^2*R*S^2 - 24*b*R*S + b^2*R^2 - 28*R) # = 0


###########
### Solver:

solve.S4Ht.V3bP2 = function(R, b, debug=TRUE, all=FALSE) {
	b1 = b[1];
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
	# - valid only for Case: x2 = x3 = x4!
	# - but there are NO distinct solutions anyway (for Power = 2);
	# - the (b*S + 1) does NOT work with this formula;
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

### Examples:

R = -5;
b = 3
sol.all = solve.S4Ht.V3bP2(R=R, b=b)
sol = sol.all$sol;
x1 = sol[,1]; x2 = sol[,2];
x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R

### Test: Case 3 equal
x1=x2=x3 = sol.all$sol3$sol[,1];
x4 = sol.all$sol3$sol[,2];
x1^2 + b*x2*x3*x4 # - R
x4^2 + b*x1*x2*x3 # - R

poly.calc(sol$sol3$S) * b^2


### Cases:

### Case x1 == x2 == x3, but != x4:
# - degenerates to a non-symmetric S2 system;
x^2 + b*x^2*x4 # - R
x4^2 + b*x^3 # - R

b^3*x^7 - b^2*R*x^4 + x^4 - 2*R*x^2 + R^2 # = 0
(x^2 + b*x^3 - R) * (b^2*x^4 - b*x^3 + x^2 - R) #= 0

### Case x1 == x2
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


############
### Workout:
debug.old = function(R, b, x, E, tol) {
	len = length(E$S)
	S1 = E$S; # 1 copy;
	S = matrix(E$S, ncol=len, nrow=4, byrow=T)
	E2 = matrix(E$E2, ncol=len, nrow=4, byrow=T)
	E3 = matrix(E$E3, ncol=len, nrow=4, byrow=T)
	E4 = matrix(E$E4, ncol=len, nrow=4, byrow=T)
	isZero = round0(E4/x - (R - x^2)/b[1], tol=tol) == 0
	isZ = apply(isZero, 2, all)
	return(list(sol=cbind(x=as.vector(x)), S=S1, isZ=isZ, isZero=isZero))
}
solve.old = function(x, E) {
	len = length(E$S);
	S = matrix(E$S, ncol=len, nrow=4, byrow=T)
	E2 = matrix(E$E2, ncol=len, nrow=4, byrow=T)
	E3 = matrix(E$E3, ncol=len, nrow=4, byrow=T)
	E4 = matrix(E$E4, ncol=len, nrow=4, byrow=T)
	SS3  = S - x; E2S3 = E2 - x*SS3; E3S3 = E4 / x;
	# x2
	x2 = sapply(seq_along(x), function(id) roots(c(1, -SS3[id], E2S3[id], -E3S3[id])))
	#
	x2 = as.vector(x2);
	x = rep(as.vector(x), each=3);
	SS2 = rep(as.vector(SS3), each=3) - x2;
	E2S2 = rep(as.vector(E3S3), each=3) / x2;
	# TODO: root[2]
	x3 = sapply(seq_along(x2), function(id) roots(c(1, -SS2[id], E2S2[id]))[1])
	x4 = rep(SS2, each=1) - x3; sol=cbind(x1=x, x2=x2, x3=x3, x4=x4);
}

b = -4:4
b = b[b != 0]
R = 2;

sapply(b, function(b) {
	sol = solve.S4(R=R, b, tol=5E-2); # tol=5E-2
	table(sol$isZ)[1]
	} )

R = 2
b = 3
sol = solve.S4(R=R, b=b, tol=5E-2)
poly.calc(sol$S[ ! sol$isZ]) * 9 *b^4

S = sol$S; # ...

### full Polynomial:
# coefficients are correct;
(- 16859136*R^5*b^4 + 602112*R^6*b^6) +
(- 18264064*R^4*b^3 - 23432192*R^5*b^5 + 344064*R^6*b^7)*S^1 +
(- 6322176*R^3*b^2 - 6547968*R^4*b^4 - 3177216*R^5*b^6 - 230400*R^6*b^8)*S^2 +
(- 702464*R^2*b + 6886656*R^3*b^3 + 20108480*R^4*b^5 + 3500096*R^5*b^7 + 27648*R^6*b^9)*S^3 +
(2728320*R^2*b^2 + 8715168*R^3*b^4 + 5966384*R^4*b^6 + 28224*R^5*b^8)*S^4 +
(219520*R*b - 1150912*R^2*b^3 - 6173432*R^3*b^5 - 2784304*R^4*b^7 - 76464*R^5*b^9)*S^5 +
(- 498624*R*b^2 - 2393144*R^2*b^4 - 1846004*R^3*b^6 - 274868*R^4*b^8 - 432*R^5*b^10)*S^6 +
(125832*R*b^3 + 1536372*R^2*b^5 + 163910*R^3*b^7 + 84660*R^4*b^9 - 10976*b)*S^7 +
(362644*R*b^4 - 294070*R^2*b^6 + 293519*R^3*b^8 + 1224*R^4*b^10 + 19992*b^2)*S^8 +
(- 296758*R*b^5 + 260163*R^2*b^7 - 37878*R^3*b^9 - 980*b^3)*S^9 +
(127503*R*b^6 - 92493*R^2*b^8 - 1251*R^3*b^10 - 16406*b^4)*S^10 +
(- 48741*R*b^7 + 6021*R^2*b^9 + 11395*b^5)*S^11 +
(8594*R*b^8 + 567*R^2*b^10 - 3746*b^6)*S^12 +
(147*R*b^9 + 645*b^7)*S^13 +
(- 117*R*b^10 + 151*b^8)*S^14 +
(- 84*b^9)*S^15 +
(9*b^10)*S^16
# can be factored into: P[10] * P[6]
# P[6]: 9*b^4 * S^6 + ...; # TODO: check: is NOT part of solution!
# P[10]: can be factored itself in P[3]*P[7];
# P[10] = P[3]*P[1]*P[2]*P[4];


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
E2*S^2 - 2*E2^2 + bd*E3*S - 4*bd*E4 - 3*R*S # = 0


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

