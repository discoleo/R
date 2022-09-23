########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
### Leading 1, NL Mixed: V3
###
### draft v.0.1a


### Type L1 NLm V3
# - Leading: 1 var;
# - Non-Leading: Mixed, 3 vars;

### Subtypes:
# V3a: x1^n + b*x1*x2*x3 = R
# V3b: x1^n + b*x2*x3*x4 = R


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

##################
### Type: V3a  ###
##################

###############
### Order 2 ###
###############

### V3a:
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

### Diff Eq[i] - Eq[i+1] =>
(x1 - x2)*(x1 + x2 - b*x3*x4) # = 0
(x2 - x3)*(x2 + x3 - b*x1*x4) # = 0
# ...
# Case: x[i] != x[j]: Sum =>
3*S - b*E2 # = 0

### Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### =>
#   b*E2 = 3*S
# - b^2*E3 = b*S^2 - 6*S - 4*b*R
# - 4*b^3*E4 = b^2*S^3 - 12*b*S^2 - b^2*R*S + 18*S + 12*b*R

### Sum(x1^2*...) =>
b^3*S^4 - 16*b^2*S^3 - b^3*R*S^2 + 34*b*S^2 + 24*b^2*R*S + 24*S + 16*b*R # = 0

### TODO:
# - may still be a false root!


###########
### Solver:

solve.S4 = function(R, b, max.perm=0, tol=1E-3, debug=TRUE, old=FALSE) {
	# tol = was used for debugging;
	b1 = b[1];
	### Case: x1 == x2 == x3, but != x4;
	coeff.3eq = c(b1^2, - b1, 1, 0, - R);
	x3 = roots(coeff.3eq);
	y = (R - x3^2) / b1 / x3^2;
	S.3eq = 3*x3 + y;
	sol3 = list(sol=cbind(x123=x3, x4=y), S=S.3eq);
	
	### Case: x1 == x2, x3 == x4, but x1 != x3;
	S2 = roots(c(b1, -1, - b1*R))
	xy2 = S2 / b1;
	x2 = sapply(seq_along(S2), function(id) roots(c(1, -S2[id], xy2[id])))
	y2 = x2[2:1, ];
	x2 = as.vector(x2); y2 = as.vector(y2);
	sol22 = cbind(x1=x2, x2=x2, x3=y2, x4=y2); # + many permutations;
	
	### Case: x1 = x2, x3 != x4; [only 1 real case]
	coeff = c(b1^2, b1, 1 - 2*b1^2*R, -2*b1*R, b1^2*R^2);
	S34 = roots(coeff);
	S34 = c(1/b1, S34); # the real case: (b*S - 1) * P[4]
	if(debug) print(S34);
	p = S34^2 - R;
	x34.d = sqrt(S34^2 - 4*p + 0i); # TODO: +/-;
	x3 = (S34 + x34.d) / 2; x4 = (S34 - x34.d) / 2;
	# robust
	x1 = ((b1*x4-1)*R + x3^2) / (b1^2*p*x4);
	# contains also the variants: x1 == x2 == x3 != x4;
	sol22 = rbind(sol22, cbind(x1=x1, x2=x1, x3=x3, x4=x4));
	
	### Case: x[i] != x[j]
	# TODO !!!
	# b^3*S^4 - 16*b^2*S^3 - b^3*R*S^2 + 34*b*S^2 + 24*b^2*R*S + 24*S + 16*b*R
	coeff = c(b1^3, - 16*b1^2, - b1^3*R + 34*b1, 24*b1^2*R + 24, 16*b1*R);
	S = roots(coeff);
	if(debug) print(S);
	E2 = 3*S / b1;
	E3 = - (S^2 - 2*E2 - 4*R) / b1;
	E4 = - (S^3 - 3*E2*S + 3*E3 - R*S) / (4*b1);

	x = sapply(seq(length(S)), function(id) roots(c(1, -S[id], E2[id], -E3[id], E4[id])));
	# debugging: true roots
	if(FALSE) {
		E = list(S=S, E2=E2, E3=E3, E4=E4)
		return(debug.old(R, b, x, E, tol=tol))
	}
	#
	sol = if(old) solve.old(as.vector(x), E=list(S=S, E2=E2, E3=E3, E4=E4))
		else solve.EnAll(x, max.perm=max.perm); # generates 5*24 = 120 roots!
	return(list(sol=sol, sol3=sol3, sol22=sol22))
}

### Examples:

R = -5;
b = 3
sol = solve.S4(R=R, b=b, max.perm=1)
sol.sol = sol$sol22; # sol$sol22; # sol$sol; #
x1 = sol.sol[,1]; x2 = sol.sol[,2];
x3 = sol.sol[,3]; x4 = sol.sol[,4];

### Test
x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R

### Test: Case 3 equal
x1=x2=x3=sol$sol3$sol[,1];
x4 = sol$sol3$sol[,2];
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

