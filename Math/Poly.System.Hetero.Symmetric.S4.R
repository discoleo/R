########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogenous Symmetric
###
### draft v.0.1c



# "The many moods of an Irish setter"
#  and the many variants of S4!



### V1: x1^n + b*x2 = R
### V2: x1^n + b*x2*x3 = R
### V3: x1^n + b*x2*x3*x4 = R
### V4: x1^n + b*x1*x2*x3*x4 = R
### V_: x1^n + b*x1*x2 = R
### V_: x1^n + b*x1*x2*x3 = R
### ...

### TODO:
# - some proper classification;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R
# - e.g. round0(), round0.p(),
#   solve.EnAll(), solveEn();


########################

###############
### Order 2 ###
###############

### V3:
### x1^2 + b*x2*x3*x4 = R

### Sum =>
S^2 - 2*E2 + b*E3 - 4*R # = 0

### Sum(x1*...) =>
(x1^3 + x2^3 + x3^3 + x4^3) + 4*b*E4 - R*S # = 0
S^3 - 3*E2*S + 3*E3 + 4*b*E4 - R*S # = 0

### Sum(x1^2*...) =>
(x1^4 + x2^4 + x3^4 + x4^4) + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2 - 4*E4 + b*E4*S - R*(S^2 - 2*E2) # = 0
S^4 - R*S^2 + 2*E2^2 - 4*E2*S^2 + 2*R*E2 + 4*E3*S + b*E4*S - 4*E4 # = 0

### Sum(x2*x3*x4*...) =>
E4*S + b*Sum((x2*x3*x4)^2 ) - R*E3 # = 0
E4*S + b*(E3^2 - 2*E4*E2) - R*E3 # = 0
E4*S - 2*b*E2*E4 + b*E3^2 - R*E3 # = 0


### Debug:
R = 2;
b = 3;
#
x1 =  1.6128981492115229;
x2 = -0.5852711593612232;
x3 = x4 = x2;

x1^2 + b*x2*x3*x4 # - R
x2^2 + b*x1*x3*x4 # - R
x3^2 + b*x1*x2*x4 # - R
x4^2 + b*x1*x2*x3 # - R


### SEq1:
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2)*E3^1 +
(3*S*b^3 - 14*b^2)*E3^2

### SEq2:
(- 8*R*S^3 + 12*R*S^4*b - 8*R^2*S - 24*R^2*S^2*b - 32*R^3*b + S^5 - S^6*b) +
(8*R - 28*R*S*b + 10*R*S^2*b^2 + 32*R^2*b^2 - 8*S^2 + 10*S^3*b - 3*S^4*b^2)*E3^1 +
(- 10*R*b^3 + 9*S*b^2 - S^2*b^3 - 8*b)*E3^2 +
(b^4)*E3^3


### E3 =>
Subst = 560*R*S^2 - 780*R*S^3*b - 224*R*S^4*b^2 + 213*R*S^5*b^3 - 9*R*S^6*b^4 + 3024*R^2*S*b - 564*R^2*S^2*b^2 - 768*R^2*S^3*b^3 + 24*R^2*S^4*b^4 + 816*R^3*S*b^3 - 16*R^3*S^2*b^4 + 3136*R^3*b^2 - 56*S^4 + 36*S^5*b + 37*S^6*b^2 - 18*S^7*b^3 + S^8*b^4;
Subst = - Subst;

E3Div = - 324*R*S*b^2 - 256*R*S^2*b^3 + 249*R*S^3*b^4 - 5*R*S^4*b^5 + 1008*R*b - 252*R^2*S*b^4 + 4*R^2*S^2*b^5 - 1408*R^2*b^3 + 336*S - 188*S^2*b - 240*S^3*b^2 + 199*S^4*b^3 - 48*S^5*b^4 + S^6*b^5;

### =>
(10*R*S^2 - 14*R*S^3*b + 40*R^2*S*b - S^4 + S^5*b) +
(- 22*R*S*b^2 + 32*R*b + 6*S - 9*S^2*b + 4*S^3*b^2) * (Subst/E3Div) +
(3*S*b^3 - 14*b^2) * (Subst/E3Div)^2

### Eq:
b^5*S^10 +
	+ 2*b^4*S^9 - (2*b^3 + 6*R*b^5)*S^8  + (19*b^2 - 110*R*b^4)*S^7 +
	+ b*(-34 + 16*R*b^2 + 9*R^2*b^4)*S^6 + (-56 - 432*R*b^2 + 523*R^2*b^4)*S^5 +
	+ b*(568*R + 2026*R^2*b^2 - 4*R^3*b^4)*S^4 + (1120*R + 2472*R^2*b^2 - 596*R^3*b^4)*S^3 +
	- (2624*b*R^2 + 6608*b^3*R^3)*S^2 + (-3584*R^2 - 13184*b^2*R^3 + 256*b^4*R^4)*S +
	+ -7168*b*R^3 + 256*b^3*R^4
# (b*S^3 + 4*S^2 - 64*R) * (b*S + 1) * P[6]
### b*S + 1: Solution to "distinct" system (but x3 == x4);
### P[6]: Solution to degenerate System
# P[6] = (b*S^2 - 2*S - 4*b*R) * P[4]
(- 28*R + R^2*b^2) +
(- 24*R*b)*S^1 +
(7 - 2*R*b^2)*S^2 +
(- b)*S^3 +
(b^2)*S^4


#############
### Solution:

solve.S4 = function(R, b, max.perm=0, tol=1E-3, debug=FALSE, old=FALSE) {
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
	### Case: still x1 == x2 == x3
	# - but with different formula & numerically unstable;
	coeff = c(b1^2, - b1, (7 - 2*R*b1^2), - 24*R*b1, (- 28*R + R^2*b1^2))

	# Numerical instability of roots!
	S = roots(coeff); S = round0(S);
	S = c(S, -1/b1); # add the remaining roots;
	# E3
	Subst = 560*R*S^2 - 780*R*S^3*b - 224*R*S^4*b^2 + 213*R*S^5*b^3 - 9*R*S^6*b^4 + 3024*R^2*S*b +
		- 564*R^2*S^2*b^2 - 768*R^2*S^3*b^3 + 24*R^2*S^4*b^4 + 816*R^3*S*b^3 - 16*R^3*S^2*b^4 +
		3136*R^3*b^2 - 56*S^4 + 36*S^5*b + 37*S^6*b^2 - 18*S^7*b^3 + S^8*b^4;
	Subst = - Subst;

	E3Div = - 324*R*S*b^2 - 256*R*S^2*b^3 + 249*R*S^3*b^4 - 5*R*S^4*b^5 + 1008*R*b - 252*R^2*S*b^4 +
		4*R^2*S^2*b^5 - 1408*R^2*b^3 + 336*S - 188*S^2*b - 240*S^3*b^2 + 199*S^4*b^3 - 48*S^5*b^4 + S^6*b^5;
	E3 = round0(Subst / E3Div);
	E2 = round0(S^2 + b*E3 - 4*R) / 2;
	E4 = - round0(S^3 - 3*E2*S + 3*E3 - R*S) / 4 / b

	x = sapply(seq(length(S)), function(id) roots(c(1, -S[id], E2[id], -E3[id], E4[id])))
	# debugging: true roots
	if(debug) {
		E = list(S=S, E2=E2, E3=E3, E4=E4)
		return(debug.old(R, b, x, E, tol=tol))
	}
	#
	sol = if(old) solve.old(as.vector(x), E=list(S=S, E2=E2, E3=E3, E4=E4))
		else solve.EnAll(x, max.perm=max.perm); # generates 5*24 = 120 roots!
	return(list(sol=sol, sol3=sol3, sol22=sol22))
}

### TODO:
# - may still contain numerically unstable roots!

R = 2
b = 3
sol = solve.S4(R=R, b=b, max.perm=1)
sol.sol = sol$sol; # sol$sol22;
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
# - degenerates to a S2 system;
x^2 + b*x^2*x4 # - R
x4^2 + b*x^3 # - R

b^3*x^7 - b^2*R*x^4 + x^4 - 2*R*x^2 + R^2 # = 0
(x^2 + b*x^3 - R) * (b^2*x^4 - b*x^3 + x^2 - R) #= 0

### Case x1 == x2
# degenerates to a S3 system:
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
	x3 = sapply(seq_along(x2), function(id) roots(c(1, -SS2[id], E2S2[id]))[2])
	x4 = SS2 - x3; sol=cbind(x1=x, x2=x2, x3=x3, x4=x4);
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

### Hetero-Systems

### s^2 + c^2 = R1
### s1*c1 + a12*s1*c2 + a13*c1*s2 = R2


# s1^2 + c1^2 = R1
# s2^2 + c2^2 = R1
# s1*c1 + a2*s1*c2 + a3*c1*s2 = R2
# s2*c2 + a2*s2*c1 + a3*c2*s1 = R2

### Special Case:
# s1 = s2; c1 = c2;
s^2 + c^2 - R1 # = 0
(a2 + a3 + 1)*s*c - R2 # = 0


### Case:
# s1 != s2, c1 != c2
# TODO


