########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###
### draft v.0.2f-types



# "The many moods of an Irish setter"
#  and the many variants of S4!



### V1: x1^n + b*x2 = R
### V2a: x1^n + b*x1*x2 = R
### V2b: x1^n + b*x2*x3 = R
### V3a: x1^n + b*x1*x2*x3 = R
### V3b: x1^n + b*x2*x3*x4 = R
### V4: x1^n + b*x1*x2*x3*x4 = R
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
	# TDODO !!!
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


###########################
###########################
###########################

##############
### Simple ###
##############

### x[i]^n + b*x[i+1] = R

# x1^n + b*x2 = R
# x2^n + b*x3 = R
# x3^n + b*x4 = R
# x4^n + b*x1 = R

### Solution:

### Case 1: x1=x2=x3=x4
# P[n]: n solutions;

### Case 2: x1=x3, x2=x4, x1 != x3
# P[2] o P[(n^2 - n)/2];

### Case 3: all distinct;
# P[4] o P[(n^4 - n^2)/4];

# - System is decomposable into 3 subsystems:
#   P[n] * (P[2] o P[(n^2 - n)/2]) * (P[4] o P[(n^4 - n^2)/4]);


###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1] = R

### Solution:

### Case 1: x1=x2=x3=x4
# - 2 solutions;

### Case 2: x1=x3, x2=x4
# - 4 solutions: 2 overlap Case 1;

### Case 3: all distinct;
# Classic Poly: P[16 - 4] = P[12];
# S: P[3];

### Derivation:
# - see file: Poly.System.Hetero.Symmetric.S4.Derivation.R;

### Sum =>
S^2 - 2*E2 + b*S - 4*R # = 0
### Diff =>
(-E3^2 + E3*E2*S - E4*S^2) + b^6 # = 0
### Eq 3:
E4^2 - b^4*E4 + b*R^3*S - b^2*R^2*E2 + b^3*R*E3 - R^4 # = 0
### Eq 4:
E4^2 - b^4*E4 + 2*R^2*E4 + 2*R*E2*E4 - R*E3^2 - 2*R^2*E3*S +
	+ 2*R^3*E2 + R^2*E2^2 - R^3*S^2 + R^4 # = 0

### Eq:
S^3 + (3*b^2 - 4*R)*S - 4*b^3 # = 0

### Solver:
solve.Simple.S4P2 = function(R, b, debug=TRUE) {
	coeff = c(1, 0, 3*b[1]^2 - 4*R, - 4*b[1]^3);
	S = roots(coeff);
	if(debug) print(S);
	E2 = (S^2 + b*S - 4*R) / 2;
	E3 = E3.helper(S, E2, R, b);
	E4 = (-E3^2 + E3*E2*S + b^6) / S^2;
	#
	len = length(S);
	x1 = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id], E4[id])));
	x1 = as.vector(x1);
	x2 = xi.f(x1, R, b, n=2);
	x3 = xi.f(x2, R, b, n=2);
	x4 = xi.f(x3, R, b, n=2);
	sol = cbind(x1=x1, x2=x2, x3=x3, x4=x4);
	id = order(abs(x1), abs(Re(x1)));
	sol = sol[id,];
	return(sol);
}
xi.f = function(x, R, b, n=2) {
	(R - x^n) / b[1];
}
E3.helper = function(S, E2, R, b) {
	# see file:
	# Poly.System.Hetero.Symmetric.S4.Derivation.R;
	E3.helper.f(S, E2, R, b);
}
test.S4P2.Simple = function(sol, R, b, n=2) {
	x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + b*x2 # - R
	err2 = x2^n + b*x3 # - R
	err3 = x3^n + b*x4 # - R
	err4 = x4^n + b*x1 # - R
	err = rbind(err1, err2, err3, err4);
	if( ! missing(R)) err = err - R;
	err = round0(err);
	return(err);
}

### Examples:

R = -1
b = 3
sol = solve.Simple.S4P2(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

test.S4P2.Simple(sol, b=b)


### Test
x1^2 + b*x2 # - R
x2^2 + b*x3 # - R
x3^2 + b*x4 # - R
x4^2 + b*x1 # - R


###############

###############
### Order 3 ###
###############

### x[i]^3 + b*x[i+1] = R

### Solution:

### Case 1: x1=x2=x3=x4
# - 3 solutions;

### Case 2: x1=x3, x2=x4
# - 9 solutions: 3 overlap Case 1;

### Case 3: all distinct;
# Classic Poly: P[81 - 9] = P[72];
# S: P[18];

### TODO;


###########################
###########################

###########################
### Mixt Hetero-Systems ###
###########################

### x^2 + y^2 = R1
### x1*y1 + a2*x1*y2 + a3*y1*x2 = R2


# x1^2 + y1^2 = R1
# x2^2 + y2^2 = R1
# x1*y1 + a2*x1*y2 + a3*x2*y1 = R2
# x2*y2 + a2*x2*y1 + a3*x1*y2 = R2

### Solution:

### Special Case:
# x1 = x2; y1 = y2;
x^2 + y^2 - R1 # = 0
(a2 + a3 + 1)*x*y - R2 # = 0
# x1 = -x2; y1 = -y2;
x^2 + y^2 - R1 # = 0
(a2 + a3 - 1)*x*y + R2 # = 0

S^2 - 2*x*y - R1 # = 0

### Case:
# x1 != +/- x2; y1 != +/- y2;
# TODO:
# - are there any such solutions ???

### Solution
solve.S4P2 = function(R, a) {
	solve.part = function(div) {
		xy = R[2] / div;
		S = roots(c(1, 0, -2*xy - R[1]));
		xy = rep(xy, each=2);
		xy.diff = sqrt(R[1] - 2*xy + 0i)
		x = (S + xy.diff) / 2;
		y = (S - xy.diff) / 2;
		sol = cbind(x, y); sol = rbind(sol, sol[,2:1])
	}
	a.s = a[1] + a[2];
	sol2 = solve.part(a.s + 1);
	sol = cbind(sol2, sol2);
	# - R[2]
	R[2] = - R[2];
	sol2 = solve.part(a.s - 1);
	sol = rbind(sol, cbind(sol2, -sol2));
	return(sol)
}

### Examples
R = c(-1, 2)
a = c(3,4)
sol = solve.S4P2(R, a);
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^2 + y1^2 # - R[1]
x2^2 + y2^2 # - R[1]
x1*y1 + a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
x2*y2 + a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


######################

### Simple P3

### x^3 + y^3 = R1
### a1*x1*y2 + a2*y1*x2 = R2

# x1^3 + y1^3 = R1
# x2^3 + y2^3 = R1
# a1*x1*y2 + a2*x2*y1 = R2
# a1*x2*y1 + a2*x1*y2 = R2

### Solution:

### Diff: Eq 3 - 4
(a1 - a2)*(x1*y2 - x2*y1) # = 0
# assumption: a1 != a2
x1*y2 - x2*y1 # = 0
# x1*y2 = x2*y1 = R2 / (a1+a2);

### Mult Eq 1r * 2r
x1^3*y2^3 - y1^3*x2^3 + R1*(x2^3 + y1^3) - R1^2 # = 0
R1*(x2^3 + y1^3) - R1^2 # = 0
# x2^3 + y1^3 = R1;
# x1^3 + y2^3 = R1;
# =>
(x1+y2)^3 - 3*x1*y2*(x1+y2) - R1 # = 0
S12^3 - 3*R2/(a1+a2)*S12 - R1 # = 0
(a1+a2)*S12^3 - 3*R2*S12 - (a1+a2)*R1 # = 0

# =>
(x2+y1)^3 - 3*x2*y1*(x2+y1) - R1 # = 0
S21^3 - 3*R2/(a1+a2)*S21 - R1 # = 0
(a1+a2)*S21^3 - 3*R2*S21 - (a1+a2)*R1 # = 0

### Solver:
solve.simple.S4P3 = function(R, a, debug=TRUE) {
	as = a[1] + a[2];
	S12 = roots(c(as, 0, - 3*R[2], -as*R[1]))
	S21 = roots(c(as, 0, - 3*R[2], -as*R[1]))
	if(debug) print(S12);
	solve.S2 = function(S, xy) {
		xy.diff = sqrt(S^2 - 4*xy + 0i);
		x = (S + xy.diff) / 2;
		y = (S - x);
		sol = cbind(x, y)
		return(rbind(sol, sol[,2:1]))
	}
	xy = R[2] / as;
	sol1 = solve.S2(S12, xy);
	sol2 = solve.S2(S21, xy);
	sol = cbind(sol1, sol2)
	# (x1, y2, x2, y1)
	return(sol[,c(1,4,3,2)])
}

### Examples
R = c(-1, 2)
a = c(3, 4)
sol = solve.simple.S4P3(R, a);
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^3 + y1^3 # - R[1]
x2^3 + y2^3 # - R[1]
a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


######################

### P3

### x^3 + y^3 = R1
### x1*y1 + a1*x1*y2 + a2*y1*x2 = R2

# x1^3 + y1^3 = R1
# x2^3 + y2^3 = R1
# x1*y1 + a1*x1*y2 + a2*x2*y1 = R2
# x2*y2 + a1*x2*y1 + a2*x1*y2 = R2

### Solution:

### Solve Pseudo-Liniar:
# (a1^2 - a2^2) * x1*y2 = (a1 - a2)*R2 - (a1*x1*y1 - a2*x2*y2)
# (a1^2 - a2^2) * x2*y1 = (a1 - a2)*R2 + (a2*x1*y1 - a1*x2*y2)

### Mult Eq 1r * 2r
x1^3*y2^3 - y1^3*x2^3 + R1*(x2^3 + y1^3) - R1^2 # = 0
# =>
# x1^3 + y2^3 # =
2*R1 - (x2^3 + y1^3)
# =>
# x2^3 = ... - y1^3
# y1^3 = R1 - x1^3 *OR*
# y2^3 = R1 - x2^3
# => Mult =>
x2^3*y1^3 + x1^3*y1^3 - R1*(...) + R1*y1^3 + (...)*x1^3 # = 0

### Diff =>
# x1^3 - y2^3 = x2^3 - y1^3

### Mult =>
x1^3*x2^3 + y1^3*y2^3 + x1^3*y2^3 + x2^3*y1^3 - R1^2 # = 0
# also:
x1^3*x2^3 - y1^3*y2^3 + R1*(y1^3 + y2^3) - R1^2 # = 0

### TODO:
# - general case;

### Special Cases:
# x1 = x2; y1 != y2;
# valid if: (a2 - a1 + 1) == 0;
(a2 - a1 + 1)*y1 - (a2 - a1 + 1)*y2 # = 0
# y1 = y2; x1 != x2;
# valid if: (a1 - a2 + 1) == 0;
(a1 - a2 + 1)*x1 - (a1 - a2 + 1)*x2 # = 0


### Solver
solve.complete.S4P3 = function(R, a) {
	a.diff = a[1] - a[2];
	if(abs(a.diff) == 1) {
		# simple solution
		m.all = unity(3, all=TRUE)
		m = m.all[2]
		if(a.diff > 0) {
			# (a1*m^j + (a2+1))*y1 = R2 / x;
			id = 1; # TODO
			div = (a[1]*m^id + (a[2]+1));
			x3 = roots(c(1, -R[1], R[2]^3 / div^3));
			x = as.vector(sapply(rootn(x3, 3), function(x) x*m.all));
			y = R[2] / x / div;
			sol = cbind(x1=x, y1=y, x2=x, y2=y*m^id);
		} else {
			# (a2*m^j + (a1+1))*x1 = R2 / y;
			id = 1; # TODO
			div = (a[2]*m^id + (a[1]+1));
			y3 = roots(c(1, -R[1], R[2]^3 / div^3));
			y = as.vector(sapply(rootn(y3, 3), function(y) y*m.all));
			x = R[2] / y / div;
			sol = cbind(x1=x, y1=y, x2=x*m^id, y2=y);
		}
	} else if(round0(a[1] + a[2]) == 0) {
		# x1*y1 + x2*y2 = 2*R2;
	}
}

### Examples:

### Special case: x1 == x2
R = c(-1, 2)
a = c(4, 3)
sol = solve.complete.S4P3(R, a)
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];


### Special case: y1 == y2
R = c(-1, 2)
a = c(3, 4)
sol = solve.complete.S4P3(R, a)
x1 = sol[,1]; x2 = sol[,3]; y1 = sol[,2]; y2 = sol[,4];

### Test
x1^3 + y1^3 # - R[1]
x2^3 + y2^3 # - R[1]
x1*y1 + a[1]*x1*y2 + a[2]*x2*y1 # - R[2]
x2*y2 + a[1]*x2*y1 + a[2]*x1*y2 # - R[2]


### Debug
R = c(-1, 2)
a = c(3, 5)
x1 = -1.3908831512 - 2.5503258859i;
y1 =  1.4119576563 + 2.5831136271i;
x2 =  0.4925974849 + 0.7567717580i;
y2 =  0.2526662764 + 0.6178135417i;


R = c(-1, 2)
a = c(3, -3)
x1 = -0.1559428928 - 2.0642921485i;
y1 =  0.2316237980 + 2.0785818464i;
x2 = -0.7699227697 + 0.3433777972i;
y2 =  0.6635372661 - 0.7460946871i;


######################
######################

### 3D
# - possible extension to 3D sphere;

# x1^2 + y1^2 + z1^2 = R1
# x2^2 + y2^2 + z2^2 = R1
# x1*y2 + a1*y1*z2 + a2*z1*x2 = R2
# x2*y1 + a1*y2*z1 + a2*z2*x1 = R2
# x1*z2 + a3*y1*x2 + a4*z1*y2 = R3
# x2*z1 + a3*y2*x1 + a4*z2*y1 = R3

### Case:
# x1 = x2; y1 = y2; z1 = z2;
x^2 + y^2 + z^2 - R1 # = 0
x*y + a1*y*z + a2*x*z - R2 # = 0
x*z + a3*x*y + a4*y*z - R3 # = 0
### Case:
# x1 = -x2; y1 = -y2; z1 = -z2;
x^2 + y^2 + z^2 - R1 # = 0
x*y + a1*y*z + a2*x*z + R2 # = 0
x*z + a3*x*y + a4*y*z + R3 # = 0


### Variants Eqs [5]&[6]:
### V1:
# x1*y1*z1 = R3
# x2*y2*z2 = R3; # must be the same R3;
# - but does NOT preserve the negative solution;
### V2:
# x1*y1*z1 + a1*x1*y1*z2 + a2*y1*z1*x2 + a3*z1*x1*y2 = R3
# x2*y2*z2 + a1*x2*y2*z1 + a2*y2*z2*x1 + a3*z2*x2*y1 = R3
# - does NOT preserve the negative solution;
### V3:
# x1*y1*y2*z2 + a1*y1*z1*z2*x2 + a2*z1*x1*x2*y2 = R3
# x2*y2*y1*z1 + a1*y2*z2*z1*x1 + a2*z2*x2*x1*y1 = R3


#####################
#####################
#####################

# and now for the many variants
# and roots of the asymmetric types!

# x1^3*x2 + b1*x1*x2*x3*x4 = R1
# x2^3*x3 + b2*x1*x2*x3*x4 = R2
# x3^3*x4 + b3*x1*x2*x3*x4 = R3
# x4^3*x1 + b4*x1*x2*x3*x4 = R4

### Solution:

### =>
# x1^3*x2 = R1 - b1*x1*x2*x3*x4
### Prod =>
(x1*x2*x3*x4)^4 - b1*b2*b3*b4*(x1*x2*x3*x4)^4 +
	+ b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4)*(x1*x2*x3*x4)^3 +
	- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2)*
		(x1*x2*x3*x4)^2 +
	+ R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4)*(x1*x2*x3*x4) +
	- R1*R2*R3*R4 # = 0

### Special Case: b1*b2*b3*b4 = 1
b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4)*(x1*x2*x3*x4)^3 +
	- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2)*
		(x1*x2*x3*x4)^2 +
	+ R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4)*(x1*x2*x3*x4) +
	- R1*R2*R3*R4 # = 0

### Solver
solve.Pr2.S4P31 = function(R, b, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
	coeff = c(1- b1*b2*b3*b4,
		b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4),
		- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2),
		R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4),
		- R1*R2*R3*R4)
	p = roots(coeff);
	if(debug) print(p);
	len = length(p)
	Xij = sapply(p, function(p) R - b*p);
	X13sq = Xij[1,]*Xij[3,] / p;
	# X24sq = p^2 / X13sq;
	X13 = sqrt(X13sq + 0i); X13 = c(X13, -X13);
	p = c(p, p); len = length(p); Xij = cbind(Xij, Xij);
	X24 = p / X13;
	# TODO: 10 roots per p;
	x1 = rootn(Xij[1,]^3 / Xij[2,] * X13, 10);
	# NOT robust!
	# sapply(seq(len),
		# function(id) rootn(Xij[1,id]^27 * Xij[3,id]^3 / Xij[2,id]^9 / Xij[4,id], 80));
	x2 = Xij[1,] / x1^3; x3 = Xij[2,] / x2^3; x4 = Xij[3,] / x3^3;
	sol = cbind(x1=as.vector(x1), x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	invisible(sol);
}

### Examples:
R = c(1,2,3,4)
b = c(1,2,-2, -1/4)
#
sol = solve.Pr2.S4P31(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
# TODO: debug vs robust ???
x1^3*x2 + b[1]*x1*x2*x3*x4 # - R1
x2^3*x3 + b[2]*x1*x2*x3*x4 # - R2
x3^3*x4 + b[3]*x1*x2*x3*x4 # - R3
x4^3*x1 + b[4]*x1*x2*x3*x4 # - R4


###############################

### Variant

# x1^2*x2*x3 + b1*x1*x2*x3*x4 = R1
# x2^2*x3*x4 + b2*x1*x2*x3*x4 = R2
# x3^2*x4*x1 + b3*x1*x2*x3*x4 = R3
# x4^2*x1*x2 + b4*x1*x2*x3*x4 = R4

### Solution:
# Step 1: same as above;


### Solver
solve.Pr3.S4P211 = function(R, b, debug=TRUE) {
	R1 = R[1]; R2 = R[2]; R3 = R[3]; R4 = R[4];
	b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
	coeff = c(1- b1*b2*b3*b4,
		b1*b2*b3*b4*(R1/b1 + R2/b2 + R3/b3 + R4/b4),
		- (b1*b2*R3*R4 + b1*b3*R2*R4 + b1*b4*R2*R3 + b2*b3*R1*R4 + b2*b4*R1*R3 + b3*b4*R1*R2),
		R1*R2*R3*R4*(b1/R1 + b2/R2 + b3/R3 + b4/R4),
		- R1*R2*R3*R4)
	p = roots(coeff);
	if(debug) print(p);
	len = length(p)
	Xij = sapply(p, function(p) R - b*p);
	X13sq = Xij[1,]*Xij[3,] / p;
	# X24sq = p^2 / X13sq;
	X13 = sqrt(X13sq + 0i); X13 = c(X13, -X13);
	p = c(p, p); len = length(p); Xij = cbind(Xij, Xij);
	X24 = p / X13; X23 = Xij[2,] / X24; X12 = Xij[1,] / X13;
	x1 = sqrt(Xij[1,] / X23);
	x1 = c(x1, -x1);
	x2 = c(X12, X12) / x1; x3 = c(X13, X13) / x1; x4 = c(X24, X24) / x2;
	sol = cbind(x1=as.vector(x1), x2=as.vector(x2), x3=as.vector(x3), x4=as.vector(x4))
	invisible(sol);
}

### Examples:
R = c(1,2,3,4)
b = c(1,2,-2, -1/4)
#
sol = solve.Pr3.S4P211(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2*x2*x3 + b[1]*x1*x2*x3*x4 # - R1
x2^2*x3*x4 + b[2]*x1*x2*x3*x4 # - R2
x3^2*x4*x1 + b[3]*x1*x2*x3*x4 # - R3
x4^2*x1*x2 + b[4]*x1*x2*x3*x4 # - R4


### Ex 2:
R = c(-1,2,3,4)
b = c(1,2,-2, -3)
#
sol = solve.Pr3.S4P211(R, b)
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2*x2*x3 + b[1]*x1*x2*x3*x4 # - R1
x2^2*x3*x4 + b[2]*x1*x2*x3*x4 # - R2
x3^2*x4*x1 + b[3]*x1*x2*x3*x4 # - R3
x4^2*x1*x2 + b[4]*x1*x2*x3*x4 # - R4

# degenerate Polynomial
round0.p(poly.calc(x1)) * prod(b) * 11

