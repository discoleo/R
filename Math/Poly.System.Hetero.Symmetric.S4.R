########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
###
### draft v.0.3b



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
	E3 = E3.helper(S, R, b);
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
E3.helper = function(S, R, b) {
	pE3 = ((87*b^16 + 57*R*b^14 + 2987*R^2*b^12 - 14958*R^3*b^10 + 16796*R^4*b^8 + 1052*R^5*b^6 +
			- 11312*R^6*b^4 + 6464*R^7*b^2 - 1280*R^8)*S^2 +
		(727*R*b^15 - 9711*R^2*b^13 + 11120*R^3*b^11 + 30776*R^4*b^9 - 68788*R^5*b^7 + 49168*R^6*b^5 +
			- 13760*R^7*b^3 + 768*R^8*b + 111*b^17)*S +
		(- 948*R*b^16 + 6132*R^2*b^14 + 6224*R^3*b^12 - 35360*R^4*b^10 + 34736*R^5*b^8 - 11904*R^6*b^6 +
			+ 768*R^7*b^4 - 196*b^18));
	pDiv = ((- 28*b^13 - 262*R*b^11 + 2450*R^2*b^9 - 4616*R^3*b^7 + 2596*R^4*b^5 + 32*R^5*b^3 - 448*R^6*b)*S^2 +
		(114*b^14 + 116*R*b^12 - 12*R^2*b^10 - 7044*R^3*b^8 + 16236*R^4*b^6 - 14768*R^5*b^4 + 6208*R^6*b^2 +
			- 1280*R^7)*S +
		(80*R*b^13 - 2256*R^2*b^11 + 8432*R^3*b^9 - 10384*R^4*b^7 + 5248*R^5*b^5 - 1280*R^6*b^3 - 88*b^15));
	return(- pE3/pDiv);
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


#############################
#############################

#############################
### Mixed Terms: Type V2a ###
#############################

###############
### Order 2 ###
###############

### x1^2 + b*x1*x2 = R

### Solution:

### Case: all x[i] different;
# - for derivation, see file:
#   Poly.System.Hetero.Symmetric.S4.Derivation.R;

### Eq:
S^3 - R*(b^4 - 2*b^3 + 4*b^2 - 4*b + 4)*S

### Auxiliary eqs:
### Special Case: S = 0
E3 = 0; E2 = -2*R; E4 = R^2 / (b^2 + 1);

### Case: S != 0
E2 = -b*R*(b^2 - b + 2);
E3 = - (b+1) * R * S;
E4 = R*E3 / ((b+1)*S); # - R^2;


### Solver:
xip.f = function(x, R, b, n=2, p=1) {
	(R - x^n) / b[1] / x^p;
}
solve.S4P2.V2a = function(R, b, be=0, all=FALSE, debug=TRUE) {
	b0 = (b[1]^4 - 2*b[1]^3 + 4*b[1]^2 - 4*b[1] + 4);
	coeff = c(1, 0, - R*b0);
	noExt = TRUE; if(length(be) < 2) be = c(be, 0);
	if(any(be != 0)) {
		noExt = FALSE;
		coeff = coeff + c(be[2]*b0, be[1]*b0, 0); # only powers: 1 & 2;
	}
	S = roots(coeff);
	if(debug) print(c(S, 0));
	#
	solve0 = function(x1, R) {
		x2 = xip.f(x1, R, b, p=1);
		x3 = xip.f(x2, R, b, p=1);
		x4 = xip.f(x3, R, b, p=1);
		sol = cbind(x1, x2, x3, x4);
		return(sol);
	}
	### Special case: S == 0
	E3 = 0; E2 = -2*R; E4 = R^2 / (b^2 + 1);
	x1 = roots(c(1, 0, E2, -E3, E4));
	sol = solve0(x1, R);
	### Case: S != 0
	R1 = if(noExt) rep(R, length(S)) else R - sapply(S, function(S) sum(be * S^seq_along(be)));
	E2 = -b*R1*(b^2 - b + 2); # E2 = rep(E2, length(S));
	E3 = - (b+1) * R1 * S;
	E4 = - R1^2;
	x1 = sapply(seq_along(S), function(id) roots(c(1, -S[id], E2[id], -E3[id], E4[id])));
	R1 = rep(R1, each=4)
	sol2 = solve0(as.vector(x1), R=R1);
	sol = rbind(sol, sol2);
	### Case: all equal | pair-wise equal:
	# ((b+1)*x^2 - R) * ((b-1)*x^2 + R)
	if(all) {
		if(b != -1 || be[1] != 0) {
			x1 = roots(c(b+1 + 16*be[2], 4*be[1], -R));
			sol2 = cbind(x1, x1, x1, x1);
			sol = rbind(sol, sol2);
		}
		### Case: x1 = x3, x2 = x4;
		x1 = roots(c(b-1, 0, R));
		# x1 + x2 = 0
		x2 = - x1;
		sol2 = cbind(x1, x2, x1, x2);
		sol = rbind(sol, sol2);
	}
	return(sol);
}

### Examples:
R = 2
b = -3
sol = solve.S4P2.V2a(R, b);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
x1^2 + b*x1*x2 # - R
x2^2 + b*x2*x3 # - R
x3^2 + b*x3*x4 # - R
x4^2 + b*x4*x1 # - R


#########
### Ex 2:
R = -4
b = 3
be = 1
sol = solve.S4P2.V2a(R, b, be=be, all=TRUE);
x1 = sol[,1]; x2 = sol[,2]; x3 = sol[,3]; x4 = sol[,4];

### Test
S = x1+x2+x3+x4; ext = be[1]*S;
x1^2 + b*x1*x2 + ext # - R
x2^2 + b*x2*x3 + ext # - R
x3^2 + b*x3*x4 + ext # - R
x4^2 + b*x4*x1 + ext # - R


### Classic polynomial
# - see Derivation;


###########################
###########################

###############
### Order 3 ###
###############

### x1^3 + b*x1*x2 = R

### Solution:

### Case 1: all equal
# - 3 solutions;

### Case 2: x1=x3, x2=x4, distinct;
# x1^3 = x2^3 (6 solutions when distinct);

### Case 3: all x[i] distinct;
# - if (x1, x2, x3, x4) is a solution,
#   then the following are also solutions:
#   (m*x1, m^2*x2, m*x3, m^2*x4) & (m^2*x1, m*x2, m^2*x3, m*x4);
# - system is decomposable into:
#   P[4] o P[3] o P[6], where P[3] is degenerate;

### TODO


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

