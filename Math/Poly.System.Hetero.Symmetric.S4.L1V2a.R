########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
### Leading 1, NL Mixed
###
### draft v.0.1b


### Type L1 V2a
# - Leading: 1 var;
# - Non-Leading: Mixed / Multiple vars - V2a


### Subtype:
# V2a: x1^n + b*x1*x2 = R


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


####################

### History

# TODO:
# - keep separate or merge into:
#   Poly.System.Hetero.Symmetric.S4.L1V2.R ?

### v.0.1b
# - renamed file to:
#   Poly.System.Hetero.Symmetric.S4.L1V2a.R
### v.0.1a
# - moved specific V2a code from file:
#   Poly.System.Hetero.Symmetric.S4.R


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

