########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S4
### Heterogeneous Symmetric
### Leading 1, NL Mixed: V2
###
### draft v.0.1a


### Type L1 V2
# - Leading: 1 var;
# - Non-Leading / Side-Chain: Mixed, 2 vars;

### Subtypes:
# V2a: x1^n + b*x1*x2 = R
# V2b: x1^n + b*x2*x3 = R
# V2c: x1^n + b*x1*x3 = R

### Note:
# - V2c: (x1, x3) and (x2, x4) are actually independent systems;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


test.S4Ht.V2b = function(sol, b, R=NULL, n=2) {
	x1 = sol[,1]; x2 = sol[,2];
	x3 = sol[,3]; x4 = sol[,4];
	err1 = x1^n + b*x2*x3;
	err2 = x2^n + b*x3*x4;
	err3 = x3^n + b*x4*x1;
	err4 = x4^n + b*x1*x2;
	err = rbind(err1, err2, err3, err4);
	if( ! is.null(R)) {
		err = err - R;
	}
	notNAN = apply(err, 2, function(x) ! any(is.nan(x)));
	err[ , notNAN] = round0(err[ , notNAN]);
	return(err);
}


##########################
##########################

##################
### Type: V2b  ###
##################

###############
### Order 2 ###
###############

### V3b:
### x1^2 + b*x2*x3 = R


### Solution:

### Eqs:

### Eq 1:
S^3 - 3*E2*S + (b + 3)*E3 - R*S # = 0
### Eq 2:
(b + 1)*E2*S - (b + 1)*(b^2 - b + 3)*E3 + (b^2 - b - 3)*R*S # = 0
### Eq 3:
4*(b^2 - 1)*E2^2 - (2*b^2 + b - 1)*E2*S^2 + (8*b^2 - 8*b - 16)*R*E2 +
	- (b^3 - 4*b - 3)*E3*S - (b^2 - 7*b - 7)*R*S^2 - 16*(b + 1)*R^2 # = 0
### Eq 4:
4*(b^4 + b^2)*E4 + (b^3 - 2*b + 3)*E3*S +
	- (6*b^2 - 8*b + 4)*E2^2 + (2*b^2 - 3*b + 1)*E2*S^2 - (12*b^2 - 24*b + 16)*R*E2 +
	+ (3*b^2 - 7*b + 7)*R*S^2 - (4*b^2 - 16*b + 16)*R^2 # = 0

### Eq S:
(b-1)*(b+1)*(b^2+1)*S^2 - (b^2 + 2*b + 2)*(3*b^2 - 2)*R # = 0


### Solver:

solve.S4Ht.V2bP2 = function(R, b, debug=TRUE, all=FALSE) {
	coeff = coeff.S4Ht.V2bP2(R, b=b);
	S = roots(coeff);
	S = rootn(S, n=2);
	S = c(S, -S);
	if(debug) print(S);
	#
	E2x0  = (b^3 + 2*b + 3)*S^2 + (2*b^2 - 8*b - 12)*R;
	E2div = (b + 1)*(3*b^2 - 4*b + 6);
	E2 = E2x0 / E2div;
	E3 = - (S^3 - R*S - 3*E2*S) / (b + 3);
	E4 = (b^3 - 2*b + 3)*E3*S +
		- (6*b^2 - 8*b + 4)*E2^2 + (2*b^2 - 3*b + 1)*E2*S^2 - (12*b^2 - 24*b + 16)*R*E2 +
		+ (3*b^2 - 7*b + 7)*R*S^2 - (4*b^2 - 16*b + 16)*R^2;
	E4 = - E4 / (4*b^2*(b^2 + 1));
	#
	len = length(S);
	x1 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], E2[id], -E3[id], E4[id]));
	})
	x1 = as.vector(x1);
	S  = rep(S, each=4);
	E4 = rep(E4, each=4);
	# Robust:
	e3  = E4 / x1;
	p23 = (R - x1^2) / b;
	x4  = e3 / p23;
	p12 = (R - x4^2) / b;
	x2 = p12 / x1;
	x3 = S - (x1 + x2 + x4);
	#
	sol = cbind(x1, x2, x3, x4);
	# Note: solution contains already all permutations!
	# if(all) sol = rbind(sol, sol[, c(2,3,4,1)]);
	return(sol);
}
coeff.S4Ht.V2bP2 = function(R, b) {
	coeff = c((b-1)*(b+1)*(b^2+1),
		- (b^2 + 2*b + 2)*(3*b^2 - 2)*R);
	return(coeff);
}
### Test:
test.S4Ht.V2bP2 = function(sol, b, R=NULL) {
	test.S4Ht.V2b(sol, b=b, R=R, n=2)
}

### Examples:

### Ex 1:
R = 1
b = -4
sol = solve.S4Ht.V2bP2(R, b)

test.S4Ht.V2bP2(sol, b=b)


### Ex 2:
R = 4
b = -5
sol = solve.S4Ht.V2bP2(R, b)

test.S4Ht.V2bP2(sol, b=b)


### Ex 3:
R = -3
b = 4
sol = solve.S4Ht.V2bP2(R, b)

test.S4Ht.V2bP2(sol, b=b)

