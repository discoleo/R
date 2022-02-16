########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Quasi-Asymmetric S2:
### Derived Polynomials
###
### draft v.0.1c-comment


####################
####################

### Helper Functions

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

source("Polynomials.Helper.R")


##############################

### Derived Polynomials
# - from Quasi-Asymmetric S2Ht;
# - see also Polynomials.Derived.P6.R
#   for another/similar approach (for P6);

###############
### Order 3 ###
###############

# - Exact solutions for the following P[9] polynomials;

### P[9]:
# b, k = parameters;
(x^3 + 4*b1*x + 8*b0)*(x^3 + b1*x - b0)^2 + 27*k^2*x^3 # = 0
x^9 + 6*b1*x^7 + 6*b0*x^6 + 9*b1^2*x^5 + 6*b0*b1*x^4 + (27*k^2 + 4*b1^3 - 15*b0^2)*x^3 - 12*b0^2*b1*x + 8*b0^3

### P[9] Variant:
# - only a scaling of b1 (when b2 == 0);
S^9 - 2*(k2^2 - 3*b1)*S^7 + 6*b0*S^6 + (k2^2 - 3*b1)^2*S^5 - 2*b0*(k2^2 - 3*b1)*S^4 +
	+ (27*k1^2 - (4*k2^3 - 18*k2*b1)*k1 - k2^2*b1^2 + 4*b1^3 - 15*b0^2)*S^3 +
	+ 4*b0^2*(k2^2 - 3*b1)*S + 8*b0^3


### Solver:
solve.P9.S2sp = function(k, b, sort=TRUE, debug=TRUE) {
	len = length(b);
	bL  = if(len == 2) c(1, 0) else if(len == 3) 1
		else if(len == 1) c(1,0,0) else stop("Only P3 => P9 supported!");
	b = c(bL, rev(b));
	if(length(k) < 2) k = c(k, 0);
	bE1 = b;
	bE1[4] = bE1[4] + k[1]; bE1[2] = bE1[2] - k[2];
	x = roots(bE1);
	#
	bE2 = b;
	bE2[4] = bE2[4] - k[1]; bE2[2] = bE2[2] + k[2];
	y = roots(bE2);
	if(debug) {print(x); print(y); }
	#
	s = expand.grid(x, y);
	S = apply(s, 1, sum);
	if(sort) S = sort.sol(matrix(S, ncol=1), ncol=1)[,1];
	sol = list(x=x, y=y, S=S);
	return(sol);
}
test.P9.S2sp = function(S, b, k) {
	if(length(b) != 2) warning("Not yet implemented!")
	b0 = b[1]; b1 = b[2];
	if(is.list(S)) S = S$S;
	if(length(k) == 1) {
		err = (S^3 + 4*b1*S + 8*b0)*(S^3 + b1*S - b0)^2 + 27*k^2*S^3;
	} else {
		k1 = k[1]; k2 = k[2];
		err = S^9 - 2*(k2^2 - 3*b1)*S^7 + 6*b0*S^6 + (k2^2 - 3*b1)^2*S^5 - 2*b0*(k2^2 - 3*b1)*S^4 +
			+ (27*k1^2 - (4*k2^3 - 18*k2*b1)*k1 - k2^2*b1^2 + 4*b1^3 - 15*b0^2)*S^3 +
			+ 4*b0^2*(k2^2 - 3*b1)*S + 8*b0^3
	}
	err = round0(err);
	return(err);
}

### Examples:
k = 3
b = c(2, -1)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))


#########
### Ex 2:
k = sqrt(2)*1i
b = c(2, -2)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))


#########
### Ex 3:
k = sqrt(1/27)*1i/2
b = c(1/2, -2)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))


#########
### Ex 4:
k = sqrt(1/3)*1i/2
b = c(1/2, -2)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))


#########
### Ex 5:
b = c(1/2, -1/3)
k = sqrt(-(4*b[2]^3 - 15*b[1]^2)/27)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))
x = S
1 + x - x^4 + x^5 + 3*x^6 - 2*x^7 + x^9


#########
### Ex 6: Variant
b = c(-2, 3)
k = c(2, 1)
sol = solve.P9.S2sp(k, b); S = sol$S;

### Test
test.P9.S2sp(S, b, k)

round0(poly.calc(S))
x = S
-64 - 128*x + 247*x^3 - 32*x^4 + 64*x^5 - 12*x^6 + 16*x^7 + x^9


### Derivation:

p1 = toPoly.pm("(S^3 - 3*xy*S + b1*S + 2*b0) * (S^2 - xy + b1) + k2*(2*k1 - k2*(S^2 - 2*xy))*S")
p2 = toPoly.pm("(S^2 - 4*xy)*(S^2 - xy + b1)^2 - (k2*(S^2 - 2*xy) - 2*k1)^2")
pDiv = toPoly.pm("(k2*S^2 + 2*k1 + 2*k2*b1)^2")
pRez = solve.pm(p1, p2, "xy");
p = pRez$Rez;
pRez = div.pm(p, pDiv, by="k1");
p = pRez$Rez;
p = sort.pm(p, "S")
top.pm(p, "S")

# P[13] = P[9] * P[4]


#####################
#####################

