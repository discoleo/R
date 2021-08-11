########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Binomial Expansions
###
### draft v.0.2b-vars


### Asymmetric Polynomial Systems: 2 Variables
### Binomial Expansions


### Base: Class 1 Polynomials


###############
### History ###
###############


### draft v.0.2b-ht - v.0.2b-vars:
# - Ht-variant for Class 1 Order 3;
# - some concrete & special cases; (v.0.2b-sp)
# - Entangled variants: Ht & Diff-type; (v.0.2b-vars)
### draft v.0.2a:
# - System derived from Class 2 polynomials;
### draft v.0.1f:
# - solver for the Order 5 / Simple system;
### draft v.0.1e:
# - Multiplicative entanglement;
### draft v.0.1d - v.0.1d-Eq2:
# - Base Order 3 variants:
#   1st equation & 2nd Eq; (v.0.1d - v.0.1d-Eq2)
# - analysis of the variants; (v.0.1d-var)
### draft v.0.1c:
# - exact solution to the Order 3 system;
### draft v.0.1b:
# - Simple Order 5: Cardano-type;
### draft v.0.1a:
# - Order 3;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

solve.Cardano = function(c, d, n=3) {
	m = unity(n=n, all=TRUE);
	xdet = rootn(d^2 - c^n, 2);
	p = d + xdet; p = rootn(p, n);
	q = d - xdet; q = rootn(q, n);
	if(n %% 2 == 0 && c < 0) q = -q;
	sol = p*m + q/m;
	return(sol);
}


##########################

##########################
### Polynomial Systems ###
##########################

#####################
### Base: Class 1 ###
#####################

# Binomial expansions of Class 1 polynomials;

###############
### Order 3 ###
###############

### Base:
# x = k^2 + b11*k
# y = k^2 + b21*k
# where k^3 = K;

### Simple System:
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - 2*K^2 - (b11^3 + b21^3)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0

### Entangled Variants:

### 1.) Ht-Variant:
x^3 + 3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b21*K*y - 7*K^2 + b21^3*K - (b11+b21)^3*K # = 0
y^3 + 3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b11*K*x - 7*K^2 + b11^3*K - (b11+b21)^3*K # = 0

### 2.) Diff-Variant
x^3 - y^3 - 3*b11*K*x + 3*b21*K*y - (b11^3 - b21^3)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
# variant: y => -y;

### 3.) Other Variants:
# - see section on Multiplicative variants;


### Derivation:

### P(x^3) + P(y^3)
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - 2*K^2 - (b11^3 + b21^3)*K # = 0

### (x+y)^3 - (x^3+y^3)
3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b11*K*x + 3*b21*K*y - 6*K^2 - 3*b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - 2*(b11+b21)*K*(x+y) + b11*K*x + b21*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0

### Note:
### Eq 1 + 3*Eq 2 => S = (x+y) =>
#   S^3 - 6*(b11+b21)*K*S + ... = 0;
# - enables solving for S = (x+y);
### Step 2:
# x + y = S;
# x*y*S - b11*K*x - b21*K*y = ... - S^3 / 3;


### Solver:
solve.DP3 = function(K, b, all=TRUE) {
	d.f = function(b) 4*K^2 + b^3*K / 2;
	# can also use the direct formulas;
	bs = b[1] + b[2];
	S = solve.Cardano(2*bs*K, d.f(bs), n=3);
	x = sapply(S, function(S) {
		roots(c(-S, S^2 + (b[1] - b[2])*K, -2*K^2 - b[1]*b[2]*bs*K - (bs + b[1])*K*S))
	})
	S = rep(S, each=2);
	y = S - x;
	sol = cbind(x=as.vector(x), y=as.vector(y));
	return(sol);
}
# simple variant of solver (non-robust);
solve.DP3.old = function(K, b, all=FALSE) {
	c.f = function(b) b*K;
	d.f = function(b) (K^2 + b^3*K) / 2;
	# can also use the direct formulas;
	x = solve.Cardano(c.f(b[1]), d.f(b[1]), n=3);
	y = solve.Cardano(c.f(b[2]), d.f(b[2]), n=3);
	if(all) {
		# is actually NOT correct
		sol = expand.grid(x, y);
	} else sol = cbind(x=x, y=y);
	return(sol);
}

### Examples:

### Ex 1:
b = c(1,-2)
K = 3
#
sol = solve.DP3(K, b);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
# concrete Example
round0(x^3 + y^3 - 3*K*x + 6*K*y - 2*K^2 + 7*K) # = 0
round0(x*y*(x+y) + 3*K*x - 2*K^2 - 2*K) # = 0

### Ht-Variant / concrete:
x^3 + 3*x*y*(x+y) + 6*K*(x+y) - 6*K*y - 7*K^2 - 7*K # = 0
y^3 + 3*x*y*(x+y) + 6*K*(x+y) + 3*K*x - 7*K^2 + 2*K # = 0


### Ex 2:
b = c(-2,3)
K = 3
#
sol = solve.DP3(K, b, all=T);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
### Test
# only 3 solutions are correct, but not necessarily the base-set;
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y # = 2*K^2 + (b11^3 + b21^3)*K
err = x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K
round0(err)


#############
### Variants:

### non-Correlated Variants

### Base
# A^3 - 3*c1*A - 2*d1 = 0
# B^3 - 3*c2*B - 2*d2 = 0
# x = A + B
# y = A - B

# Note:
# - variants seem to be reducible;
#   (at least for base-order 3)

### System
x^4 - y^4 - x*y*(x^2-y^2) - 3*(c1+3*c2)*(x^2-y^2) - 4*(d1+3*d2)*x + 4*(d1-3*d2)*y # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 3*(3*c1+c2)*(x^2-y^2) - 4*(3*d1+d2)*x + 4*(3*d1-d2)*y # = 0

### Transformed
x^4 - y^4 - 6*(c1+c2)*(x^2-y^2) - 8*(d1+d2)*x + 8*(d1-d2)*y # = 0
x*y*(x^2-y^2) - 3*(c1-c2)*(x^2-y^2) - 4*(d1-d2)*x + 4*(d1+d2)*y # = 0
# Case: (x,y) != 0 => (x != y) =>
x^3 + y^3 + 3*x*y*(x+y) - 12*c1*(x+y) - 16*d1 # = 0

### Transformed: Variant
x^4 - y^4 - 6*(c1+c2)*(x^2-y^2) - 8*(d1+d2)*x + 8*(d1-d2)*y # = 0
x^3 + y^3 + 3*x*y*(x+y) + x*y*(x^2-y^2) - 3*(c1-c2)*(x^2-y^2) - 4*(d1-d2+3*c1)*x + 4*(d1+d2-3*c1)*y - 16*d1 # = 0
# Sum(Eq 1 + 2*Eq 2) =>
(x - y + 2)*(x^3 + y^3 + 3*x*y*(x+y) - 12*c1*(x+y) - 16*d1) # = 0


### Solution:
### Trivial solution
# x = y = 0;
### Non-trivial solution:
# see above & Derivation;

### Derivation:
### Eq 1: x^3 + Y^3:
x^3 + y^3 - 2*A^3 - 6*A*B^2 # = 0
x^3 + y^3 - 2*(3*c1*A + 2*d1) - 6*A*B^2 # = 0
# [simple/redundant]
4*x^3 + 4*y^3 - 12*c1*(x+y) - 16*d1 - 3*(x+y)*(x-y)^2 # = 0
4*x^3 + 4*y^3 - 12*c1*(x+y) - 16*d1 - 3*(x+y)*(x^2 + y^2 - 2*x*y) # = 0
x^3 + y^3 + 3*x*y*(x + y) - 12*c1*(x+y) - 16*d1 # = 0
# variant:
B*(x^3 + y^3) - 2*B*(3*c1*A + 2*d1) - 6*A*B^3 # = 0
(x-y)*(x^3 + y^3) - (x-y)*(3*c1*(x+y) + 4*d1) - 3*(x+y)*(3*c2*(x-y) + 4*d2) # = 0
x^4 - y^4 - x*y*(x^2-y^2) - 3*c1*(x^2-y^2) - 4*d1*(x-y) - 9*c2*(x^2-y^2) - 12*d2*(x+y) # = 0
x^4 - y^4 - x*y*(x^2-y^2) - 3*(c1+3*c2)*(x^2-y^2) - 4*(d1+3*d2)*x + 4*(d1-3*d2)*y # = 0

### Eq 2: x^3 - y^3
x^3 - y^3 - 2*B^3 - 6*A^2*B # = 0
x^3 - y^3 - 2*(3*c2*B + 2*d2) - 6*A^2*B # = 0
A*(x^3 - y^3) - 2*A*(3*c2*B + 2*d2) - 6*A^3*B # = 0
(x+y)*(x^3 - y^3) - 2*(x+y)*(3*c2*B + 2*d2) - 6*(3*c1*A + 2*d1)*(x-y) # = 0
(x+y)*(x^3 - y^3) - 2*(x+y)*(3*c2*B + 2*d2) - 3*(3*c1*(x+y) + 4*d1)*(x-y) # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 9*c1*(x^2-y^2) - (x+y)*(3*c2*(x-y) + 4*d2) - 12*d1*(x-y) # = 0
x^4 - y^4 + x*y*(x^2-y^2) - 3*(3*c1+c2)*(x^2-y^2) - 4*(3*d1+d2)*x + 4*(3*d1-d2)*y # = 0


### Examples
c = c(1, 2)
d  = c(2, 1)
#
c1 = c[1]; c2 = c[2]; d1 = d[1]; d2 = d[2];
A = solve.Cardano(c[1], d[1], n=3)
B = solve.Cardano(c[2], d[2], n=3)
x = A + B; y = A - B;
### "all" roots
sol = expand.grid(A, B);
x = sol[,1] + sol[,2]; y = sol[,1] - sol[,2];


### Test


##################

### Multiplicative Variants

### Base:
# x = k^2 + b11*k
# y = k^2 + b21*k
# where k^3 = K;

### System:
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - (2*K^2 + (b11^3 + b21^3)*K) # = 0
(b11+b21)*(x*y)^2 - K*((b11+b21)^2 + 2*b11*b21)*x*y +
	- b11*K*(K + b21^3)*x - b21*K*(K + b11^3)*y # = 0

### Special Case:
# b21 = - b11
x^3 + y^3 - 3*b11*K*x + 3*b11*K*y - 2*K^2 # = 0
2*b11*x*y - (K - b11^3)*x + (K + b11^3)*y # = 0


### Derivation:

### Prod =>
# (x*y) = b11*b21*k^2 + K*k + (b11+b21)*K;
# =>
(x*y)^3 - 3*(b11+b21)*K*(x*y)^2 + 3*(b11+b21)^2*K^2*x*y - 3*b11*b21*K^2*x*y +
	- K^4 + 3*b11*b21*(b11+b21)*K^3 - (b11+b21)^3*K^3 - (b11*b21)^3*K^2 # = 0
(3*b11*K*x + K^2 + b11^3*K)*(3*b21*K*y + K^2 + b21^3*K) - 3*(b11+b21)*K*(x*y)^2 + 3*(b11+b21)^2*K^2*x*y - 3*b11*b21*K^2*x*y +
	- K^4 + 3*b11*b21*(b11+b21)*K^3 - (b11+b21)^3*K^3 - (b11*b21)^3*K^2
K*(3*b11*x + K + b11^3)*(3*b21*y + K + b21^3) - 3*(b11+b21)*(x*y)^2 + 3*(b11+b21)^2*K*x*y - 3*b11*b21*K*x*y +
	- K^3 + 3*b11*b21*(b11+b21)*K^2 - (b11+b21)^3*K^2 - (b11*b21)^3*K
K*(9*b11*b21*x*y + 3*b11*K*x + 3*b11*b21^3*x + 3*b21*K*y + 3*b21*b11^3*y + K^2 + (b11^3+b21^3)*K + b11^3*b21^3) +
	- 3*(b11+b21)*(x*y)^2 + 3*((b11+b21)^2 - b11*b21)*K*x*y +
	- K^3 + 3*b11*b21*(b11+b21)*K^2 - (b11+b21)^3*K^2 - (b11*b21)^3*K
(b11+b21)*(x*y)^2 - ((b11+b21)^2 + 2*b11*b21)*K*x*y +
	- b11*K*(K + b21^3)*x - b21*K*(K + b11^3)*y # = 0


### Examples:

### Ex 1:
b = c(1,-2)
K = 3
#
b11 = b[1]; b21 = b[2];
k = rootn(K, 3);
x = k^2 + b[1]*k;
y = k^2 + b[2]*k;


### Ex 2: Special case
b = c(2,-2)
K = 3
#
b11 = b[1]; b21 = b[2];
k = rootn(K, 3);
x = k^2 + b[1]*k;
y = k^2 + b[2]*k;


#####################
#####################

###############
### Order 5 ###
###############

### Simple case:
### Base:
# x = k^4 + b11*k
# y = k^4 + b21*k
# where k^5 = K;

### System:
# x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y = 2*K^4 + (b11^5 + b21^5)*K
# Eq 2: see derivation;


### Derivation:

### P(x^3) + P(y^3)
x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y - 2*K^4 - (b11^5 + b21^5)*K # = 0

### (x+y)^5 - (x^5+y^5)
5*x*y*(x^3 + y^3) + 10*x^2*y^2*(x+y) - 10*(b11+b21)*K*(x+y)^3 + 20*(b11+b21)^2*K^2*(x+y) +
	+  5*b11*K*x^3 + 5*b21*K*y^3 - 5*b11^2*K^2*x - 5*b21^2*K^2*y - 30*K^4 - (b11+b21)^5*K + (b11^5 + b21^5)*K # = 0
5*x*y*(x^3 + y^3) + 10*x^2*y^2*(x+y) - 10*(b11+b21)*K*(x^3+y^3) - 30*(b11+b21)*K*x*y*(x+y) + 20*(b11+b21)^2*K^2*(x+y) +
	+  5*b11*K*x^3 + 5*b21*K*y^3 - 5*b11^2*K^2*x - 5*b21^2*K^2*y - 30*K^4 - (b11+b21)^5*K + (b11^5 + b21^5)*K # = 0
x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*(x+y) +
	+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
	- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K # = 0

### Solution:
x*y*(S^3 - 3*x*y*S) + 2*x^2*y^2*S - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*S +
	+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
	- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K # = 0
S*x^4 - (2*S^2 + (b11-b21)*K)*x^3 + (2*S^3 - 3*b21*K*S)*x^2 +
	- S^4*x + (b11^2 - b21^2)*K^2*x + 3*b21*K*S^2*x +
	- (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*S +
	+ 6*K^4 + b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K + (2*b11+b21)*K*S^3 # = 0


### Solver:
solve.DP5 = function(K, b) {
	bs = b[1] + b[2]; bd = b[1] - b[2]; bp = b[1]*b[2];
	S = solve.Cardano(c=2*bs*K, d=(16*K^4 + bs^5*K/2), n=5);
	x = sapply(S, function(S) {
		coeff = c(S, - (2*S^2 + bd*K), (2*S^3 - 3*b[2]*K*S),
			- S^4 + bs*bd*K^2 + 3*b[2]*K*S^2,
			- (b[1]^2 + 3*bs^2 + 2*bp)*K^2*S +
				+ 6*K^4 + bp*(bs^2 - bp)*bs*K + (b[1]+bs)*K*S^3);
			return(roots(coeff));
	})
	S = rep(S, each=4);
	y = S - x;
	sol = cbind(x = as.vector(x), y = as.vector(y));
	return(sol);
}
test.DP5 = function(x, y, K, b) {
	b11 = b[1]; b21 = b[2];
	err1 = x^5 + y^5 - 5*b11*K*x^3 - 5*b21*K*y^3 + 5*b11^2*K^2*x + 5*b21^2*K^2*y - 2*K^4 - (b11^5 + b21^5)*K;
	err2 = x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - (b11+2*b21)*K*x^3 - (2*b11+b21)*K*y^3 - 6*(b11+b21)*K*x*y*(x+y) +
		+ (3*b11^2 + 4*b21^2 + 8*b11*b21)*K^2*x + (4*b11^2 + 3*b21^2 + 8*b11*b21)*K^2*y +
		- 6*K^4 - b11*b21*((b11+b21)^2 - b11*b21)*(b11+b21)*K;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples
b = c(-1, 2)
K = 3
#
n = 5
b11 = b[1]; b21 = b[2];
k = rootn(K, n);
x = k^4 + b[1]*k;
y = k^4 + b[2]*k;
S = x + y;
### concrete Example:
x^5 + y^5 + 5*K*x^3 - 10*K*y^3 + 5*K^2*x + 20*K^2*y - 2*K^4 - 31*K # = 0
x*y*(x^3 + y^3) + 2*x^2*y^2*(x+y) - 3*K*x^3 - 6*K*x*y*(x+y) + 3*K^2*x - 6*K^4 + 6*K # = 0


### Ex 1:
b = c(-1, 2)
K = 3
sol = solve.DP5(K, b);
x = sol[,1]; y = sol[,2];

### Test
test.DP5(x, y, K, b)


### Classic Polynomial
p1 = toPoly.pm("x^5 + y^5 + 5*K*x^3 - 10*K*y^3 + 5*K^2*x + 20*K^2*y - 2*K^4 - 31*K");
p2 = toPoly.pm("x^4*y + x*y^4 + 2*x^3*y^2 + 2*x^2*y^3 - 3*K*x^3 - 6*K*x^2*y - 6*K*x*y^2 + 3*K^2*x - 6*K^4 + 6*K");
# Big-Integers
library(gmp)
p1$coeff = as.bigz(p1$coeff);
p2$coeff = as.bigz(p2$coeff);
pR = solve.pm(p1, p2, "y", asBigNum=TRUE);
str(pR) # pR$Rez = polynomial with 221 terms;
max(pR$Rez$x) # the question to the universe!
# Q: Are there also other roots?
# Or: Are the remaining False-roots?
# Solver: only 20 roots;


round0.p(poly.calc(x) * (32*K^3 + 1))


#####################
#####################

#####################
### Base: Class 2 ###
#####################

# Binomial expansions of Class 2 polynomials;

###############
### Order 4 ###
###############

### Base:
# x = s2*m^2 + s1*m
# y = s3*m^3 - s1*m
# where m^5 = 1;

### Simple System:
# x^4 + y^4 = ...;
# (x + y)^4 = ...;

### Entangled Variant:
# x^4 + 4*x*y*(x^2+y^2) + 6*x^2*y^2 = ...;
# y^4 + 4*x*y*(x^2+y^2) + 6*x^2*y^2 = ...;


### Solution:

### Simple system:
# Step 1: compute S = x + y;
# Step 2:
# - solve:
#   x + y = S;
#   the 1st Eq;


### Examples:

m = unity(5, all=FALSE);

s = c(3, 2, -1)
s1 = s[1]; s2 = s[2]; s3 = s[3];
x = s2*m^2 + s1*m;
y = s3*m^3 - s1*m;

### Test
### Eq 1:
x^4 + (s1 + s2)*x^3 + (s1^2 + 2*s1*s2 + s2^2)*x^2 + (s1^3 + 3*s1^2*s2 - 2*s1*s2^2 + s2^3)*x +
	+ s1^4 - s1^3*s2 + s1^2*s2^2 - s1*s2^3 + s2^4 +
y^4 - (s1 - s3)*y^3 + (s1^2 - 2*s1*s3 + s3^2)*y^2 - (s1^3 + 2*s1^2*s3 + 3*s1*s3^2 - s3^3)*y +
	+ s1^4 + s1^3*s3 + s1^2*s3^2 + s1*s3^3 + s3^4 # = 0
### Eq 2:
x^4 + y^4 + 4*x^3*y + 4*x*y^3 + 6*x^2*y^2 + 3*(s2+s3)*x*y*(x+y) + 2*(s2^2 - 3*s2*s3 + s3^2)*y*x +
	+ (s2 + s3) * (x^3 + y^3) + (s2^2 - 3*s2*s3 + s3^2) * (x^2 + y^2) +
	+ (s2^3 - 2*s2^2*s3 - 2*s2*s3^2 + s3^3) * (x + y) +
	+ s2^4 + s3^4 - s2^3*s3 + s2^2*s3^2 - s2*s3^3 # = 0


### Derivation:
p1 = toPoly.Class2.pm(4, s.id=c(1,2));
print.p(p1, "x")

p2 = toPoly.Class2.pm(4, s.id=c(1,3), xn="y");
p2 = replace.pm(p2, toPoly.pm("-s1"), "s1", 1);
print.p(p2, "y")

p3 = toPoly.Class2.pm(4, s.id=c(2,3));
p3 = replace.pm(p3, toPoly.pm("x+y"), "x", 1);
print.p(p3, "x")

