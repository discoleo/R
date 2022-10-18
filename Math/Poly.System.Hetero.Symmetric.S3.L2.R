########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S3
### Heterogeneous Symmetric
### Type: 2 & 3 Leading Terms
###
### draft v.0.1a-clean

### Leading Terms: 2 or more;
# - simple/univariate monomials;

# Note:
# - Special Cases (Diff) are in file:
#   Poly.System.Hetero.Symmetric.S3.Diff.R;
# - Mixed Leading Terms are in files:
#   Poly.System.Hetero.Symmetric.S2.Lnn.R;
#   Poly.System.Hetero.Symmetric.S3.Leading.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


#######################
#######################

#######################
### 2 Leading Terms ###
#######################

### Theory:

# a1*x^n + a2*y^n + P(x,y,z) = R
# a1*y^n + a2*z^n + P(y,z,x) = R
# a1*z^n + a2*x^n + P(z,x,y) = R

### Transform:
# - simplification of leading term;

### Eq 0:
# Sum =>
# (a1 + a2)*(x^n + y^n + z^n) + Ps - 3*R # = 0
# (x^n + y^n + z^n) + Ps/as - 3*R/as # = 0
# - where Ps = P(x,y,z) + P(y,z,x) + P(z,x,y);
#   as = a1 + a2;
# Note: Ps = symmetric polynomial;

### Transformed System:
# 1/a2 * Eq 1 + 1/a1 * Eq 3 - (Eq 0) =>
aa*x^n + P(x,y,z)/a2 + P(z,x,y)/a1 - Ps/as - ar*R # = 0
aa*y^n + P(y,z,x)/a2 + P(x,y,z)/a1 - Ps/as - ar*R # = 0
aa*z^n + P(z,x,y)/a2 + P(y,z,x)/a1 - Ps/as - ar*R # = 0
# - with:
#   aa = a1/a2 + a2/a1 - 1;
#   ar = (1/a1 + 1/a2 - 3/as) = (as^2 - 3*a1*a2) / (a1*a2*as);


###############

###############
### Order 2 ###
###############


### x[i]^2 + x[j]^2 + b*x[j]

# x^2 + y^2 + b1*y = R
# y^2 + z^2 + b1*z = R
# z^2 + x^2 + b1*x = R

### Solution

### Sum =>
2*(x^2 + y^2 + z^2) + b1*(x+y+z) - 3*R # = 0
2*S^2 - 4*E2 + b1*S - 3*R
# 4*E2 = 2*S^2 + b1*S - 3*R;

### Sum(z*...) =>
x^2*z+y^2*z + y^2*x+z^2*x + x^2*y+z^2*y + b1*E2 - R*S # = 0
E2*S - 3*E3 + b1*E2 - R*S # = 0

### Diff =>
# x^2 - z^2 = -b1*(y - z)
# Note: excludes x == y == z;
### Prod =>
(x+y)*(x+z)*(y+z) - b1^3 # = 0
x^2*z+y^2*z + y^2*x+z^2*x + x^2*y+z^2*y + 2*x*y*z - b1^3 # = 0
E2*S - E3 - b1^3 # = 0
# E3 = E2*S - b1^3

### =>
E2*S - 3*E3 + b1*E2 - R*S # = 0
8*E2*S - 4*b1*E2 + 4*R*S - 12*b1^3 # = 0
2*(2*S^2 + b1*S - 3*R)*S - b1*(2*S^2 + b1*S - 3*R) + 4*R*S - 12*b1^3 # = 0
4*S^3 - (2*R + b1^2)*S - 12*b1^3 + 3*b1*R # = 0

### Eq S:
(2*S - 3*b1)*(2*S^2 + 3*b1*S + 4*b1^2 - R)

### Alternatives:
### Redundant:
# Sum((x+y)*...), Sum(x*y*...);

### Alternative:
### Eq 3:
# Sum(y^2*...) =>
(x^2*y^2+y^2*z^2+z^2*x^2) + (x^4+y^4+z^4) + b1*(x^3+y^3+z^3) - R*(x^2+y^2+z^2) # = 0
3*E2^2 - 4*E2*S^2 + 2*E3*S +
	+ S^4 + b1*(S^3 - 3*E2*S + 3*E3) - R*(S^2 - 2*E2) # = 0


### Solver:

solve.S3Ht.L2sP2 = function(R, b, b.ext=0, debug=TRUE) {
	be1 = b.ext[1];
	be2 = if(length(b.ext) < 2) 0 else b.ext[2];
	# coeff = c(4, 0, - (2*R[1] + b[1]^2), -12*b[1]^3 + 3*b[1]*R[1])
	coeff = c(2, 3*b[1], 4*b[1]^2 - R[1])
	coeff = coeff + c(be2, be1, 0)
	S = round0(roots(coeff)) # numerical stability
	if(debug) print(S)
	R1 = R[1] - be1*S - be2*S^2;
	E2 = (2*S^2 + b[1]*S - 3*R1) / 4
	E3 = E2*S - b[1]^3;
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	x = as.vector(x);
	### fast prototype: was based on permutations;
	### robust
	S  = rep(S, each=3);
	E2 = rep(E2, each=3);
	R1 = rep(R1, each=3);
	#
	yz.s = S - x; yz = E2 - x*yz.s;
	y = 2*(R1 - x^2) - b[1]*x - yz.s^2 + 2*yz;
	y = y / b[1];
	z = yz.s - y;
	sol = cbind(x=as.vector(x), y=as.vector(y), z=as.vector(z))
	return(sol)
}
poly.S3Ht.L2sP2 = function(R, b, b.ext=0, max.leading=FALSE) {
	b1 = b[1]; b2 = b.ext[1];
	b3 = if(length(b.ext) > 1) b.ext[2] else 0;
	# ext: + b2*S + b3*S^2;
	coeff = c((b3 + 2)^3, (b3 + 2)^2*(3*b1 + b2),
		(b3 + 2)*(- 3*(b3 + 2)*R - 3*b1^2*(2*b3^2 + 2*b3 - 1) + 2*b1*b2*(b3 - 1) - b2^2),
		(- 2*(3*b1 + b2)*(b3 + 2)*R + b1^3*(2*b3^3 - 6*b3^2 - 6*b3 + 1) +
			- (6*b3^2 + 4*b3 + 11)*b1^2*b2 + b1*b2^2*(2*b3 - 5) - b2^3),
		(6*R^2 - 7*b1^2*R - 2*b1^4 + 6*b1*b2*R - 8*b1^3*b2 + b2^2*R + 2*b1^2*b2^2 + 3*b3*R^2 +
			8*b1^2*b3*R - 5*b1^4*b3 - 2*b1^3*b2*b3 + 3*b1^2*b2^2*b3 + 8*b1^2*b3^2*R + 7*b1^4*b3^2 +
			- 5*b1^3*b2*b3^2 + 9*b1^4*b3^3),
		(3*b1*R^2 - 2*b1^3*R - b1^5 + b2*R^2 + 8*b1^2*b2*R - 7*b1^4*b2 + 2*b1*b2^2*R + 5*b1^3*b2^2 + 3*b1^2*b2^3 +
			- 2*b1^3*b3*R + 5*b1^5*b3 + 2*b1^2*b2*b3*R - 6*b1^4*b2*b3 - 7*b1^3*b2^2*b3 - 2*b1^3*b3^2*R +
			5*b1^5*b3^2 + 11*b1^4*b2*b3^2 - 6*b1^5*b3^3),
		(-R^3 + 2*b1^2*R^2 - 2*b1*b2*R^2 - 2*b1^2*b3*R^2 +
			+ 3*b1^4*R + 4*b1^3*b2*R - 3*b1^2*b2^2*R + 5*b1^4*b3*R + 3*b1^3*b2*b3*R +
			+ b1^6*b3^3 + b1^6 - 12*b1^5*b2*b3 + b1^5*b2 + 7*b1^4*b2^2 - b1^3*b2^3 + 6*b1^6*b3 +
			+ 2*b1^4*b2^2*b3 - 5*b1^4*b3^2*R + 17*b1^6*b3^2 - 3*b1^5*b2*b3^2)
	);
	if(max.leading) coeff else rev(coeff);
}

### Examples:

R = -2
b = 4
sol = solve.S3Ht.L2sP2(R, b)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y # - R
y^2 + z^2 + b[1]*z # - R
z^2 + x^2 + b[1]*x # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.S3Ht.L2sP2(R, b) / 8;


### Extensions:

R = -2;
b = 1;
b.ext = c(1)
#
sol = solve.S3Ht.L2sP2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
x^2 + y^2 + b[1]*y + b.ext[1]*(x+y+z) # - R
y^2 + z^2 + b[1]*z + b.ext[1]*(x+y+z) # - R
z^2 + x^2 + b[1]*x + b.ext[1]*(x+y+z) # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.S3Ht.L2sP2(R, b, b.ext) / 8
err = 1 + 2*x^2 + 2*x^3 + 3*x^4 + 2*x^5 + x^6
round0(err)


### Ext 2:
R = 1;
b = 1;
b.ext = c(0, 1)
#
sol = solve.S3Ht.L2sP2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
S = x+y+z; ext = b.ext[1]*S + b.ext[2]*S^2;
x^2 + y^2 + b[1]*y + ext # - R
y^2 + z^2 + b[1]*z + ext # - R
z^2 + x^2 + b[1]*x + ext # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.S3Ht.L2sP2(R, b, b.ext) / 27
err = 1 + x^2 - x^3 - 2*x^4 + x^5 + x^6
round0(err)


### Ext 2 Ex 2:
R = 1;
b = 2;
b.ext = c(1, -1)
#
sol = solve.S3Ht.L2sP2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
S = x+y+z; ext = b.ext[1]*S + b.ext[2]*S^2;
x^2 + y^2 + b[1]*y + ext # - R
y^2 + z^2 + b[1]*z + ext # - R
z^2 + x^2 + b[1]*x + ext # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.S3Ht.L2sP2(R, b, b.ext)
err = 991 + 447*x - 88*x^2 - 89*x^3 + 7*x^5 + x^6
round0(err)


### Ext 2 Ex 3:
R = 6;
b = -1;
b.ext = c(3, -3)
#
sol = solve.S3Ht.L2sP2(R, b, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
S = x+y+z; ext = b.ext[1]*S + b.ext[2]*S^2;
x^2 + y^2 + b[1]*y + ext # - R
y^2 + z^2 + b[1]*z + ext # - R
z^2 + x^2 + b[1]*x + ext # - R

### Classic Polynomial
round0.p(poly.calc(x))
poly.S3Ht.L2sP2(R, b, b.ext)
err = 11 + 2*x + 5*x^2 - 2*x^3 + x^6
round0(err)


### Special Cases:
### b.ext[2] = -3
b1 = b[1]; b2 = b.ext[1]; b3 = b.ext[2];
(3*R*b1^2*b2^2 + 5*R*b1^3*b2 + 57*R*b1^4 + 2*R^2*b1*b2 - 8*R^2*b1^2 + R^3 + b1^3*b2^3 +
		- b1^4*b2^2 - 10*b1^5*b2 - 109*b1^6) +
	(- 2*R*b1*b2^2 - 2*R*b1^2*b2 + 14*R*b1^3 - 3*R^2*b1 - R^2*b2 - 3*b1^2*b2^3 - 26*b1^3*b2^2 +
		- 110*b1^4*b2 - 191*b1^5)*x +
	(- 6*R*b1*b2 - 41*R*b1^2 - R*b2^2 + 3*R^2 + 7*b1^2*b2^2 + 47*b1^3*b2 + 167*b1^4)*x^2 +
	(- 6*R*b1 - 2*R*b2 + 11*b1*b2^2 + 53*b1^2*b2 + 89*b1^3 + b2^3)*x^3 +
	(3*R - 8*b1*b2 - 33*b1^2 - b2^2)*x^4 +
	- (3*b1 + b2)*x^5 + x^6
### shifted
b1 = b[1] / 2; b2 = b.ext[1] / 6; x = x - b1 - b2;
(- 324*R*b1*b2^3 - 158*R*b1^2*b2^2 - 132*R*b1^3*b2 + 851*R*b1^4 - 45*R*b2^4 + 18*R^2*b1*b2 +
		- 35*R^2*b1^2 - 3*R^2*b2^2 + R^3 + 1170*b1*b2^5 + 1905*b1^2*b2^4 + 1692*b1^3*b2^3 +
		- 4975*b1^4*b2^2 - 8238*b1^5*b2 - 9841*b1^6 + 175*b2^6) +
	16*(- 27*R*b1*b2^2 - 37*R*b1^2*b2 - 15*R*b1^3 - 6*R*b2^3 + 171*b1*b2^4 + 393*b1^2*b2^3 +
		504*b1^3*b2^2 + 331*b1^4*b2 + 51*b1^5 + 30*b2^5)*x + # (2*b2 - 3*b1)*(...)
	(- 108*R*b1*b2 - 182*R*b1^2 - 54*R*b2^2 + 3*R^2 + 1836*b1*b2^3 + 4770*b1^2*b2^2 +
		5868*b1^3*b2 + 3971*b1^4 + 387*b2^4)*x^2 +
	16*(9*b1*b2^2 + 15*b1^2*b2 + 9*b1^3 + 2*b2^3)*x^3 + # (2*b2 - 3*b1)*(...)
	3*(R - 42*b1*b2 - 49*b1^2 - 17*b2^2)*x^4 + x^6


### Debug
R = 2
b = 3
x = -3.0643873807 + 0.5677216544i
y = -0.7500000000 + 2.3196254315i
z =  1.5643873807 + 0.5677216544i
S = x+y+z; E2 = x*(y+z)+y*z; E3 = x*y*z;


#########################

### Variant:

### x[i]^2 + x[j]^2 + b*(x[i] + x[j])

# x^2 + y^2 + b1*(x+y) = R
# y^2 + z^2 + b1*(y+z) = R
# x^2 + z^2 + b1*(x+z) = R

### Solution

### trivial solution: x = y = z;

### Diff =>
# x^2 - z^2 = -b1*(x - z)
# (x-z)*(x + z + b1) = 0

### Case: x != y != z
# - has NO solutions;

### Case: x = y && x != z
# - trivial;
# x^2 + b1*x = R/2
# z^2 + b1*z = R/2 # x & z are the conjugate roots;
# also:
# x + z = - b1
# x*z = - R / 2

### Example
b = 3
R = 1
#
x = roots(c(1, b[1], -R/2))
y = x
z = x[c(2,1)]
sol = cbind(x,y,z)
sol

### Test
x^2 + y^2 + b[1]*(x+y)
y^2 + z^2 + b[1]*(y+z)
x^2 + z^2 + b[1]*(x+z)

### TODO:
# - extensions;


########################
########################

### High-Power Terms: 2
### Variant

### x[i]^2 + x[j]^2 + b*x[k]

# x^2 + y^2 + b1*z = R
# y^2 + z^2 + b1*x = R
# x^2 + z^2 + b1*y = R

# - trivial solution: x = y = z;
# - trivial system;
# - equivalent to the previous variant: + b1*(x+y);

### Solution

### Diff =>
# x^2 - z^2 = b1*(x - z)
# (x-z)*(x + z - b1) = 0

### Case: x != y != z
# - has NO solutions;

### Case x = y && x != z
# z = -x + b1;
# =>
# 2*x^2 + b1*z - R = 0
# 2*x^2 - b1*x + b1^2 - R

### Example

b = 3
R = 1
#
x = roots(c(2, - b[1], b[1]^2 - R))
y = x
z = -x + b[1]
sol = cbind(x, y, z)
sol

### Test
x^2 + y^2 + b[1]*z
y^2 + z^2 + b[1]*x
x^2 + z^2 + b[1]*y


########################
########################
########################

########################
### Leading Terms: 3 ###
########################

###############
### Order 2 ###
###############

# x^2 + a1*y^2 + a2*z^2 = R
# y^2 + a1*z^2 + a2*x^2 = R
# z^2 + a1*x^2 + a2*y^2 = R

# - trivial type: without side chain;
#   [early experimental]

### Solution:
# - complicated solution;

### Sum =>
(a1 + a2 + 1)*(x^2 + y^2 + z^2) - 3*R # = 0
(a1 + a2 + 1)*(S^2 - 2*E2) - 3*R # = 0
# 2*(a1 + a2 + 1)*E2 = (a1 + a2 + 1)*S^2 - 3*R;

### Sum(x^2*...) =>
(x^4 + y^4 + z^4) + (a1 + a2)*((x*y)^2 + (x*z)^2 + (y*z)^2) - R*(x^2 + y^2 + z^2) # = 0
(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (a1 + a2)*(E2^2 - 2*E3*S) - R*(S^2 - 2*E2)
S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 - 2*(a1 + a2 - 2)*E3*S + 2*R*E2
# 2*(a1 + a2 - 2)*E3*S =
#   S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 + 2*R*E2

### Sum(z*...) =>
a2*(x^3 + y^3 + z^3) + (x^2*z + y^2*x + z^2*y) + a1*(y^2*z + z^2*x + x^2*y) - R*S # = 0
a2*(S^3 - 3*E2*S + 3*E3) + (E2*S - 3*E3) + (a1 - 1)*(y^2*z + z^2*x + x^2*y) - R*S # = 0
(a1 - 1)*(y^2*z + z^2*x + x^2*y) + a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S # = 0 # Eq 3a
# *(x^2*z + y^2*x + z^2*y) =>
(a1 - 1)*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*(x^2*z + y^2*x + z^2*y)
### Sum: Eq 3a + Eq 3b =>
(a1 - 1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*
	((a1 - 1)*(E2*S - 3*E3) + (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S))
(a1 - 1)^2*(E3*S^3 + E2^3 - 6*E3*E2*S + 9*E3^2) +
	+ (a2*S^3 - 3*a2*E2*S + E2*S + 3*a2*E3 - 3*E3 - R*S)*
	(a2*S^3 - 3*a2*E2*S + a1*E2*S + 3*a2*E3 - 3*a1*E3 - R*S)



### Eq S:
((a1 + a2 + 1)*S^2 - 9*R)^2 * ((a1 + a2 + 1)*S^2 - R)^2 # * P0;

### P[0]
(4 - 8*(a1 + a2) + 6*a1*a2 - 3*a1*a2*(a1 + a2) + a1*a2*(a1^2 + a2^2) +
	+ 9*(a1^2 + a2^2) - 5*(a1^3 + a2^3) + a1^4 + a2^4) 


### Q:
# - Do A1-type extensions have additional roots?
# - It seems NO additional roots possible!
### Technique: Sequential factorization
# - NO additional factors in this case;


### Solver:
solve.HP3.S3P2 = function(R, a, b.ext=0, debug=TRUE) {
	a.s = (a[1] + a[2] + 1);
	coeff = c(a.s, 0, -R) # only Non-equal roots!
	len = max(length(coeff), length(b.ext) + 1)
	coeff = c(rep(0, len - length(coeff)), coeff)
	b.all = c(rep(0, len - length(b.ext) - 1), rev(b.ext), 0)
	coeff = coeff + b.all;
	S = roots(coeff);
	if(debug) print(S);
	#
	pow = seq(length(b.ext));
	R1 = R[1] - sapply(S, function(S) sum(b.ext * (S^pow)));
	E2 = (a.s*S^2 - 3*R1) / a.s / 2;
	E3 = (S^4 - R1*S^2 - 4*E2*S^2 + (a.s + 1)*E2^2 + 2*R1*E2) /
		(2*(a.s - 3)*S) # TODO: a.s == 3
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	#
	sol = solve.EnAll(x, n = 3)
	return(sol)
}
test.HP3.S3P2 = function(sol, R, a, b.ext=0) {
	a = c(1, a) # TODO: include 1 in a;
	S = apply(sol, 1, sum)
	len = length(b.ext)
	ext = sapply(S, function(S) sum(b.ext * S^seq(len)));
	err = apply(sol, 1, function(sol) sum(a*sol^2))
	err = err + ext;
	return(round0(err))
}

### Examples:

R = -2
a = c(2, 3)
#
sol = solve.HP3.S3P2(R, a)
test.HP3.S3P2(sol, R, a)


### Ex 2:
R = -2
a = c(2, 3)
b.ext = c(-1)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 3:
R = -2
a = c(2, 3)
b.ext = c(-1, -1)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 4:
R = -5
a = c(2, 3)
b.ext = c(0, 10)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Ex 5:
R = -2
a = c(2, 3)
b.ext = c(-1, -1, 2)
#
sol = solve.HP3.S3P2(R, a, b.ext=b.ext)
test.HP3.S3P2(sol, R, a, b.ext=b.ext)


### Test
x^2 + a[2]*y^2 + a[3]*z^2 # - R
y^2 + a[2]*z^2 + a[3]*x^2 # - R
z^2 + a[2]*x^2 + a[3]*y^2 # - R

perm.gen = function(x) {
	len = length(x)
	id = seq(len)
	id.m = outer(id, id, function(i, j) ((i+j+1) %% len + 1))
	p.m = x[id.m]
	dim(p.m) = dim(id.m)
	p.m
}

R = 1;
a = c(1,2,3)
a1 = a[2]; a2 = a[3];
p.m = perm.gen(a)
d = det(p.m)

sol = solve(p.m, rep(R, 3))
x = sqrt(sol[1]); y = -x; z = -x;
S = x+y+z; E2 = x*y+x*z+y*z; E3 = x*y*z;


### Sum(y^2*...) =>
a1*(x^4 + y^4 + z^4) + (a2 + 1)*((x*y)^2 + (x*z)^2 + (y*z)^2) - R*(x^2 + y^2 + z^2) # = 0
### Diff: Eq2 - Eq3 =>
(a1 - 1)*((x*y)^2 + (x*z)^2 + (y*z)^2) - (a1 - 1)*(x^4 + y^4 + z^4) # = 0
### Case: a1 != 1
((x*y)^2 + (x*z)^2 + (y*z)^2) - (x^4 + y^4 + z^4) # = 0
(E2^2 - 2*E3*S) - (S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) # = 0
E2^2 + 6*E3*S + S^4 - 4*E2*S^2
# =>
(a1 + a2 - 2)*E2^2 + 3*(S^4 - R*S^2 - 4*E2*S^2 + (a1 + a2 + 2)*E2^2 + 2*R*E2) +
	+ (a1 + a2 - 2)*S^4 - 4*(a1 + a2 - 2)*E2*S^2
4*(a1 + a2 + 1)*E2^2 - 4*(a1 + a2 + 1)*E2*S^2 +
	+ (a1 + a2 + 1)*S^4 + 6*R*E2 - 3*R*S^2
- 2*(a1 + a2 + 1)*E2*S^2 + (a1 + a2 + 1)*S^4 - 3*R*S^2 # redundant


########################
########################

###############
### Order 3 ###
###############

# x^3 + a1*y^3 + a2*z^3 = R
# y^3 + a1*z^3 + a2*x^3 = R
# z^3 + a1*x^3 + a2*y^3 = R

# - trivial variant: without side chain;
#   [early experimental]

### Solution

### Sum =>
(a1 + a2 + 1)*(x^3 + y^3 + z^3) - 3*R # = 0
(a1 + a2 + 1)*(S^3 - 3*E2*S + 3*E3) - 3*R

### Sum(x^3*...) =>
(x^6 + y^6 + z^6) + (a1 + a2)*(x^3*y^3 + x^3*z^3 + y^3*z^3) - R*(x^3 + y^3 + z^3) # = 0
(- 2*E2^3 + 3*E3^2 - 12*E2*E3*S + 9*E2^2*S^2 + 6*E3*S^3 - 6*E2*S^4 + S^6) +
	(a1 + a2)*(E2^3 - 3*E3*E2*S + 3*E3^2) - R*(S^3 - 3*E2*S + 3*E3) # = 0
(a1 + a2 - 2)*E2^3 + 3*(a1 + a2 + 1)*E3^2 - 3*(a1 + a2 + 4)*E2*E3*S + 9*E2^2*S^2 - 6*E2*S^4 +
	+ 3*R*E2*S + 6*E3*S^3 - 3*R*E3 + S^6 - R*S^3 # = 0

### Sum(y*...) =>
a1*(x^4 + y^4 + z^4) + (x^3*y + y^3*z + z^3*x) + a2*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (x^3*y + y^3*z + z^3*x) + a2*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*(S^4 - 4*E2*S^2 + 4*E3*S + 2*E2^2) +
	+ (E2*(S^2 - 2*E2) - E3*S) + (a2 - 1)*(x*y^3 + y*z^3 + z*x^3) - R*S # = 0
a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S +
	+ (a2 - 1)*(x*y^3 + y*z^3 + z*x^3) # = 0 ### Eq 3a
### * (x^3*y + y^3*z + z^3*x) =>
(a2 - 1)*(E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4) +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S) *
	(x^3*y + y^3*z + z^3*x) ### Eq 3b
### Eq 3a + Eq 3b =>
(a2 - 1)^2*(E3*S^5 - 5*E2*E3*S^3 + 7*E3^2*S^2 + E2^2*E3*S + E2^4) +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S)^2 +
	(a1*S^4 - 4*a1*E2*S^2 + E2*S^2 + 4*a1*E3*S - E3*S + 2*a1*E2^2 - 2*E2^2 - R*S) *
	(a2 - 1)*(E2*S^2 - 2*E2^2 - E3*S)


### TODO:

### Eq:
S * ((a1 + a2 + 1)*S^3 - 27*R) * ((a1 + a2 + 1)^2 * S^6 + 27*R^2)

### Solver:
solve.3HT.S3P3 = function(R, a, b.ext=0, max.perm=1, debug=TRUE) {
	a.s = (a[1] + a[2] + 1);
	coeff = c(a.s^2, 0,0,0,0,0, 27*R^2) # only Non-equal roots!
	len = max(length(coeff), 2*length(b.ext) + 1)
	coeff = c(rep(0, len - length(coeff)), coeff)
	b.all = c(rep(0, len - length(b.ext) - 1), -2*27*R*rev(b.ext), 0)
	id2 = len - rev(2 * seq(1, length(b.ext)))
	b.all[id2] = b.all[id2] + 27*b.ext^2; # TODO: cross-products
	coeff = coeff + b.all;
	S = roots(coeff); S = c(S, 0)
	if(debug) print(S);
	#
	pow = seq(length(b.ext));
	R1 = R[1] - sapply(S, function(S) sum(b.ext * (S^pow)));
	E2 = round0(calc.E2(S, R1, a)) # !!! # TODO: has 6000 monoms & overflows !!!
	# E2[non-linearized] = P[E2^4, S^8];
	E3 = (3*R1 - a.s*(S^3 - 3*E2*S)) / 3 / a.s;
	#
	len = length(S)
	x = sapply(seq(len), function(id) roots(c(1, -S[id], E2[id], -E3[id])))
	#
	sol = solve.EnAll(x, n = 3, max.perm=max.perm)
	return(sol)
}
calc.E2 = function(S, R, a) {
	if(a[1] != -2 || a[2] != 2) stop("works only with a = c(-2, 2)!");
	# possible coefficients: BUT possible NOT exact!
	# TODO: find exact coefficients;
	E2Subst = 108*R*S^2 - 3240*R*S^6 + 11700*R*S^10 + 56610*R^2*S^7 - 112860*R^3*S^4 - 18*S^5 +
		+ 540*S^9 - 3000*S^13;
	E2Div = - 135*R + 4050*R*S^4 - 450*R*S^8 - 1404*R^2*S - 71055*R^2*S^5 + 10395*R^3*S^2 +
		+ 45*S^3 - 1350*S^7 + 7600*S^11;
	- E2Subst / E2Div;
}

### Examples:
R = -1
a = c(-2, 2)
b.ext = c(1)
sol = solve.3HT.S3P3(R, a, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; z = sol[,3];

### Test
ext = b.ext[1] * (x+y+z)
x^3 + a[1]*y^3 + a[2]*z^3 + ext # - R
y^3 + a[1]*z^3 + a[2]*x^3 + ext # - R
z^3 + a[1]*x^3 + a[2]*y^3 + ext # - R


### Debug
R = -1
a = c(-2, 2)
m = complex(re=cos(2*pi/3), im=sin(2*pi/3))
x = rootn(R / (a[1] + a[2] + 1), 3)
y = x*m; z = x*m^2;
a1 = a[1]; a2 = a[2];
S = x + y + z; E2 = x*(y+z) + y*z; E3 = x*y*z;

poly.calc(c(x+y+y, x+z+z, x+x+y, x+x+z, y+y+z, y+z+z)) * (a[1]+a[2]+1)^2

