########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S2:
### Mixed Type: Composite
###
### draft v.0.1e


### Heterogeneous Symmetric
### Polynomial Systems: 2 Variables
### Mixed/Composite: Hetero + Symmetric

### Example:
# x^p + b*y = Ru # where Ru = unknown
# y^p + b*x = Ru
# x*y = R2

### Variants:
# x^n + y^n = R2
# x*y*(x+y) = R2
# or
# another Hetero-Symmetric system;


######################
######################

###############
### History ###
###############


### draft v.0.1e:
# - experimental systems: difference of powers;
### draft v.0.1d - v.0.1d-poly:
# - Variant: sum of powers;
#   x^3 + y^3 = R;
# - parametric P[6];
# - example:
#   3 - 3*x + 3*x^2 - 2*x^3 + 3*x^4 + x^6 = 0;
### draft v.0.1b - v.0.1c:
# - Order 3 & Order 4 variants:
#   x*y*(x+y) = R; [v.0.1b]
#   x*y = R; [v.0.1c]
# - example P[6]:
#   1 - 2*x + x^2 + 4*x^4 + 3*x^6 = 0; [v.0.1b-ord4]
# - full parametric polynomial: P[6]; [v.0.1b-poly]
### draft v.0.1a:
# - solved: double Ht-Symmetric Composite;
# - example P[6]:
#   3 - 2*x + x^2 + 4*x^3 - x^4 + x^6 = 0;


######################
######################

################
### Order: 3 ###
################

# x^3 + b2*y^2 + b1*y = Ru1
# y^3 + b2*x^2 + b1*x = Ru1
# x*y*(x+y) = R

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### =>
S^3 - b2*S^2 - b1*S - R # = 0


### Solver:
solve.HtComposite.S2P3 = function(R, b, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	coeff = c(1, -b[2], -b[1], -R[1]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

R = -1
b = c(1, -2)
sol = solve.HtComposite.S2P3(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^3 + b[2]*y^2 + b[1]*y, y^3 + b[2]*x^2 + b[1]*x)
x*y*(x+y) # - R[1]

###
round0.p(poly.calc(x))
err = 1 + 6*x - x^2 - 4*x^3 - 2*x^4 + 2*x^5 + x^6
round0(err)


######################

################
### Order: 4 ###
################

# x^4 + b2*y^2 + b1*y = Ru1
# y^4 + b2*x^2 + b1*x = Ru1
# x*y*(x+y) = R

### Simple Variant:
# x*y = R;

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^3 - 2*x*y*S - b2*S - b1 # = 0
### =>
S^3 - b2*S - b1 - 2*R # = 0

### Simple Variant:
S^3 - (2*R + b2)*S - b1 # = 0

### Solver:
solve.HtComposite.S2P4 = function(R, b, type=c("PSum", "Simple"), k=2, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	type = pmatch(type[1], eval(formals(solve.HtComposite.S2P4)$type));
	if(type == 1) {
		coeff = c(1, 0, -b[2], -b[1] - k*R[1]);
	} else {
		coeff = c(1, 0, -b[2] - k*R[1], -b[1]);
	}
	S = roots(coeff);
	if(debug) print(S);
	S = S[S != 0]; # exclude 0!
	xy = if(type == 1) R[1] / S else R[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

R = -1
b = c(-1, -1)
sol = solve.HtComposite.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x*y*(x+y) # - R[1]

###
round0.p(poly.calc(x) * (2*R + b[1]))
err = 1 - 2*x + x^2 + 4*x^4 + 3*x^6
round0(err)


#########
### Ex 2:
R = -1
b = c(-1, 2)
sol = solve.HtComposite.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x*y*(x+y) # - R[1]

###
round0.p(poly.calc(x) * (2*R + b[1]))
err = 1 + 4*x + 4*x^2 - 8*x^4 + 3*x^6
round0(err)

### parametric
R^3 - 2*b[2]*R^2*x + R*b[2]^2*x^2 + (2*R^2 - b[1]*R - b[1]^2)*x^3 +
	- b[2]*(3*R + b[1])*x^4 + (2*R + b[1])*x^6
# - special case: b1 = 1 - 2*R
R^3 - 2*b[2]*R^2*x + R*b[2]^2*x^2 + (3*R - 1)*x^3 +
	- b[2]*(R + 1)*x^4 + x^6


#########
### Ex 3: Simple type
R = -1
b = c(1, -3)
sol = solve.HtComposite.S2P4(R, b, type="Simple")
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x*y # - R[1]

###
round0.p(poly.calc(x))
err = -1 - 2*x^2 - x^3 + 2*x^4 + x^6
round0(err)
# (x^3 - 1/2)^2 + 2*(x^2 - 1/2)^2 - 7/4


#########
### Ex 4:
R = -1
b = c(-1, -1)
k = 3
sol = solve.HtComposite.S2P4(R, b, k=k)
x = sol[,1]; y = sol[,2];

### Test
# [fails] k != 2

###
round0.p(poly.calc(x) * -(k*R + b[1]))
err = 1 - 2*x + x^2 + 4*x^3 + 5*x^4 + 4*x^6
round0(err)


######################
######################

#################
### Variant:  ###
### Power Sum ###
#################

# x^4 + b2*y^2 + b1*y = Ru1
# y^4 + b2*x^2 + b1*x = Ru1
# x^3 + y^3 = R

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^3 - 2*x*y*S - b2*S - b1 # = 0

### Eq 3:
S^3 - 3*x*y*S - R # = 0
# 3*x*y*S = S^3 - R;
# =>
3*S^3 - 2*(S^3 - R) - 3*b2*S - 3*b1 # = 0
S^3 - 3*b2*S - 3*b1 + 2*R # = 0

### Solver:
solve.HtComposite.S2P4 = function(R, b, k=2, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	coeff = c(1, 0, -3*b[2], -3*b[1] + k*R[1]);
	S = roots(coeff);
	if(debug) print(S);
	S = S[S != 0]; # exclude 0!
	xy = (S^3 - R) / (3*S);
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

R = 2
b = c(1, 3)
sol = solve.HtComposite.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x^3 + y^3 # - R[1]

###
round0.p(poly.calc(x) * (2*R - 3*b[1]))
err = -53 + 9*x + 27*x^2 - 2*x^3 - 9*x^4 + x^6
round0(err)

### parametric P[6]
(R - b[1])^3 - R*b[2]^3 + 3*b[1]*b[2]*(R - b[1])*x +
	+ 3*b[2]^2*(R - b[1])*x^2 - R*(2*R - 3*b[1])*x^3 +
	- 3*b[2]*(R - b[1])*x^4 + (2*R - 3*b[1])*x^6

#########
### Ex 2:
R = 2
b = c(1, -1)
sol = solve.HtComposite.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x^3 + y^3 # - R[1]

###
round0.p(poly.calc(x) * (2*R - 3*b[1]))
err = 3 - 3*x + 3*x^2 - 2*x^3 + 3*x^4 + x^6
round0(err)
# x^3*(x^3+1)*(x+1) + 3*(x^5 + 1)

#########
### Ex 3:
R = 2
b = c(-1, 1)
k = -R # k = -1 # also interesting
sol = solve.HtComposite.S2P4(R, b, k=k)
x = sol[,1]; y = sol[,2];

### Test
# [fails] k != 2

###
round0.p(poly.calc(x) * (k*R - 3*b[1]) * -27)
err = 53 - 45*x - 27*x^2 - 54*x^3 + 27*x^4 + 27*x^6
round0(err)


#######################
#######################

#######################
### Hidden Symmetry ###
#######################

# x^4 + b2*y^2 + b1*y = Ru1
# y^4 + b2*x^2 + b1*x = Ru1
# x^2 - y^2 = R

### Solution:

# Note:
# - we win an sqrt("R"), but (x-y) is directly computable from S;

### Diff Eqs 1, 2 =>
S^3 - 2*x*y*S - b2*S - b1 # = 0

### Eq 3:
(x-y)*S # = R
(S^2 - 4*x*y)*S^2 - R^2 # = 0
# 4*S^2*x*y = S^4 - R^2
# =>
2*S^4 - (S^4 - R^2) - 2*b2*S^2 - 2*b1*S # = 0
S^4 - 2*b2*S^2 - 2*b1*S + R^2 # = 0

### Solver:
solve.HtCompositeDiff.S2P4 = function(R, b, k=2, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	coeff = c(1, 0, -k*b[2], -k*b[1], R[1]^2);
	S = roots(coeff);
	if(debug) print(S);
	S = S[S != 0]; # exclude 0!
	xy.d = R[1] / S;
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	return(sol);
}

### Examples:

R = 2
b = c(1, 3)
sol = solve.HtCompositeDiff.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x^2 - y^2 # - R[1]

###
round0.p(poly.calc(x) * R[1]^2)


###########################

####################
### Diff Order 5 ###
####################

# x^4 + b2*y^2 + b1*y = Ru1
# y^4 + b2*x^2 + b1*x = Ru1
# x^5 - y^5 = 0

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^3 - 2*x*y*S - b2*S - b1 # = 0

### Eq 3:
# y = x * {m, m^2, m^3, m^4}, where m^5 = 1;
# =>
(m+1)^3*x^3 - 2*m*(m+1)*x^3 - b2*(m+1)*x - b1

### Solver:
solve.HtCompositeDiff.S2P4 = function(R, b, sort=TRUE, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	m = unity(5, all=T);
	m = m[-1]; # exclude 1;
	solve.P5 = function(id) {
		m = m[id];
		coeff = c((m+1)^3 - 2*m*(m+1), 0, - b[2]*(m+1), - b[1]);
		roots(coeff);
	}
	x = sapply(1:4, solve.P5);
	if(debug) print(x);
	x = as.vector(x);
	y = x * rep(m, each=3);
	sol = cbind(x=x, y=y);
	if(sort) sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = 0 # not used
b = c(1, 3)
sol = solve.HtCompositeDiff.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x^5 - y^5 # - R[1]

### P[12]
round0.p(poly.calc(x))


######################

######################
### Full Composite ###
######################

####################
### Order: 3 & 4 ###
####################

# x^3 + b1*y = Ru1
# y^3 + b1*x = Ru1
# x^4 + b2*y = Ru2
# y^4 + b2*x = Ru2

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b1 # = 0
# x*y = S^2 - b1;

### Diff Eqs 3, 4 =>
S*(S^2 - 2*x*y) - b2 # = 0

### =>
S*(S^2 - 2*b1) + b2 # = 0
S^3 - 2*b1*S + b2 # = 0

### Solver:
solve.HtComposite2.S2P34 = function(b, debug=TRUE) {
	coeff = c(1, 0, -2*b[1], b[2]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	return(sol);
}

### Examples:

b = c(1, -2)
sol = solve.HtComposite2.S2P34(b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^3 + b[1]*y, y^3 + b[1]*x)
cbind(x^4 + b[2]*y, y^4 + b[2]*x)

###
round0.p(poly.calc(x))
err = 3 - 2*x + x^2 + 4*x^3 - x^4 + x^6
round0(err)

