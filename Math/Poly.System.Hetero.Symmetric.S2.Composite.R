########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S2:
### Mixed Type: Composite
###
### draft v.0.2b


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


### draft v.0.2b:
# - unknown coefficients, variant:
#   Bu + b3*S = R2;
# - P[6] examples:
#   -4 - 10*x - 3*x^2 - x^3 + x^5 + x^6 = 0;
### draft v.0.2a:
# - unknown coefficients:
#   x^3 + Bu*x*y + b2*y^2 + b1*y = R1;
#   Bu * S = R2;
### draft v.0.1g - v.0.1h:
# - inter-connected system:
#   b3*x*y + Ru = R; [v.0.1g]
#   x*y*S + b3*Ru = R; [v.0.1h]
# - parametric P[6]; [v.0.1g-poly]
# - examples of P[6]:
#   2 + 6*x + 6*x^2 + 4*x^4 + x^6 = 0;
#   32*(20 + 6*x + 3*x^2) + x^6 = 0;
#   -375 + 80*x + 65*x^2 - 2*x^5 + x^6 = 0;
### draft v.0.1f:
# - other sums of powers:
#   x^3 + y^3 = R*(x^2 + y^2);
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

#################
### Variant:  ###
### Power Sum ###
#################

# x^4 + b2*y^2 + b1*y = Ru1
# y^4 + b2*x^2 + b1*x = Ru1
# x^3 + y^3 = R*(x^2 + y^2)

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^3 - 2*x*y*S - b2*S - b1 # = 0

### Eq 3:
S^3 - 3*x*y*S - R*S^2 + 2*R*x*y # = 0
# (3*S - 2*R)*x*y = S^3 - R*S^2
# =>
(3*S - 2*R)*S^3 - 2*(3*S - 2*R)*x*y*S - b2*(3*S - 2*R)*S - b1*(3*S - 2*R) # = 0
S^4 - 3*b2*S^2 + 2*b2*R*S - 3*b1*S + 2*b1*R # = 0


### Solver:
solve.HtCompositePows.S2P4 = function(R, b, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	coeff = c(1, 0, -3*b[2], 2*b[2]*R[1] - 3*b[1], 2*b[1]*R[1]);
	S = roots(coeff);
	if(debug) print(S);
	xy = (S^3 - R*S^2) / (3*S - 2*R[1]);
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = -1
b = c(1, 3)
sol = solve.HtCompositePows.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
cbind(x^4 + b[2]*y^2 + b[1]*y, y^4 + b[2]*x^2 + b[1]*x)
x^3 + y^3 - R[1]*(x^2 + y^2) # - R[1]

###
round0.p(poly.calc(x) * 2^2)
err = -1 - 4*x - 8*x^2 + 24*x^3 + 49*x^4 + 10*x^5 - 15*x^6 + 4*x^8
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


#######################
#######################

#######################
### Inter-Connected ###
#######################

# x^3 + b2*y^2 + b1*y = Ru
# y^3 + b2*x^2 + b1*x = Ru
# b3*x*y + Ru = R

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### Sum =>
x^3 + y^3 + b2*(x^2+y^2) + b1*(x+y) - 2*Ru # = 0
S^3 - 3*x*y*S + b2*(S^2 - 2*x*y) + b1*S - 2*Ru # = 0
# =>
S^3 - b2*S^2 - 2*b1*S - b2^2*S - b1*b2 + Ru # = 0
S^3 - b2*S^2 - 2*b1*S - b2^2*S - b1*b2 + R - b3*x*y
S^3 - (b2+b3)*S^2 - (2*b1 + b2^2 - b2*b3)*S + R - b1*(b2 - b3)


### Solver:
solve.HtCompositeDep.S2P4 = function(R, b, debug=TRUE) {
	if(length(b) < 2) stop("At least b1 and b3 are needed!")
	coeff = c(1, - (b[2]+b[3]), - (2*b[1] + b[2]^2 - b[2]*b[3]),
		- b[1]*(b[2] - b[3]) + R[1]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = 1
b = c(-2, 1, -1)
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
Ru  = x^3 + b[2]*y^2 + b[1]*y;
Ru2 = y^3 + b[2]*x^2 + b[1]*x;
round0(Ru - Ru2)
b[3]*x*y + Ru # - R[1]

###
round0.p(poly.calc(x))
err = 14 + 3*x + 6*x^2 - 6*x^3 + 4*x^4 + x^6
round0(err)

### parametric P[6]
b1 = b[1]; b2 = b[2]; b3 = b[3];
x^6 - (b2+b3)*x^5 + (- b1 + b2^2 + b3^2)*x^4 +
	- (2*R - 2*b1*(b2+b3) + b2*(b2-b3)^2)*x^3 +
	+ (R*(b2+b3) + b1^2 + b2^4 - b1*b2*b3 - b2^3*b3 - b1*b3^2)*x^2 +
	+ (R*b1 - R*b2*b3 + b1*b2^3 - 2*b1^2*b3 - b1*b2^2*b3)*x +
	+ R^2 - 3*R*b1*b2 - R*b2^3 - b1^3


###
b = c(-2, 1, -1)
sapply(-4:4, function(x) print(round0.p(
	poly.calc(solve.HtCompositeDep.S2P4(x, b, debug=F)[,1]))))

### Ex 2:
R = -2
b = c(-2, 1, -1)
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

round0.p(poly.calc(x))
err = 2 + 6*x + 6*x^2 + 4*x^4 + x^6
round0(err)


#########
### Ex 3:
b = c(0, 2, -1);
b[1] = b[2]^2 + b[3]^2;
R = b[1]*(b[2]+b[3]) - (b[2]*(b[2]-b[3])^2) / 2;
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

round0.p(poly.calc(x))
err = 43 + 82*x + 50*x^2 - x^5 + x^6
round0(err)


#########
### Ex 4:
b = c(0, 2, -2);
b[1] = b[2]^2 + b[3]^2;
R = b[1]*(b[2]+b[3]) - (b[2]*(b[2]-b[3])^2) / 2;
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

round0.p(poly.calc(x))
err = 32*(20 + 6*x + 3*x^2) + x^6
round0(err)


#######################

### Variant

# x^3 + b2*y^2 + b1*y = Ru
# y^3 + b2*x^2 + b1*x = Ru
# x*y*S + b3*Ru = R

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### Eq 3 =>
x*y*S + b3*Ru - R # = 0
S^3 - b2*S^2 - b1*S + b3*Ru - R # = 0

### Sum =>
S^3 - 3*x*y*S + b2*(S^2 - 2*x*y) + b1*S - 2*Ru # = 0
S^3 - b2*S^2 - 2*b1*S - b2^2*S - b1*b2 + Ru # = 0
# =>
b3*S^3 - b2*b3*S^2 - 2*b1*b3*S - b2^2*b3*S - b1*b2*b3 + b3*Ru # = 0
(b3-1)*S^3 - b2*(b3-1)*S^2 - b1*(2*b3-1)*S - b2^2*b3*S - b1*b2*b3 + R # = 0


### Solver:
solve.HtCompositeDep.S2P4 = function(R, b, debug=TRUE) {
	if(length(b) < 2) stop("At least b1 and b3 are needed!")
	coeff = c(b[3]-1, - b[2]*(b[3]-1), - b[1]*(2*b[3]-1) - b[2]^2*b[3],
		- b[1]*b[2]*b[3] + R[1]);
	S = roots(coeff);
	if(debug) print(S);
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	sol = cbind(x=x, y=y);
	sol = rbind(sol, sol[,2:1]);
	sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = 1
b = c(-2, 1, 2)
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
Ru  = x^3 + b[2]*y^2 + b[1]*y;
Ru2 = y^3 + b[2]*x^2 + b[1]*x;
round0(Ru - Ru2)
x*y*(x+y) + b[3]*Ru # - R[1]

###
round0.p(poly.calc(x))
err = 43 - x + 9*x^2 - 10*x^3 + 2*x^4 - x^5 + x^6
round0(err)

### parametric P[6]:
b1 = b[1]; b2 = b[2]; b3 = b[3];
(b3-1)^2*x^6 - b2*(b3-1)^2*x^5 +
	+ (b3-1)*(2*b1 + b3*(b2^2 - b1))*x^4 +
	- (b3-1)*(2*R + b1*b2 - 2*b1*b2*b3 + b2^3*b3)*x^3 +
	+ (R*b2*(b3-1) - b1^2*(b3-1) - 2*b1*b2^2*b3*(b3-1) + (b1+b2^2)^2*b3^2)*x^2 +
	+ (R*b1*(b3-1) - R*(b1+b2^2) + b1^2*b2*b3 + b1*b2^3*b3^2)*x +
	+ R^2 - 3*R*b1*b2*b3 - R*b2^3*b3 - b1^3*b3^2


###
R = -1
b = c(-2, 1, 2)
sapply(-4:4, function(x) print(round0.p(
	poly.calc(solve.HtCompositeDep.S2P4(R, c(x, b[2], b[3]), debug=F)[,1]))))

#########
### Ex 2:
R = -1
b = c(0, 1, 2)
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
Ru  = x^3 + b[2]*y^2 + b[1]*y;
Ru2 = y^3 + b[2]*x^2 + b[1]*x;
round0(Ru - Ru2)
x*y*(x+y) + b[3]*Ru # - R[1]

###
round0.p(poly.calc(x))
err = 3 + x + 3*x^2 + 2*x^4 - x^5 + x^6
round0(err)


#########
### Ex 3:
b = c(3, 1, NA); # b = c(5, 2, NA);
b[3] = -2*b[1] / (b[2]^2 - b[1])
R = - (b[1]*b[2] - 2*b[1]*b[2]*b[3] + b[2]^3*b[3])/2;
sol = solve.HtCompositeDep.S2P4(R, b)
x = sol[,1]; y = sol[,2];

### Test
Ru  = x^3 + b[2]*y^2 + b[1]*y;
Ru2 = y^3 + b[2]*x^2 + b[1]*x;
round0(Ru - Ru2)
x*y*(x+y) + b[3]*Ru # - R[1]

###
round0.p(poly.calc(x) * (b[3]-1)^2)
err = -387 + 66*x + 102*x^2 - 4*x^5 + 4*x^6
round0(err)


############################
############################

############################
### Unknown Coefficients ###
############################

# x^3 + Bu*x*y + b2*y^2 + b1*y = R1
# y^3 + Bu*x*y + b2*x^2 + b1*x = R1
# Bu + b3*S = R2

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### Sum =>
S^3 - 3*x*y*S + 2*Bu*x*y + b2*(S^2 - 2*x*y) + b1*S - 2*R1 # = 0
S^3 - b2*S^2 - 2*b1*S - b2^2*S - Bu*(S^2 - b2*S - b1) - b1*b2 + R1 # = 0
# =>
S^3 - b2*S^2 - 2*b1*S - b2^2*S - (R2 - b3*S)*(S^2 - b2*S - b1) - b1*b2 + R1 # = 0
(b3+1)*S^3 - (R2 + b2 + b2*b3)*S^2 - 2*b1*S - b1*b3*S - b2^2*S + b2*R2*S +
	- b1*b2 + b1*R2 + R1 # = 0


### Solver:
solve.HtCompositeDep.S2P3 = function(R, b, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	if(length(R) < 2) stop("R-parameter: 2 values needed!");
	b31 = b[3] + 1;
	coeff = c(b31, - (R[2] + b[2]*b31), - b[1] - b[1]*b31 - b[2]^2 + b[2]*R[2],
		- b[1]*b[2] + b[1]*R[2] + R[1]);
	S = roots(coeff);
	if(debug) print(S);
	# S = S[S != 0]; # exclude 0;
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	Bu = R[2] - b[3]*S;
	sol = cbind(x=x, y=y, Bu=Bu);
	sol = rbind(sol, sol[,c(2:1, 3)]);
	sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = c(0, -1)
b = c(-1, 2, -2)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

### Test
x^3 + Bu*x*y + b[2]*y^2 + b[1]*y # - R[1]
y^3 + Bu*x*y + b[2]*x^2 + b[1]*x # - R[1]
Bu + b[3]*(x+y) # - R[2]

###
round0.p(poly.calc(x) * abs(b[3]+1))
err = 1 - 6*x + 6*x^2 + 6*x^3 - 3*x^5 + x^6
round0(err)


#########
### Ex 2:
R = c(1, -1)
b = c(1, -2, -2)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

### Test
x^3 + Bu*x*y + b[2]*y^2 + b[1]*y # - R[1]
y^3 + Bu*x*y + b[2]*x^2 + b[1]*x # - R[1]
Bu + b3*(x+y) # - R[2]

###
round0.p(poly.calc(x) * abs(b[3]+1))
err = 14 - 19*x + 21*x^2 - 4*x^3 - 6*x^4 + x^5 + x^6
round0(err)


#########
### Ex 3:
R = c(1, 2)
b = c(NA, -2, -2); b[1] = (R[2]^2 + b[2]^2 + b[2]^2*b[3]) / (1 + 3*b[3] + 2*b[3]^2)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

###
round0.p(poly.calc(x) * abs(b[3]+1))
err = 9 - 4*x + 36*x^2 + 2*x^3 + 4*x^5 + x^6
round0(err)


#########
### Ex 4:
R = c(1, 2)
b = c(NA, 1, -2); b[1] = (R[2]^2 + b[2]^2 + b[2]^2*b[3]) / (1 + 3*b[3] + 2*b[3]^2)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

###
round0.p(poly.calc(x) * abs(b[3]+1))
err = -4 - 10*x - 3*x^2 - x^3 + x^5 + x^6
round0(err)


#########
### Ex 4:
R = c(NA, 2)
b = c(NA, 1, -2); b[1] = (R[2]^2 + b[2]^2 + b[2]^2*b[3]) / (1 + 3*b[3] + 2*b[3]^2)
R[1] = 3/2;
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

###
round0.p(poly.calc(x) * abs(b[3]+1) * 4)
err = -19 - 54*x - 10*x^2 + 4*x^5 + 4*x^6
round0(err)


### P[6]
R1 = R[1]; R2 = R[2]; b1 = b[1]; b2 = b[2]; b3 = b[3];
(b3+1)^2*x^6 - (R2 + b2 + R2*b3 + 2*b2*b3 + b2*b3^2)*x^5 +
	+ (R2^2 - b1 + b2^2 - 3*b1*b3 + b2^2*b3 - 2*b1*b3^2)*x^4 +
	+ (- 2*R1*(b3+1) + 2*R2*b1 - R2^2*b2 + 2*b1*b2 + 2*R2*b2^2 - b2^3 +
		+ R2*b1*b3 + 3*b1*b2*b3 + R2*b2^2*b3 - b2^3*b3 + b1*b2*b3^2)*x^3 +
	+ (R1*R2 - R2^2*b1 + b1^2 + R1*b2 - R2*b1*b2 - R2*b2^3 + b2^4 +
		+ b1^2*b3 + R1*b2*b3 + R2*b1*b2*b3 - 2*b1*b2^2*b3 + b1^2*b3^2)*x^2 +
	+ (R1*b1 - 2*R2*b1^2 - R1*R2*b2 - R2*b1*b2^2 + b1*b2^3 + 2*R1*b1*b3 - b1^2*b2*b3 +
		+ R1*b2^2*b3)*x +
	+ R1^2 - b1^3 - 3*R1*b1*b2 - R1*b2^3


####################

### Variant:

# x^3 + Bu*x*y + b2*y^2 + b1*y = R1
# y^3 + Bu*x*y + b2*x^2 + b1*x = R1
# Bu*S = R2

### Solution:

### Case: x != y

### Diff Eqs 1, 2 =>
S^2 - x*y - b2*S - b1 # = 0
# x*y = S^2 - b2*S - b1;

### Sum =>
S^3 - 3*x*y*S + 2*Bu*x*y + b2*(S^2 - 2*x*y) + b1*S - 2*R1 # = 0
S^3 - b2*S^2 - 2*b1*S - b2^2*S - Bu*(S^2 - b2*S - b1) - b1*b2 + R1 # = 0
S^4 - b2*S^3 - 2*b1*S^2 - b2^2*S^2 - R2*S^2 + R2*b2*S - b1*b2*S + R1*S + R2*b1 # = 0

### Solver:
solve.HtCompositeDep.S2P3 = function(R, b, debug=TRUE) {
	if(length(b) < 2) b = c(b, 0);
	if(length(R) < 2) stop("R-parameter: 2 values needed!")
	coeff = c(1, - b[2], - (2*b[1] + b[2]^2 + R[2]),
		R[2]*b[2] - b[1]*b[2] + R[1], R[2]*b[1]);
	S = roots(coeff);
	if(debug) print(S);
	S = S[S != 0]; # exclude 0;
	xy = S^2 - b[2]*S - b[1];
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2;
	y = S - x;
	Bu = R[2] / S;
	sol = cbind(x=x, y=y, Bu=Bu);
	sol = rbind(sol, sol[,c(2:1, 3)]);
	sol = sort.sol(sol);
	return(sol);
}

### Examples:

R = c(0, -1)
b = c(-1, 2)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

### Test
x^3 + Bu*x*y + b[2]*y^2 + b[1]*y # - R[1]
y^3 + Bu*x*y + b[2]*x^2 + b[1]*x # - R[1]
Bu*(x+y) # - R[2]

###
round0.p(poly.calc(x))
err = 1 - 10*x + 20*x^2 - 14*x^3 + 10*x^4 - 8*x^5 + 5*x^6 - 2*x^7 + x^8
round0(err)


#########
### Ex 2:
R = c(2, -1)
b = c(0, -1)
sol = solve.HtCompositeDep.S2P3(R, b)
x = sol[,1]; y = sol[,2]; Bu = sol[,3];

### Test
x^3 + Bu*x*y + b[2]*y^2 + b[1]*y # - R[1]
y^3 + Bu*x*y + b[2]*x^2 + b[1]*x # - R[1]
Bu*(x+y) # - R[2]

###
round0.p(poly.calc(x))
err = 9 + 3*x - 3*x^2 - 6*x^3 + x^5 + x^6
round0(err)

