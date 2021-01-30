########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Mixt Variable
### Hetero-Symmetric S2 + Symmetric
###
### draft v.0.1e


### Mixt Polynomial Systems:
### Hetero-Symmetric S2
### + Symmetric 3rd Var

### Example 1:
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R1
# x^n + y^n = R2

### Example 2:
# (x + d)^n + y^n = R1
# x^n + (y + d)^n = R1
# (x - d)^n + (y - d)^n = R2

### Example 3:
# x^n + d*y = R1
# y^n + d*x = R1
# Variants: e.g. x*y*d = R2


###############
### History ###
###############


### draft v.0.1e:
# - Eq 3 variant for the Mixt-Symmetric S3P3;
### draft v.0.1c - v.0.1d:
# - Mixt Hetero-Symmetric system: Order 3;
#  -- first/basic variant: x*y; (=> generalized symmetric P[6])
#  -- variant: x*y*d (=> P[12]);
### draft v.0.1b - v.0.1b-var:
# - Mixt Hetero-Symmetric system: Order 2;
# - more variants for Eq 3;
### draft v.0.1a:
# - moved to new file
#   from Poly.System.Asymmetric.S2.R;
### [previous file]
# - some initial work on the Mixt Hetero-Symmetric S3P3:
#   (x+d)^3 + y^3 = R;
# - Mixt Symmetric variants (in v.0.2a-sym);


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;


##########################

##########################
### Polynomial Systems ###
##########################

########################
### Hetero-Symmetric ###
########################

###############
### Order 2 ###
###############

# x^2 + d*y = R1
# y^2 + d*x = R1
# x*y = R2

# see next sections for various variants of Eq 3;

### Solution:

### Sum Eq 1 + Eq 2 =>
x^2 + y^2 + d*(x + y) - 2*R1 # = 0
S^2 + d*S - 2*x*y - 2*R1 # = 0
# =>
S^2 + d*S - 2*R1 - 2*R2 # = 0
# d*S = - S^2 + 2*R1 + 2*R2

### Diff =>
x^2 - y^2 - d*(x - y) # = 0
(x - y)*(x + y - d) # = 0
(x - y)*(S - d) # = 0
### Case: x != y =>
S - d # = 0
# d = S

### =>
2*S^2 - 2*R1 - 2*R2 # = 0
S^2 - R1 - R2 # = 0

### Solver:
solve.Ht2V.S3P2 = function(R, debug=TRUE) {
	coeff = c(1, 0, - R[1] - R[2])
	S = roots(coeff)
	if(debug) print(S);
	d = S;
	R2 = R[2];
	xy.d = sqrt(S^2 - 4*R2 + 0i);
	x = (S + xy.d) / 2
	y = (S - xy.d) / 2
	sol = cbind(as.vector(x), as.vector(y), as.vector(d));
	sol = rbind(sol, sol[,c(2,1,3)]);
	return(sol)
}

### Examples:

R = c(1, 1)
sol = solve.Ht2V.S3P2(R)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
x^2 + d*y # - R[1]
y^2 + d*x # - R[1]
x*y # - R[2]

##################

### Variant:
# x^2 + d*y = R1
# y^2 + d*x = R1
# x*y*(x + y) = R2

### Solution:

### Sum Eq 1 + Eq 2 =>
x^2 + y^2 + d*(x + y) - 2*R1 # = 0
S^2 + d*S - 2*x*y - 2*R1 # = 0
# =>
S^3 + d*S^2 - 2*R1*S - 2*R2

### Diff =>
### Case: x != y =>
S - d # = 0
# d = S

### =>
2*S^3 - 2*R1*S - 2*R2 # = 0
S^3 - R1*S - R2

### Solver:
solve.Ht2V.S3P2 = function(R, debug=TRUE) {
	coeff = c(1, 0, - R[1], - R[2])
	S = roots(coeff)
	if(debug) print(S);
	d = S;
	R2 = R[2];
	xy = R2 / S;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2
	y = (S - xy.d) / 2
	sol = cbind(as.vector(x), as.vector(y), as.vector(d));
	sol = rbind(sol, sol[,c(2,1,3)]);
	return(sol)
}

### Examples:

R = c(3, -1)
sol = solve.Ht2V.S3P2(R)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
x^2 + d*y # - R[1]
y^2 + d*x # - R[1]
x*y*(x + y) # - R[2]


##################

### Variant:
# x^2 + d*y = R1
# y^2 + d*x = R1
# x*y*d = R2

### Solution:

### Sum Eq 1 + Eq 2 =>
x^2 + y^2 + d*(x + y) - 2*R1 # = 0
S^2 + d*S - 2*x*y - 2*R1 # = 0
# =>
d*S^2 + d^2*S - 2*R1*d - 2*R2

### Diff =>
### Case: x != y =>
S - d # = 0
# d = S

### =>
2*S^3 - 2*R1*S - 2*R2
S^3 - R1*S - R2
# same as previous case;
# Note: S == d (for order 2);


##################

### Variant:
# x^2 + d*y = R1
# y^2 + d*x = R1
# x^2 + y^2 + d^2 = R2

### Solution:

### Diff =>
### Case: x != y =>
S - d # = 0
# d = S

### Sum Eq 1 + Eq 2 =>
x^2 + y^2 + d*(x + y) - 2*R1 # = 0
S^2 + d*S - 2*x*y - 2*R1 # = 0
### =>
2*S^2 - 2*x*y - 2*R1 # = 0
S^2 - x*y - R1 # = 0
# x*y = S^2 - R1

### Eq 3 =>
S^2 - 2*x*y + d^2 - R2 # = 0
2*S^2 - 2*x*y - R2
2*R1 - R2 # = 0
# - NO direct solutions, only for extensions;

### Solver:
solve.Ht2V.S3P2 = function(R, b1.ext=0, b2.ext=0, debug=TRUE) {
	len = max(length(b1.ext), length(b2.ext));
	b1.ext = c(rep(0, len - length(b1.ext)), rev(b1.ext));
	b2.ext = c(rep(0, len - length(b2.ext)), rev(b2.ext));
	b = - 2*b1.ext + b2.ext;
	if(all(b == 0)) stop("NO solution possible!")
	coeff = c(b, 2*R[1] - R[2])
	S = roots(coeff)
	if(debug) print(S);
	d = S;
	R1 = R[1] - sapply(S, function(S) sum(rev(b1.ext)*S^seq(len)));
	xy = S^2 - R1;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2
	y = (S - xy.d) / 2
	sol = cbind(as.vector(x), as.vector(y), as.vector(d));
	sol = rbind(sol, sol[,c(2,1,3)]);
	return(sol)
}

### Examples:

R = c(3, -1)
b1.ext = c(1,-1, 1)
sol = solve.Ht2V.S3P2(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
S = x+y; len = length(b1.ext); ext1 = sapply(S, function(S) sum(b1.ext*S^seq(len)));
x^2 + d*y + ext1 # - R[1]
y^2 + d*x + ext1 # - R[1]
x^2 + y^2 + d^2 # - R[2]

### Classic Polynomial
round0.p(poly.calc(x)) * 8


##################

### Variant:
# x^2 + d*y = R1
# y^2 + d*x = R1
# x^2 + y^2 = R2

### Solution:

### Diff =>
### Case: x != y =>
S - d # = 0
# d = S

### Sum Eq 1 + Eq 2 =>
S^2 - x*y - R1 # = 0
# x*y = S^2 - R1

### Eq 3 =>
S^2 - 2*x*y - R2 # = 0
# =>
x*y - R1 + R2 # = 0
S^2 - 2*R1 + R2 # = 0


### Solver:
solve.Ht2V.S3P2 = function(R, b1.ext=0, b2.ext=0, debug=TRUE) {
	len = max(2, length(b1.ext), length(b2.ext));
	b1.ext = c(rep(0, len - length(b1.ext)), rev(b1.ext));
	b2.ext = c(rep(0, len - length(b2.ext)), rev(b2.ext));
	b = 2*b1.ext - b2.ext;
	coeff = c(rep(0, len - 2), 1, 0, -2*R[1] + R[2]) + c(b, 0)
	S = roots(coeff)
	if(debug) print(S);
	d = S;
	R1 = R[1] - sapply(S, function(S) sum(rev(b1.ext)*S^seq(len)));
	xy = S^2 - R1;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d) / 2
	y = (S - xy.d) / 2
	sol = cbind(as.vector(x), as.vector(y), as.vector(d));
	sol = rbind(sol, sol[,c(2,1,3)]);
	return(sol)
}

### Examples:

R = c(3, -1)
b1.ext = c(1, -1, 1/2)
sol = solve.Ht2V.S3P2(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
S = x+y; len = length(b1.ext); ext1 = sapply(S, function(S) sum(b1.ext*S^seq(len)));
x^2 + d*y + ext1 # - R[1]
y^2 + d*x + ext1 # - R[1]
x^2 + y^2 # - R[2]

### Classic Polynomial
round0.p(poly.calc(x)) * 8
err = 37 + 8*x + 10*x^2 + 12*x^3 + 16*x^4 - 8*x^5 + 8*x^6
round0(err)

####################
####################

###############
### Order 3 ###
###############

### Hetero-Symmetric

# x^3 + d*y = R1
# y^3 + d*x = R1
# x*y = R2

### Diff =>
x^3 - y^3 - d*(x - y) # = 0
(x - y)*(x^2 + y^2 + x*y - d) # = 0
(x - y)*(S^2 - x*y - d) # = 0
### Case: x != y =>
S^2 - R2 - d # = 0
# d = S^2 - R2

### Sum =>
x^3 + y^3 + d*(x + y) - 2*R1 # = 0
S^3 - 3*x*y*S + d*S - 2*R1
S^3 - 3*R2*S + (S^2 - R2)*S - 2*R1
2*S^3 - 4*R2*S - 2*R1
S^3 - 2*R2*S - R1

### Solver:
solve.MxHtSym.S3P3 = function(R, b1.ext=0, b2.ext=0, debug=TRUE) {
	len = max(3, length(b1.ext), 1 + length(b2.ext));
	b1m.ext = c(rep(0, len - length(b1.ext)), rev(b1.ext));
	b2m.ext = c(rep(0, len - length(b2.ext) - 1), rev(b2.ext));
	coeff = c(rep(0, len - 3), 1, 0, -2*R[2], -R[1]) +
		+ c(b1m.ext, 0) + c(2*b2m.ext, 0, 0)
	S = roots(coeff)
	if(debug) print(S);
	#
	len1 = length(b1.ext); len2 = length(b2.ext);
	# R1 = R[1] - sapply(S, function(S) sum(b1.ext * S^seq(len1)))
	R2 = R[2] - sapply(S, function(S) sum(b2.ext * S^seq(len2)));
	d = S^2 - R2;
	xy = R2;
	#
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), d=as.vector(d))
	sol = rbind(sol, sol[,c(2,1,3)])
	return(sol)
}
test.MxHt.S3P3 = function(sol, R, b1.ext=0, b2.ext=0, type="x*y") {
	x = sol[,1]; y = sol[,2]; d = sol[,3]; S = x+y;
	lenPow1 = length(b1.ext); ext1 = sapply(S, function(S) sum(b1.ext * S^seq(lenPow1)))
	lenPow2 = length(b2.ext); ext3 = sapply(S, function(S) sum(b2.ext * S^seq(lenPow2)))
	err1 = x^3 + d*y + ext1 # - R[1]
	err2 = y^3 + d*x + ext1 # - R[1]
	type = match(type, c("x*y", "xyd", "xyS"))
	if(is.na(type)) { stop("Type NOT yet supported!"); }
	else if(type == 1) {
		err3 = x*y # - R[2]
	} else if(type == 2) {
		err3 = x*y * d;
	} else if(type == 3) {
		err3 = x*y * S;
	}
	err3 = err3 + ext3;
	err = rbind(err1, err2, err3)
	return(round0(err))
}

### Examples:

R = c(-1, 2)
b1.ext = c(0)
#
sol = solve.MxHtSym.S3P3(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
test.MxHt.S3P3(sol, R, b1.ext=b1.ext)

### Classic Polynomial
round0.p(poly.calc(x))
err = 8 + 4*x^2 + x^3 + 2*x^4 + x^6
round0(err)

#########
### Ex 2:
R = c(-1, -1)
b1.ext = c(0)
#
sol = solve.MxHtSym.S3P3(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
test.MxHt.S3P3(sol, R, b1.ext=b1.ext)

### Classic Polynomial
round0.p(poly.calc(x))
err = -1 + x^2 + x^3 - x^4 + x^6
round0(err)


#########
### Ex 3: generalized symmetric polynom
R = c(1, -2)
b1.ext = c(-1, 1, -2)
#
sol = solve.MxHtSym.S3P3(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
test.MxHt.S3P3(sol, R, b1.ext=b1.ext)

### Classic Polynomial
round0.p(poly.calc(x))
err = -8 - 4*x + 18*x^2 + 5*x^3 - 9*x^4 - x^5 + x^6
round0(err)

################
### Variants ###

# x*y*d = R2

### Solution

### Diff =>
### Case: x != y =>
S^2 - x*y - d # = 0
# d = S^2 - x*y

### Sum =>
x^3 + y^3 + d*(x + y) - 2*R1 # = 0
S^3 - 3*x*y*S + d*S - 2*R1
# =>
S^3 - 3*x*y*S + (S^2 - x*y)*S - 2*R1
S^3 - 2*x*y*S - R1
# 2*x*y*S = S^3 - R1

### Eq 3:
x*y*(S^2 - x*y) - R2 # = 0
x*y*S^2 - (x*y)^2 - R2 # = 0
4*x*y*S^4 - 4*(x*y)^2*S^2 - 4*R2*S^2 # = 0
2*(S^3 - R1)*S^3 - (S^3 - R1)^2 - 4*R2*S^2 # = 0
S^6 - 4*R2*S^2 - R1^2 # = 0

### Solver:
solve.MxHtSym.S3P3 = function(R, b1.ext=0, b2.ext=0, debug=TRUE) {
	len = max(6, 2*length(b1.ext), 2 + length(b2.ext));
	b1m.ext = c(rep(0, len - length(b1.ext)), rev(b1.ext));
	b1sq.ext = c(rep(0, len - 2*length(b1.ext)), rev(b1.ext)^2);
	b2m.ext = c(rep(0, len - length(b2.ext) - 2), rev(b2.ext));
	# TODO: cross-terms for ()^2;
	coeff = c(rep(0, len - 6), 1, 0,0,0, -4*R[2], 0, -R[1]^2) +
		+ c(2*R[1]*b1m.ext, 0) - c(b1sq.ext, 0, 0) +
		+ c(4*b2m.ext, 0, 0, 0)
	S = roots(coeff)
	if(debug) print(S);
	#
	len1 = length(b1.ext); len2 = length(b2.ext);
	R1 = R[1] - sapply(S, function(S) sum(b1.ext * S^seq(len1)));
	# R2 = R[2] - sapply(S, function(S) sum(b2.ext * S^seq(len2)));
	xy = (S^3 - R1) / 2 / S;
	d = S^2 - xy;
	#
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), d=as.vector(d))
	sol = rbind(sol, sol[,c(2,1,3)])
	return(sol)
}

### Examples:

R = c(-1, 2)
b1.ext = c(0)
#
sol = solve.MxHtSym.S3P3(R, b1.ext=b1.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
test.MxHt.S3P3(sol, R, b1.ext=b1.ext, type="xyd")

### Classic Polynomial
round0.p(poly.calc(x[c(1,2,4,7)]))
err = 2 + x + x^4
round0(err)

### Ex 2:
R = c(-1, 2)
b1.ext = c(-1)
b2.ext = c(0)
#
sol = solve.MxHtSym.S3P3(R, b1.ext=b1.ext, b2.ext=b2.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
test.MxHt.S3P3(sol, R, b1.ext=b1.ext, b2.ext=b2.ext, type="xyd")

### Classic Polynomial
round0.p(poly.calc(x)) # P[12]


####################
####################

###############
### Order 3 ###
###############


### S3P3:
# (x + d)^3 + y^3 = R1
# x^3 + (y + d)^3 = R1
# x^3 + y^3 = R2

### Variants Eq 3:
# x*y = R2
# x*y*(x + y) = R2
# x*y*d = R2
# x^2 + y^2 = R2
# (x - d)^3 + (y - d)^3 = R2

### Solution:
# - trivial system;
# - only Cases:
#   d = -S OR d = 0 OR x = y;
#   [independent of 3rd Eq]

### Sum Eq 1 + Eq 2:
(x + d)^3 + (y + d)^3 + x^3 + y^3 - 2*R1 # = 0
2*(x^3 + y^3) + 3*d*(x^2 + y^2) + 3*d^2*(x + y) + 2*d^3 - 2*R1 # = 0
2*R2 + 3*d*(S^2 - 2*x*y) + 3*d^2*S + 2*d^3 - 2*R1 # = 0
3*d*S^2 - 6*d*x*y + 3*d^2*S + 2*d^3 - 2*R1 + 2*R2 # = 0

### Diff: Eq 1 - Eq 2
(x - y)*((x+d)^2 + (y+d)^2 + (x+d)*(y+d) - x^2 - y^2 - x*y) # = 0
(x - y)*(3*d*S + 3*d^2) # = 0
### Case: x != y =>
# d = 0 OR
# d = -S
### =>
3*x*y*S - S^3 - R1 + R2 # = 0
# 3*x*y*S = S^3 + R1 - R2

### Eq 3 =>
S^3 - 3*x*y*S - R2 # = 0
- R1 # = 0
# has Solution only if R1 == 0
# OR various extensions;


### Solver:
solve.MxHtSym.S3P3 = function(R, b.ext=0, debug=TRUE) {
	if(all(b.ext == 0)) stop("NO solutions: x != y")
	S = roots(c(rev(b.ext), -R[1]))
	if(debug) print(S);
	#
	lenPow = length(b.ext)
	R1 = R[1] - sapply(S, function(S) sum(b.ext * S^seq(lenPow)))
	R2 = R[2];
	d = -S;
	xy = (S^3 + R1 - R2) / (3*S)
	#
	xy.d = sqrt(S^2 - 4*xy + 0i)
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	return(cbind(x=as.vector(x), y=as.vector(y), d=as.vector(d)))
}

### Examples:

R = c(-1, 2)
b.ext = c(1, -1)
#
sol = solve.MxHtSym.S3P3(R, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
lenPow = nrow(sol); ext = sapply(x+y, function(S) sum(b.ext * S^seq(lenPow)))
(x + d)^3 + y^3 + ext # - R[1]
x^3 + (y + d)^3 + ext # - R[1]
x^3 + y^3 # - R[2]


### Ex 2:
R = c(-1, 2)
b.ext = c(1, -2, 1)
#
sol = solve.MxHtSym.S3P3(R, b.ext=b.ext)
x = sol[,1]; y = sol[,2]; d = sol[,3];


###########################
###########################

#################
### Symmetric ###
#################

###############
### Order 3 ###
###############

### S3P3:
# (x + d)^3 + (y + d)^3 = R1
# (x - d)^3 + (y - d)^3 = R2
# x^3 + y^3 = R3

### Solution:

### TODO: robust!

### Diff: Eq 1 - Eq 2 =>
6*d*(x^2 + y^2) + 4*d^3 - R1 + R2 # = 0
6*d*(S^2 - 2*x*y) + 4*d^3 - R1 + R2
6*d*S^2 - 12*d*x*y + 4*d^3 - R1 + R2

### Sum: Eq 1 + Eq 2 - 2*Eq3 =>
6*d^2*(x + y) - R1 - R2 + 2*R3 # = 0
6*d^2*S - R1 - R2 + 2*R3
# 6*d^2*S = R1 + R2 - 2*R3

### Eq 3:
S^3 - 3*x*y*S - R3 # = 0
# 3*x*y*S = S^3 - R3

### =>
6*d*S^3 - 12*d*x*y*S + 4*d^3*S - (R1 - R2)*S
2*d*S^3 + 4*d*R3 + 4*d^3*S - (R1 - R2)*S
6*d*S^3 + 2*d*(R1 + R2 + 4*R3) - 3*(R1 - R2)*S
2*d*(3*S^3 + (R1 + R2 + 4*R3)) - 3*(R1 - R2)*S
12*d^2*S*(3*S^3 + (R1 + R2 + 4*R3))^2 - 27*(R1 - R2)^2*S^3
2*(R1 + R2 - 2*R3)*(3*S^3 + (R1 + R2 + 4*R3))^2 - 27*(R1 - R2)^2*S^3


### Solver:
solve.MSym.S3P3 = function(R, debug=TRUE) {
	Rd = (R[1] - R[2])
	Rd3 = (R[1] + R[2] - 2*R[3])
	Rs3 = (R[1] + R[2] + 4*R[3])
	coeff = c(18*Rd3, 0,0, 12*Rd3*Rs3 - 27*Rd^2, 0,0, 2*Rd3*Rs3^2)
	S = roots(coeff)
	#
	R1 = R[1]; R2 = R[2]; R3 = R[3];
	xy = (S^3 - R3) / (3*S);
	### d: robust!
	# d = -sqrt((R1 + R2 - 2*R3) / 6 / S + 0i);
	d2 = (R1 + R2 - 2*R3) /6 / S;
	d  = 3*(R1 - R2)*S /2 / (3*S^3 + (R1 + R2 + 4*R3));
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), d=as.vector(d))
	return(sol)
}

### Examples:

R = c(1, -1, 2)
sol = solve.MSym.S3P3(R)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]
x^3 + y^3 # - R[3]


#########
### Ex 2:
R = c(1, 2, -2)
sol = solve.MSym.S3P3(R)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]
x^3 + y^3 # - R[3]

### Classic Polynomial:
round0.p(poly.calc(c(x,y))) # ???
round0.p(poly.calc(d)) * 360 # trivial


#############
### Variants:

### Eq 3:
### x^2 + b^2 + b*d^2 = R3

### Eq 1 - Eq 2:
6*d*(x^2 + y^2) + 4*d^3 - R1 + R2 # = 0
6*d*(S^2 - 2*x*y) + 4*d^3 - R1 + R2 # = 0

### Eq 1 + Eq 2:
2*(x^3 + y^3) + 6*d^2*(x+y) - R1 - R2 # = 0
2*S^3 - 6*x*y*S + 6*d^2*S - R1 - R2 # = 0

### Eq 3:
# x^2 + y^2 = R3 - b*d^2
### Eq 1 - Eq 2 =>
6*d*(R3 - b*d^2) + 4*d^3 - R1 + R2 # = 0
### Eq:
(6*b - 4)*d^3 - 6*R3*d + R1 - R2 # = 0

# Eq 1-bis * S - Eq 2-bis * 3*d =>
6*d*x*y*S + (4*d^3 - R1 + R2)*S - 3*d*(6*d^2*S - R1 - R2) # = 0
6*d*x*y*S - 14*d^3*S - (R1 - R2)*S + 3*(R1 + R2)*d # = 0
# 6*d*x*y*S = 14*d^3*S + (R1 - R2)*S - 3*(R1 + R2)*d # = 0
# Eq 3 =>
3*d*S^3 - 6*d*x*y*S + 3*b*d^3*S - 3*R3*d*S # = 0
3*d*S^3 - (14*d^3*S + (R1 - R2)*S - 3*(R1 + R2)*d) + 3*b*d^3*S - 3*R3*d*S # = 0
3*d*S^3 + (3*b*d^3 - 14*d^3 - 3*R3*d - (R1 - R2))*S + 3*(R1 + R2)*d # = 0


### Solver:
solve.MSym.S3P3 = function(R, b=1, debug=TRUE) {
	Rd = R[1] - R[2];
	d = roots(c((6*b[1] - 4), 0, - 6*R[3], Rd))
	S = sapply(d, function(d) {
		d3 = d^3;
		coeff = c(3*d, 0, (3*b[1]*d3 - 14*d3 - 3*R[3]*d - Rd), 3*(R[1] + R[2])*d)
		S = roots(coeff)
	})
	#
	if(debug) { print(d); print(S); }
	d = matrix(d, ncol=length(d), nrow=3, byrow=TRUE)
	R1 = R[1]; R2 = R[2]; R3 = R[3];
	xy = (S^2 + b[1]*d^2 - R3) / 2;
	xy.d = sqrt(S^2 - 4*xy + 0i);
	x = (S + xy.d)/2;
	y = (S - xy.d)/2;
	sol = cbind(x=as.vector(x), y=as.vector(y), d=as.vector(d))
	return(sol)
}

### Examples:

R = c(1, -2, 2)
b = 1
sol = solve.MSym.S3P3(R, b=b)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]
x^2 + y^2 + b[1]*d^2 # - R[3]

round0.p(poly.calc(c(x,y))) * 32

### Ex 2:
R = c(1, -2, -1)
b = 1/2
sol = solve.MSym.S3P3(R, b=b)
x = sol[,1]; y = sol[,2]; d = sol[,3];

### Test
(x + d)^3 + (y + d)^3 # - R[1]
(x - d)^3 + (y - d)^3 # - R[2]
x^2 + y^2 + b[1]*d^2 # - R[3]

print(round0.p(poly.calc(c(x,y))) * 64, digits=16)
round0.p(poly.calc(c(x,y))) * 64*64


