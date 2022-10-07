########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### with Mixed Leading Term
###
### draft v.0.1i


### Heterogeneous Symmetric Polynomial Systems
# - cyclic permutation of roots;

### 2 Variables:
### Simple System:
# x^n*y^m + P(x, y) = R
# y^n*x^m + P(y, x) = R


###############
### Systems ###
###############


### Mixed-Leading Term: Order 2+1
# M21.1) x^2*y + b*x: NO solutions (x != y);
# M21.2) x^2*y + b*y: trivial;
# M21.3) x^2*y + b2*x^2 + b1*x: trivial (trivial P2);
# M21.4) x^2*y + b3*x*y + b2*x^2 + b1*x: trivial (trivial P2);
### Mixed: Order n+1:
# M31.1) x^3*y + b*x: trivial P4;
# M31.2) x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R; (P3 => P6; some nice)
# M41.1) x^4*y + b*x; (P5 => P10)
### Mixed: Order n+2:
# M42.1) x^4*y^2 + b*y = R; (P[5] => P[10])
# M42.2) x^4*y^2 + b2*y^2 + b1*y = R; (P[7] => P[14])
# M52.1) x^5*y^2 + b*y = R; (P[9] => P[18])
### Mixed: Order n+3:
# M43.1) x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R; (trivial P2; base P7)
# M43.2) x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
# M43.3) x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R; (P[3] => P[6])
# M53.1) x^5*y^3 + b*y = R; (P[7] => P[14])
# M53.2) x^5*y^3 + b2*y^2 + b1*y = R; (P[10] => P[20])


###############

###############
### History ###
###############


### v.0.1a:
# - started to move sections with Mixed Leading terms:
#   from Poly.System.Hetero.Symmetric.R
#   to this file;


######################
######################

### Helper Functions

# library(polynom)
# library(pracma)

source("Polynomials.Helper.R")

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

##########################
##########################

##########################
### Polynomial Systems ###
##########################

##########################
### Mixed Leading Term ###
##########################

#########################
### x^j*y^k + P(x, y) ###
#########################

### Variant: + b*x
### x^2*y + b*x

# x^2*y + b1*x = R
# y^2*x + b1*y = R

### Solution:
# *NO* solution x != y;

### Diff =>
# x*y = -b1

### Sum =>
(x+y)*(x*y + b1) - 2*R # = 0
# but: x*y + b1 = 0
# => NO solution (x != y);

### NO Solution!

### Extension: A1
# x^2*y + b1*x + be1*(x+y) = R
be1*S - 2*R # = 0



###################

### Variant: + b*y
### x^2*y + b*y

# x^2*y + b1*y = R
# y^2*x + b1*x = R

# very simple system

### Solution:

### Diff =>
# x*y = b1

### Sum =>
# Z = R / b1

### Solution:
solve.mx.S2P2 = function(R, b) {
	x.sum = R[1]/b[1]
	xy = b[1]
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	return(rbind(sol, sol[,2:1]))
}

### Example
R = 1
b = 3
#
sol = solve.mx.S2P2(R, b)
x = sol[,1]; y = sol[,2];
sol


### Test
x^2*y + b[1]*y
y^2*x + b[1]*x


#####################

### Variant: + b2*x^2

### x^2*y + b2*x^2 + b1*x
# x^2*y + b2*x^2 + b1*x = R
# y^2*x + b2*y^2 + b1*y = R

### Extension:
# x^2*y + b3*x*y + b2*x^2 + b1*x = R

### Solution:

### Diff =>
# x*y = -b2*Z - b1

### Sum =>
# Z = (R - b1*b2) / b2^2

### Extensions:
### E1: x^2*y + b3*x*y + b2*x^2 + b1*x = R
### Sum =>
(b2^2 - b2*b3)*Z + b1*b2 - b1*b3 - R # = 0


solve.ht21.S2P2 = function(b, R) {
	if(length(b) == 2) {
		r.sum = (R - b[1]*b[2]) / b[2]^2
	} else {
		r.sum = (R - b[1]*b[2] + b[1]*b[3]) / (b[2]^2 - b[2]*b[3])
	}
	xy = -b[2]*r.sum - b[1]
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example:
R = -1
b = c(1, 3)
#
sol = solve.ht21.S2P2(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y + b[2]*x^2 + b[1]*x
y^2*x + b[2]*y^2 + b[1]*y

### Classical Polynomial:
# trivial P2;
b[2]*x^2 - (R/b[2] - b[1])*x - R


### Example 2: Extension
R = -1
b = c(1, 3, 1)
#
sol = solve.ht21.S2P2(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y + b[3]*x*y + b[2]*x^2 + b[1]*x
y^2*x + b[3]*x*y + b[2]*y^2 + b[1]*y

### Classic Polynomial:
b[2]*(b[2] - b[3])*x^2 + (b[1]*(b[2] - b[3]) - R)*x - R*b[2]


############################
############################
############################

######################
### Leading: x^3*y ###
######################

### x^3*y + b3*(x*y)^2 + b2*x*y + b1*y

# x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^3*x + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

### Diff =>
# x*y = b1 / S;

### Sum =>
b1*S^3 - R*S^2 + b1*b2*S + b1^2*b3 - b1^2 # = 0


### Solver:

solve.htx3y = function(R, b, debug=TRUE) {
	if(R == 0 && b[1] == 0) return(list(sol=NA, p=NA))
	r.sum = roots(c(b[1], - R, b[1]*b[2], b[1]^2*b[3] - b[1]^2));
	r.sum = round0(r.sum)
	r.sum = r.sum[ r.sum != 0 ] # avoid division by 0
	if(debug) print(r.sum);
	xy = b[1] / r.sum
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

test.p = function(R, b) {
	b1 = b[1]; b2 = b[2]; b3 = b[3];
	#
	err = - b1^3 - (2*b1^2*b2 - R^2)*x +
		- b1*(2*R*b3 + b2^2 - R)*x^2 +
		- (- R*b2 + 2*b1^2 - b1^2*b3 - b1^2*b3^2)*x^3 +
		+ b1*b2*(b3 - 2)*x^4 - R*(b3 - 1)*x^5 +
		+ b1*(b3 - 1)*x^6;
	return(round0(err))
}

### Classic Polynomial:
# - for P[6]: see at the end of this section;
((1 + b[3])*x^4 + b[2]*x^2 + b[1]*x - R) * P6


### Example:
R = 1
b = c(1/2, -1, 2)
#
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x
#
err = -0.25 + 3*x - 4*x^2 - 2*x^5 + x^6
round0(err)
test.p(R, b)


### Example 2:
R = 1
b = c(1,-1,-2)
#
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x

test.p(R, b)


### Example 3:
R = 1
b = c(1, -3, 2)
#
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
err = -1 + 7*x - 12*x^2 + x^3 - x^5 + x^6
round0(err)


### Example 4:
R = 2
b = c(1,2,2)
#
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x

err = -1 - 10*x^2 + 8*x^3 - 2*x^5 + x^6
round0(err)
test.p(R, b)


### Example 5:
R = -108
b = c(2*9, 18, 2)
#
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x

err = -324 - 36*x^3 + 6*x^5 + x^6
round0(err)
test.p(R, b)


### Example 6:
R = -2
b = c(1, 2, 2)
# - but trivial solution;
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
y^3*x + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Examples:
R = 1
p = sapply(-6:6, function(r) print(solve.htx3y(R, c(1, r, 2))$p))
#
b = c(1, -4, 2)
sol = solve.htx3y(R, b)
x = sol$sol[,1]; y = sol$sol[,2];
sol
-1 + 9*x - 19*x^2 - x^5 + x^6

###
p = sapply(-6:6, function(r) print(solve.htx3y(r, c(r, -3, 2))$p))
#
r = -3; sol = solve.htx3y(r, c(r, -3, 2)); x = sol$sol[,1]
-9 - 21*x - 15*x^3 - x^5 + x^6


### Classic Polynomial:
b1 = b[1]; b2 = b[2]; b3 = b[3];
b1*(b3 - 1)*x^6 - R*(b3 - 1)*x^5 + b1*b2*(b3 - 2)*x^4 +
	+ (b1^2*(b3^2 + b3 - 2) + b2*R)*x^3 +
	- (b1*b2^2 + 2*b1*b3*R - b1*R)*x^2 +
	+ (R^2 - 2*b1^2*b2)*x - b1^3 # = 0


###################
### Simple variant:
### x^3*y + b*x

# x^3*y + b1*x = R
# y^3*x + b1*y = R

# trivial P4;

### Solution:

### Diff =>
# x*y = -b1/S

### Sum =>
# S^2 = -b1^2 / R;

solve.htx3y.S2P3 = function(R, b) {
	x.sum = sqrt(-b[1]^2 / R + 0i)
	x.sum = c(x.sum, -x.sum)
	xy = -b[1]/x.sum
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x=x, y=y)
	return(rbind(sol, sol[,2:1]))
}

### Example:
R = 1
b = 3
#
sol = solve.htx3y.S2P3(R, b)
x = sol[,1]; y = sol[,2];
sol


### Test
x^3*y + b[1]*x
y^3*x + b[1]*y


##########################
##########################
##########################

######################
### x^4*y^3 Series ###
######################

### Sections:
# 1. Chain: X
# 2. Chain: Y

### 1. Simple Side-Chain
### x^4*y^3 + b*x

### Variants:
### x^4*y^3 + b2*(x+y)^2 + b1*x

# x^4*y^3 + b2*(x+y)^2 + b1*x = R
# y^4*x^3 + b2*(x+y)^2 + b1*y = R

# trivial polynomial;

### Solution:

### Diff =>
# x*y = (-b1)^(1/3) * m3.all

### Sum =>
# S^2 = R / b2

m3.all = unity(3, all=T)

### Example
b = c(3, 1)
R = 1
#
xy = (-b[1] + 0i)^(1/3) * m3.all
x.sum = sqrt(R / b[2] + 0i)
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^4*y^3 + b[2]*(x+y)^2 + b[1]*x
y^4*x^3 + b[2]*(x+y)^2 + b[1]*y

### Classical Polynomial:
round0.p(poly.calc(sol[,1]))
# trivial
err = -3 - x^3 + 3*x^4 - 3*x^5 + x^6
round0(err)

# - trivial:
-b[1] + x^3*(x - sqrt(R/b[2]))^3
# - degenerate P[12] without sqrt;
b2^3*x^12 - 3*R*b2^2*x^10 + 3*R^2*b2*x^8 - (R^3 + 2*b1*b2^3)*x^6 - 6*b1*b2^2*R*x^4 + b1^2*b2^3 # = 0
# - although does NOT factorize;
(b2*x^4 - R*x^2)^3 - b1*b2^2*(2*b2*x^6 + 6*R*x^4 - b1*b2) # = 0


##################
### Extensions ###

### Extension 1:

### x^4*y^3 + b3*x*y + b2*x^2 + b1*x

# x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R
# y^4*x^3 + b3*x*y + b2*y^2 + b1*y = R

### Solution:

# "Trivial" solution:: x = y
# x^7 + (b2+b3)*x^2 + b1*x - R = 0
# &
# Trivial remaining P2;

### Diff =>
(x*y)^3*(x - y) + b2*(x-y)*(x+y) + b1*(x-y) # = 0
(x - y)*((x*y)^3 + b2*(x+y) + b1) # = 0
# Case: x != y
# (x*y)^3 = - (b2*S + b1);

### Sum =>
(x*y)^3*(x+y) + 2*b3*x*y + b2*(x^2 + y^2) + b1*(x+y) - 2*R # = 0
# - (b2*S + b1)*S + 2*b3*x*y + b2*(S^2 - 2*x*y) + b1*S - 2*R = 0
# 2*b3*x*y - 2*b2*x*y - 2*R = 0
# (b3-b2)*x*y = R
# x*y = R / (b3-b2)
# =>
# (x*y)^3 = R^3 / (b3-b2)^3
# b2*S + b1 + R^3 / (b3-b2)^3 = 0

### Example
b = c(3,2,1)
R = 1
#
xy = R / (b[3] - b[2])
x.sum = - (b[1] + R^3 / (b[3] - b[2])^3) / b[2]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
x = sol[,1]; y = sol[,2]
sol

### Test
x^4*y^3 + b[3]*x*y + b[2]*x^2 + b[1]*x
x^3*y^4 + b[3]*x*y + b[2]*y^2 + b[1]*y

# trivial P2 polynomial;

# all roots
x2 = roots(c(1, 0,0,0,0, (b[2]+b[3]), b[1], - R))
x = c(x, x2); y = c(y, x2)


###################
### 2. Chain: Y ###

### Extension 2
### x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y

# x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^4*x^3 + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
x^7 + b3*x^4 + b2*x^2 + b1*x - R # = 0

### Diff =>
(x*y)^3*(x - y) - b1*(x-y) # = 0
(x - y)*((x*y)^3 - b1) # = 0
# Case: x != y
# (x*y)^3 = b1;

### Sum =>
(x*y)^3*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) - 2*R # = 0
# ((x*y)^3 + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# 2*b1*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# b1*S = R - b3*(x*y)^2 - b2*x*y


### Solver:

solve.htMixt = function(b, R) {
	if(length(b) == 3) {
		xy = roots(c(1,0,0, -b[1]));
		r.sum = (R - b[3]*xy^2 - b[2]*xy) / b[1]
	} else {
		if(length(b) == 4) b = c(b, 0);
		xy = roots(c(1,0, b[5]-b[4], -b[1]))
		div = (b[4]*xy + b[1])
		r.sum = (R - b[3]*xy^2 - b[2]*xy) / div
	}
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	x = round0(x)
	sol = cbind(x, y)
	sol = sol[x != 0, ] # x == 0 => y -> Inf!
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}


### Example:
# b1 == 1 is trivial;
b = c(1, 2, 3)
R = 1
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1/2, 2, -3)
R = 1
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 3:
b = c(2,2,-2)
R = 2
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Classic Polynomial
round0.p(poly.calc(sol[,1]))

# TODO


###############
### Extension 3

### x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y

# x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^4*x^3 + b5*x*y^2 + b4*x^2*y + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
x^7 + b3*x^4 + (b4+b5)*x^3 + b2*x^2 + b1*x - R # = 0

### Diff =>
(x*y)^3*(x - y) + (b5-b4)*x*y*(x-y) - b1*(x-y) # = 0
(x - y)*((x*y)^3 + (b5-b4)*x*y - b1) # = 0
# Case: x != y
(x*y)^3 + (b5-b4)*x*y - b1 # = 0;
(x*y)^3 + db54*x*y - b1 # = 0;

### Diff:
# y*Eq 1 - x*Eq 2 =>
b4*x*y*(y^2 - x^2) + b3*(x*y)^2*(y - x) + b2*x*y*(y - x) + b1*(y^2 - x^2) - R*(y - x) # = 0
b4*x*y*S + b3*(x*y)^2 + b2*x*y + b1*S - R # = 0

### Sum =>
(x*y)^3*(x+y) + (b4+b5)*x*y*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) - 2*R # = 0
# ((x*y)^3 + (b4+b5)*x*y + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# (b4*x*y + b1)*S = R - b3*(x*y)^2 - b2*x*y


### Example:
b = c(1, 2, 3, -1, -2)
R = 1
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[5]*x^2*y + b[4]*x*y^2 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[5]*x*y^2 + b[4]*x^2*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1, 0, -1, -1, 1)
R = 1
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[5]*x^2*y + b[4]*x*y^2 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[5]*x*y^2 + b[4]*x^2*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x
#
err = 1 - 4*x - 2*x^2 + 2*x^3 - 2*x^5 + x^6
round0(err)


### Example 3:
b = c(1,1,1,0,-1)
R = 1
#
sol = solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol$p

### Test
err = 1 + 4*x + 6*x^2 - 4*x^4 - x^5 + x^6
round0(err)


### Classic Polynomial
b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4]; b5 = b[5]; db54 = b5 - b4;
b1*(b1^2 + db54*b4^2 + b4^3)*x^6 +
	- (3*b1^2*R + db54*b4^2*R + 2*b1^2*db54*b3 - 2*b1*db54*b2*b4 + 3*b1^2*b3*b4 - 3*b1*b2*b4^2)*x^5 +
	+ (3*b1*R^2 + b1*db54*b2^2 + 4*b1*db54*R*b3 - 3*b1^2*b2*b3 + b1*db54^2*b3^2 +
		- 2*db54*R*b2*b4 + 3*b1*b2^2*b4 + 3*b1*R*b3*b4 + b1*db54*b3^2*b4)*x^4 +
	+ (- R^3 + 2*b1^2*db54*b2 - db54*R*b2^2 + b1*b2^3 - 3*b1^3*b3 - 2*db54*R^2*b3 +
		+ 3*b1*R*b2*b3 - db54^2*R*b3^2 + b1*db54*b2*b3^2 + b1^2*b3^3 + 2*b1*db54*R*b4 +
		+ 3*b1^2*b2*b4 + 2*b1*db54^2*b3*b4 + 3*b1*R*b4^2 + 2*b1*db54*b3*b4^2)*x^3 +
	+ (b1^3*db54 - 2*b1*db54*R*b2 + 3*b1^2*b2^2 + 3*b1^2*R*b3 + b1^2*db54*b3^2 - 2*db54*R^2*b4 +
		+ 3*b1*R*b2*b4 - 2*db54^2*R*b3*b4 + 2*b1*db54*b2*b3*b4 + 3*b1^2*b3^2*b4 +
		+ b1*db54^2*b4^2 + b1*db54*b4^3)*x^2 +
	+ (3*b1^2*b4*R - b1^2*db54*R - db54^2*b4^2*R + 3*b1^3*b2 + 2*b1^2*db54*b3*b4 +
		+ b1*db54*b2*b4^2 + 3*b1^2*b3*b4^2)*x +
	+ b1^2*(b1^2 + db54*b4^2 + b4^3) # = 0

# P[8] => P[6]
pR = div.pm(pR, toPoly.pm("(b3*x + b4)^2"), "x")


##########################
##########################

######################
### x^4*y^2 Series ###
######################

### Simple Side-Chain
### x^4*y^2 + b*y

### Solution:
# - Case: distinct roots;

### Diff =>
(x*y)^2*S - b # = 0

### Diff:
# y^2*Eq 1 - x^2*Eq 2 =>
b*(y^3 - x^3) - R*(y^2 - x^2) # = 0
b*(S^2 - x*y) - R*S # = 0

### Eq S:
b^2*S^5 - 2*R*b*S^4 + R^2*S^3 - b^3 # = 0


### Solver:

solve.S2Ht.L42ChY = function(R, b, debug=TRUE, all=TRUE) {
	coeff = c(b^2, - 2*R*b, R^2, 0, 0, - b^3);
	S = roots(coeff);
	if(debug) print(S);
	xy = (b*S^2 - R*S) / b;
	len = length(S);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], xy[id]));
	})
	x12 = t(x12);
	x = x12[,1]; y = x12[,2];
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.L42ChY = function(sol, b, R = NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^4*y^2 + b*y;
	err2 = y^4*x^2 + b*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = -3
sol = solve.S2Ht.L42ChY(R, b)

test.S2Ht.L42ChY(sol, b=b)


### Ex 2:
R = 5
b = -3
sol = solve.S2Ht.L42ChY(R, b)

test.S2Ht.L42ChY(sol, b=b)


### Classic Polynomial:
x = sol[,1];
b^2*x^10 - 2*b*R*x^9 + R^2*x^8 - b^3*x^5 + 3*b^2*R*x^4 - b*R^2*x^3 - R^3*x^2 + b^4 # = 0


##############
### Extension:

### Side-Chain: y^2
### x^4*y^2 + b2*y^2 + b1*y

### Solution:
# - Case: distinct roots;

### Diff =>
(x*y)^2*S - b2*S - b1 # = 0

### Diff:
# y^2*Eq 1 - x^2*Eq 2 =>
b2*(y^4 - x^4) + b1*(y^3 - x^3) - R*(y^2 - x^2) # = 0
b2*(S^3 - 2*x*y*S) + b1*(S^2 - x*y) - R*S # = 0

### Eq S:
b2^2*S^7 + 2*b1*b2*S^6 + (b1^2 - 2*R*b2)*S^5 - 2*R*b1*S^4 + (R^2 - 4*b2^3)*S^3 +
	- 8*b1*b2^2*S^2 - 5*b1^2*b2*S - b1^3 # = 0


### Solver:

solve.S2Ht.L42ChY2 = function(R, b, debug=TRUE, all=TRUE) {
	b1 = b[1]; b2 = b[2];
	coeff = c(b2^2, 2*b1*b2, b1^2 - 2*b2*R, - 2*b1*R, R^2 - 4*b2^3,
		- 8*b1*b2^2, - 5*b1^2*b2, - b1^3);
	S = roots(coeff);
	if(debug) print(S);
	xy = (b2*S^3 + b1*S^2 - R*S) / (2*b2*S + b1);
	len = length(S);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], xy[id]));
	})
	x12 = t(x12);
	x = x12[,1]; y = x12[,2];
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.L42ChY2 = function(sol, b, R = NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^4*y^2 + b[2]*y^2 + b[1]*y;
	err2 = y^4*x^2 + b[2]*x^2 + b[1]*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = c(-1, 3)
sol = solve.S2Ht.L42ChY2(R, b)

test.S2Ht.L42ChY2(sol, b=b)


### Ex 2:
R = 5
b = c(-3, -1)
sol = solve.S2Ht.L42ChY2(R, b)

test.S2Ht.L42ChY2(sol, b=b)


### Ex 3:
R = 2
b = c(-1, -1)
sol = solve.S2Ht.L42ChY2(R, b)

test.S2Ht.L42ChY2(sol, b=b)


### Classic Polynomial:
x = sol[,1]; b1 = b[1]; b2 = b[2];
b2^2*x^14 + 2*b1*b2*x^13 + (b1^2 - 2*b2*R)*x^12 - 2*b1*R*x^11 + (3*b2^3 + R^2)*x^10 +
	+ 5*b1*b2^2*x^9 + (b1^2*b2 - 5*b2^2*R)*x^8 - b1*(b1^2 + 2*b2*R)*x^7 +
	+ 3*(b2^4 + b1^2*R + b2*R^2)*x^6 + b1*(4*b2^3 - R^2)*x^5 +
	- (4*b2^3*R + R^3)*x^4 + (b1^4 + b2^5 + 4*b1^2*b2*R + 2*b2^2*R^2)*x^2 +
	+ b1*b2^4*x - b2^4*R # = 0


##########################
##########################

####################
### x^4*y Series ###
####################

# x^4*y + b1*x = R
# y^4*x + b1*y = R

### Solution:

### Diff =>
x*y*(x^3 - y^3) + b1*(x-y) # = 0
(x - y)*(x*y*(x^2 + y^2 + x*y) + b1) # = 0
(x - y)*(x*y*(x+y)^2 - (x*y)^2 + b1) # = 0
# Case: x != y
x*y*(x+y)^2 - (x*y)^2 + b1 # = 0
# (x*y)^2 - x*y*S^2 = b1

### Sum =>
x*y*(x^3 + y^3) + b1*(x+y) - 2*R # = 0
(x+y)*x*y*(x^2 + y^2 - x*y) + b1*(x+y) - 2*R # = 0
S*x*y*(S^2 - 3*x*y) + b1*S - 2*R # = 0
S^3*x*y - 3*S*(x*y)^2 + b1*S - 2*R
S^3*x*y - 3*S*(x*y*S^2 + b1) + b1*S - 2*R
-2*S^3*x*y - 2*b1*S - 2*R
S^3*x*y + b1*S + R
# x*y = -(b1*S + R) / S^3;
# =>
(b1*S + R)^2 / S^6 + (b1*S + R) / S^3 * S^2 - b1 # = 0
(b1*S + R)^2 + (b1*S + R)*S^5 - b1*S^6 # = 0
(b1*S + R)^2 + R*S^5
R*S^5 + b1^2*S^2 + 2*b1*R*S + R^2 # = 0


### Solver:

solve.S2Ht.L41ChX = function(R, b, debug=TRUE, all=TRUE) {
	S = roots(c(R, 0, 0, b[1]^2, 2*b[1]*R, R^2));
	if(debug) print(S);
	#
	xy = -(b[1]*S + R) / S^3;
	x.diff = sqrt(S^2 - 4*xy + 0i);
	x = (S + x.diff)/2;
	y = (S - x.diff)/2;
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[,2:1]);
	return(sol);
}
test.S2Ht.L41ChX = function(sol, b, R=NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^4*y + b[1]*x;
	err2 = y^4*x + b[1]*y;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
b = 3
R = 1
#
sol = solve.S2Ht.L41ChX(R, b)

test.S2Ht.L41ChX(sol, b=b)


### Ex 2:
b = 3
R = -2
#
sol = solve.S2Ht.L41ChX(R, b)

test.S2Ht.L41ChX(sol, b=b)


### Test
x = sol[,1]; y = sol[,2];
x^4*y + b[1]*x
y^4*x + b[1]*y

### Classical Polynomial: P10
round0.p(poly.calc(sol[,1]))

b = b[1];
err = R*x^10 + b^2*x^7 - 2*R*b*x^6 + R^2*x^5 - b^3*x^3 + 3*R*b^2*x^2 - 3*R^2*b*x + R^3
round0(err)


### Derivation:
# y = (R - b[1]*x)/x^4
(R - b[1]*x)^4/x^16 * x + b[1]*(R - b[1]*x)/x^4 - R # = 0
(R - b[1]*x)^4/x^15 + b[1]*(R - b[1]*x)/x^4 - R
(R - b[1]*x)^4 + b[1]*(R - b[1]*x)*x^11 - R*x^15
R*x^15 - b[1]*(R - b[1]*x)*x^11 - (R - b[1]*x)^4
R*x^15 + b[1]^2*x^12 - b[1]*R*x^11 - (R - b[1]*x)^4

R*x^15 + b[1]^2*x^12 - R*b[1]*x^11 - b[1]^4*x^4 + 4*R*b[1]^3*x^3 - 6*R^2*b[1]^2*x^2 + 4*R^3*b[1]*x^1 - R^4
(x^5 + b[1]*x - R)*(R*x^10 + b[1]^2*x^7 - 2*R*b[1]*x^6 + R^2*x^5 - b[1]^3*x^3 + 3*R*b[1]^2*x^2 - 3*R^2*b[1]*x + R^3)


##########################
##########################
##########################

######################
### x^5*y^3 Series ###
######################

### Simple Side-Chain
### x^5*y^3 + b*y

### Solution:
# - Case: distinct roots;

### Diff =>
(x*y)^3*S - b # = 0

### Diff:
# y^2*Eq 1 - x^2*Eq 2 =>
b*(y^3 - x^3) - R*(y^2 - x^2) # = 0
b*(S^2 - x*y) - R*S # = 0
b*x*y - b*S^2 + R*S # = 0

### Eq S:
b^3*S^7 - 3*R*b^2*S^6 + 3*R^2*b*S^5 - R^3*S^4 - b^4 # = 0


### Solver:

solve.S2Ht.L53ChY = function(R, b, debug=TRUE, all=TRUE) {
	coeff = c(b^3, - 3*R*b^2, 3*R^2*b, - R^3,
		0, 0, 0, - b^4);
	S = roots(coeff);
	if(debug) print(S);
	xy = (b*S^2 - R*S) / b;
	len = length(S);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], xy[id]));
	})
	x12 = t(x12);
	x = x12[,1]; y = x12[,2];
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.L53ChY = function(sol, b, R = NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^5*y^3 + b*y;
	err2 = y^5*x^3 + b*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = -3
sol = solve.S2Ht.L53ChY(R, b)

test.S2Ht.L53ChY(sol, b=b)


### Ex 2:
R = 5
b = -3
sol = solve.S2Ht.L53ChY(R, b)

test.S2Ht.L53ChY(sol, b=b)


### Classic Polynomial:
b^3*x^14 - 3*R*b^2*x^13 + 3*R^2*b*x^12 - R^3*x^11 - b^4*x^7 - R*b^3*x^6 +
	+ 4*R^2*b^2*x^5 - R^3*b*x^4 - R^4*x^3 + b^5 # = 0


##############
### Extension:

### Side-Chain: y^2
### x^5*y^3 + b2*y^2 + b1*y

### Solution:
# - Case: distinct roots;

### Diff =>
(x*y)^3*S - b2*S - b1 # = 0

### Diff:
# y^2*Eq 1 - x^2*Eq 2 =>
b2*(y^4 - x^4) + b1*(y^3 - x^3) - R*(y^2 - x^2) # = 0
b2*(S^2 - 2*x*y)*S + b1*(S^2 - x*y) - R*S # = 0
(2*b2*S + b1)*x*y - b2*S^3 - b1*S^2 + R*S # = 0

### Eq S:
b2^3*S^10 + 3*b1*b2^2*S^9 + (3*b1^2*b2 - 3*b2^2*R)*S^8 + (b1^3 - 6*b1*b2*R)*S^7 +
	- 3*b1^2*R*S^6 + 3*b2*R^2*S^6 + 3*b1*R^2*S^5 - R^3*S^4 - 8*b2^4*S^4 - 20*b1*b2^3*S^3 +
	- 18*b1^2*b2^2*S^2 - 7*b1^3*b2*S - b1^4 # = 0


### Solver:

solve.S2Ht.L53ChY2 = function(R, b, debug=TRUE, all=TRUE) {
	b1 = b[1]; b2 = b[2];
	coeff = c(b2^3, 3*b1*b2^2, 3*b1^2*b2 - 3*b2^2*R, b1^3 - 6*b1*b2*R,
		- 3*b1^2*R + 3*b2*R^2, 3*b1*R^2, - (R^3 + 8*b2^4), - 20*b1*b2^3,
		- 18*b1^2*b2^2, - 7*b1^3*b2, - b1^4);
	S = roots(coeff);
	if(debug) print(S);
	xy  = (b2*S^3 + b1*S^2 - R*S) / (2*b2*S + b1);
	len = length(S);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], xy[id]));
	})
	x12 = t(x12);
	x = x12[,1]; y = x12[,2];
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.L53ChY2 = function(sol, b, R = NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^5*y^3 + b[2]*y^2 + b[1]*y;
	err2 = y^5*x^3 + b[2]*x^2 + b[1]*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = c(-1, -3)
sol = solve.S2Ht.L53ChY2(R, b)

test.S2Ht.L53ChY2(sol, b=b)

poly.calc(sol[,1])


### Ex 2:
R = 5
b = c(-3, 1)
sol = solve.S2Ht.L53ChY2(R, b)

test.S2Ht.L53ChY2(sol, b=b)


### Ex 3:
R = 4
b = c(-2, 1)
sol = solve.S2Ht.L53ChY2(R, b)

test.S2Ht.L53ChY2(sol, b=b)


### Classic Polynomial:
x = sol[,1]; b1 = b[1]; b2 = b[2];
b2^3*x^20 + 3*b1*b2^2*x^19 + 3*b2*(b1^2 - b2*R)*x^18 + (b1^3 - 6*b1*b2*R)*x^17 +
	- 3*(b1^2*R - b2*R^2)*x^16 + 3*b1*R^2*x^15 - (b2^4 + R^3)*x^14 - 4*b1*b2^3*x^13 +
	- b2^2*(6*b1^2 - 4*b2*R)*x^12 - b1*b2*(4*b1^2 - 7*b2*R)*x^11 +
	- (b1^4 - 2*b1^2*b2*R + 6*b2^2*R^2)*x^10 - (b1^3*R + 2*b1*b2*R^2)*x^9 +
	+ (b2^5 + 4*b1^2*R^2 + 4*b2*R^3)*x^8 + b1*(5*b2^4 - R^3)*x^7 +
	+ (5*b1^2*b2^3 - R^4)*x^6 - 5*b1*b2^3*R*x^5 + (b1^5 + 5*b1^3*b2*R + 5*b1*b2^2*R^2)*x^3 +
	- b2^6*x^2 - b1*b2^5*x + b2^5*R # = 0

# P[24] => P[20]
pR = div.pm(pR, "(b2^2*x^2 - b1^2 - b2*R)^2", "x")


##########################
##########################

######################
### x^5*y^2 Series ###
######################

### Simple Side-Chain
### x^5*y^2 + b*y

### Solution:
# - Case: distinct roots;

### Diff =>
(x*y)^2*(S^2 - x*y) - b # = 0

### Diff:
# y^3*Eq 1 - x^3*Eq 2 =>
b*(y^4 - x^4) - R*(y^3 - x^3) # = 0
b*(S^3 - 2*x*y*S) - R*(S^2 - x*y) # = 0

### Eq S:
b^2*S^9 - 2*R*b*S^8 + R^2*S^7 - 8*b^3*S^3 + 12*R*b^2*S^2 - 6*R^2*b*S + R^3 # = 0


### Solver:

solve.S2Ht.L52ChY = function(R, b, debug=TRUE, all=TRUE) {
	coeff = c(b^2, - 2*R*b, R^2, 0, 0, 0,
		- 8*b^3, 12*R*b^2, - 6*R^2*b, R^3);
	S = roots(coeff);
	if(debug) print(S);
	xy = (b*S^3 - R*S^2) / (2*b*S - R);
	len = length(S);
	x12 = sapply(seq(len), function(id) {
		roots(c(1, -S[id], xy[id]));
	})
	x12 = t(x12);
	x = x12[,1]; y = x12[,2];
	sol = cbind(x, y);
	if(all) sol = rbind(sol, sol[, c(2,1)]);
	return(sol);
}
test.S2Ht.L52ChY = function(sol, b, R = NULL) {
	x = sol[,1]; y = sol[,2];
	err1 = x^5*y^2 + b*y;
	err2 = y^5*x^2 + b*x;
	err = rbind(err1, err2);
	err = round0(err);
	return(err);
}

### Examples:

### Ex 1:
R = 2
b = -3
sol = solve.S2Ht.L52ChY(R, b)

test.S2Ht.L52ChY(sol, b=b)


### Ex 2:
R = 5
b = -3
sol = solve.S2Ht.L52ChY(R, b)

test.S2Ht.L52ChY(sol, b=b)


### Classic Polynomial:
b^2*x^18 - 2*b*R*x^17 + R^2*x^16 - b^3*x^12 + 3*b^2*R*x^11 - 3*b*R^2*x^10 + R^3*x^9 +
	+ b^4*x^6 - 4*b^3*R*x^5 + b^2*R^2*x^4 + b*R^3*x^3 + R^4*x^2 - b^5 # = 0


