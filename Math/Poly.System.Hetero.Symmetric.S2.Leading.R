########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### with Mixed Leading Term
###
### draft v.0.1a


### Heterogeneous Symmetric Polynomial Systems

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
# M31.2) x^3*y + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6; some nice)
# M41.1) x^4*y + b*x; (P5 => P10)
### Mixed: Order n+3:
# M43.1) x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R; (trivial P2; base P7)
# M43.2) x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)
# M43.3) x^4*y^3 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*y = R; (TODO: P3 => P6)


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

library(polynom)
library(pracma)

### helper Functions

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


solve.htx3y = function(R, b) {
	if(R == 0 && b[1] == 0) return(list(sol=NA, p=NA))
	r.sum = roots(c(b[1], - R, b[1]*b[2], b[1]^2*b[3] - b[1]^2))
	r.sum = round0(r.sum)
	r.sum = r.sum[ r.sum != 0 ] # avoid division by 0
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
	err = b1^3*(1 - b3^2) +
	(1 - b3^2)*(- R^2 + 2*b1^2*b2)*x +
	b1*(1 - b3^2)*(- R + 2*R*b3 + b2^2)*x^2 +
	(1 - b3^2)*(- R*b2 + 2*b1^2 - b1^2*b3 - b1^2*b3^2)*x^3 +
	b1*b2*(b3^2 - 1)*(b3 - 2)*x^4 +
	- R*(b3 - 1)*(b3^2 - 1)*x^5 +
	b1*(b3 - 1)*(b3^2 - 1)*x^6
	return(round0(err))
}

### Classic Polynomial:
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

######################
### x^4*y^3 Series ###
######################

### Simple Side-Chain
### x^4*y^3 + b*x

### Variants:
### x^4*y^3 + b2*(x+y)^2 + b1*x

# x^4*y^3 + b2*(x+y)^2 + b1*x = R
# y^3*x^4 + b2*(x+y)^2 + b1*y = R

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

# trivial:
-b[1] + x^3*(x - sqrt(R/b[2]))^3


##################
### Extensions ###

### Extension 1:

### x^4*y^3 + b3*x*y + b2*x^2 + b1*x

# x^4*y^3 + b3*x*y + b2*x^2 + b1*x = R
# y^3*x^4 + b3*x*y + b2*y^2 + b1*y = R

### Solution:

# "Trivial" solution:: x = y
# x^7 + (b2+b3)*x^2 + b1*x - R = 0
# &
# Trivial remaining P2;

### Diff =>
# (x*y)^3*(x - y) + b2*(x-y)*(x+y) + b1*(x-y) = 0
# (x - y)*((x*y)^3 + b2*(x+y) + b1) = 0
# Case: x != y
# (x*y)^3 = - (b2*S + b1);

### Sum =>
# (x*y)^3*(x+y) + 2*b3*x*y + b2*(x^2 + y^2) + b1*(x+y) = 2*R
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


###############
### Extension 2

### x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y

# x^4*y^3 + b3*(x*y)^2 + b2*x*y + b1*y = R
# y^3*x^4 + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
# x^7 + b3*x^4 + b2*x^2 + b1*x - R = 0

### Diff =>
# (x*y)^3*(x - y) - b1*(x-y) = 0
# (x - y)*((x*y)^3 - b1) = 0
# Case: x != y
# (x*y)^3 = b1;

### Sum =>
# (x*y)^3*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) = 2*R
# ((x*y)^3 + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# 2*b1*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# b1*S = R - b3*(x*y)^2 - b2*x*y

solve.htMixt = function(b, R) {
	if(length(b) == 3) {
		xy = roots(c(1,0,0, -b[1]))
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
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1/2, 2, -3)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 3:
b = c(2,2,-2)
R = 2
#
sol =solve.htMixt(b, R)
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
# y^3*x^4 + b5*x^2*y + b4*x*y^2 + b3*(x*y)^2 + b2*x*y + b1*x = R

### Solution:

# "Trivial" solution: x = y
# x^7 + b3*x^4 + (b4+b5)*x^3 + b2*x^2 + b1*x - R = 0

### Diff =>
# (x*y)^3*(x - y) + (b5-b4)*x*y*(x-y) - b1*(x-y) = 0
# (x - y)*((x*y)^3 + (b5-b4)*x*y - b1) = 0
# Case: x != y
# (x*y)^3 + (b5-b4)*x*y - b1 = 0;

### Sum =>
# (x*y)^3*(x+y) + (b4+b5)*x*y*(x+y) + 2*b3*(x*y)^2 + 2*b2*x*y + b1*(x+y) = 2*R
# ((x*y)^3 + (b4+b5)*x*y + b1)*S = 2*R - 2*b3*(x*y)^2 - 2*b2*x*y
# (b4*x*y + b1)*S = R - b3*(x*y)^2 - b2*x*y


### Example:
b = c(1, 2, 3, -1, -2)
R = 1
#
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
x^4*y^3 + b[5]*x^2*y + b[4]*x*y^2 + b[3]*(x*y)^2 + b[2]*x*y + b[1]*y
x^3*y^4 + b[5]*x*y^2 + b[4]*x^2*y + b[3]*(x*y)^2 + b[2]*x*y + b[1]*x


### Example 2:
b = c(1, 0, -1, -1, 1)
R = 1
#
sol =solve.htMixt(b, R)
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
sol =solve.htMixt(b, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
err = 1 + 4*x + 6*x^2 - 4*x^4 - x^5 + x^6
round0(err)


### Classic Polynomial
# TODO

