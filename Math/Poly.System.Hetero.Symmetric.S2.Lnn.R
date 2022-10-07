########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### Mixed/Trivial Leading Term
###
### draft v.0.1a

### Polynomial Systems
### Heterogeneous Symmetric
### Leading Terms: Mixed/Trivial


### Types:

### Mixed-Trivial: (x*y)^n
### Order 2:
# Mt2.1.) x^2*y^2 + b2*x^2 + b1*x: simple (P2 => P4);
# Mt2.2.) x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x; (TODO: P2 => P4)
# Mt2.3.) TODO: x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x
### Mixed: (x*y)^n
# Mt3.1.) TODO: (x*y)^3 + ...:
#      (x*y)^3 + b*x^3 = R;
#      (x*y)^3 + b*x^4 + b2*(x*y)^2 + b1*x*y = R;
# Mt5.1) TODO: (x*y)^5 + b*x^3 = R


###############

###############
### History ###
###############


### v.0.1a:
# - moved sections with Mixed/Trivial Leading term:
#   from Poly.System.Hetero.Symmetric.R
#   to this file;

### Note:
### Non-Trivial Systems: (x^n*y^m)
# - are in file:
#   Poly.System.Hetero.Symmetric.S2.Leading.R;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


##########################
##########################

##########################
### Mixed Leading Term ###
### Trivial: (x*y)^n   ###
##########################

### Trivial Systems: (x*y)^n
# - show similar behavior to that of the remaining side chain;


####################
### x^2*y^2 Term ###
####################

### x^2*y^2 + b2*x^2 + b1*x

# x^2*y^2 + b2*x^2 + b1*x = R
# y^2*x^2 + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# S = -b1/b2

### Sum =>
(x*y)^2 - b2*x*y - R # = 0


### Example:
b = c(3, 1)
R = 1
#
x.sum = - b[1]/b[2]
xy = roots(c(1, - b[2], - R))
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol


### Test
x^2*y^2 + b[2]*x^2 + b[1]*x
y^2*x^2 + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

- R + b[1]*x + (b[2] + b[1]^2/b[2]^2)*x^2 + 2*b[1]/b[2]*x^3 + x^4


######################################


######################################
### x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x

# x^2*y^2 + b3*x^2*y + b2*x^2 + b1*x = R
# y^2*x^2 + b3*y^2*x + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# b3*x*y*(x-y) + b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*(b3*x*y + b2*(x+y) + b1) = 0
# Case: x != y
# b3*x*y = - (b2*Z + b1)
# Z = -(b3*x*y + b1) / b2

### Sum =>
# 2*(x*y)^2 + b3*x*y*(x+y) + b2*(x^2 + y^2) + b1*(x+y) = 2*R
# 2*(x*y)^2 - (b2*Z + b1)*Z + b2*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# 2*(x*y)^2 - 2*b2*x*y - 2*R
# (x*y)^2 - b2*x*y - R

### Example:
b = c(3, 1, 1)
R = 1
#
xy = roots(c(1, -b[2], -R))
x.sum = -(b[3] * xy + b[1]) / b[2]
x.diff = sqrt(x.sum^2 - 4*xy + 0i)
x = (x.sum + x.diff)/2
y = (x.sum - x.diff)/2
sol = cbind(x, y)
sol = rbind(sol, sol[,2:1])
sol

### Test
x^2*y^2 + b[3]*x^2*y + b[2]*x^2 + b[1]*x
y^2*x^2 + b[3]*y^2*x + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))

### TODO: classical Polynomial: P4;


#################################


#################################
### x^2*y^2 + b4*x^2*y + b3*x*y^2
###  + b2*x^2 + b1*x

# x^2*y^2 + b4*x^2*y + b3*x*y^2 + b2*x^2 + b1*x = R
# y^2*x^2 + b4*y^2*x + b3*y*x^2 + b2*y^2 + b1*y = R

### Solution:

### Diff =>
# (b4 - b3)*x*y*(x-y) + b2*(x^2 - y^2) + b1*(x-y) = 0
# (x - y)*((b4 - b3)*x*y + b2*(x+y) + b1) = 0
# Case: x != y

# Sub-Case 1: b4 = b3
# x + y = - b1/b2;
# Z = - b1/b2;
#
# Sub-Case 2: b4 != b3
# (b4 - b3)*x*y = - (b2*Z + b1)
# Z = -((b4 - b3)*x*y + b1) / b2

### Sum =>
# 2*(x*y)^2 + (b3+b4)*x*y*(x+y) + b2*(x^2 + y^2) + b1*(x+y) = 2*R
#  2*(x*y)^2 + (b3+b4)*x*y*Z + b2*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# Sub-Case 1:
#  2*(x*y)^2 - b1/b2*(b3+b4)*x*y + b2*(b1^2/b2^2 - 2*x*y) - b1*b1/b2 - 2*R = 0
#  2*(x*y)^2 - b1/b2*(b3+b4)*x*y - 2*b2*x*y - 2*R = 0
#  2*b2*(x*y)^2 - (b1*(b3+b4) + 2*b2^2)*x*y - 2*b2*R = 0
# Sub-Case 2:
# ...

### TODO: Sub-Case 2;

solve.htMixt22 = function(b, R) {
	if(b[3] == b[4]) {
		x.sum = - b[1]/b[2];
		xy = roots(c(2*b[2], - (b[1]*(b[3]+b[4]) + 2*b[2]^2), - 2*b[2]*R))
	}
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	sol
}

### Example SC-1:
b = c(3,2,1,1)
R = 1
#
sol = solve.htMixt22(b, R)
x = sol[,1]; y = sol[,2];
sol

### Test
x^2*y^2 + b[4]*x^2*y + b[3]*x*y^2 + b[2]*x^2 + b[1]*x
y^2*x^2 + b[4]*y^2*x + b[3]*y*x^2 + b[2]*y^2 + b[1]*y

### Classical Polynomial
round0.p(poly.calc(sol[,1]))


### Example SC-2:
### TODO


##########################
##########################

######################
### (x*y)^3 Series ###
######################

### (x*y)^3 + b*x^3

# (x*y)^3 + b1*x^3 = R
# (x*y)^3 + b1*y^3 = R


m3 = unity(3, all=F)

### Solution:

### Diff =>
x^3 - y^3 # = 0
# y = x * m3

### =>
# x^6 + b1*x^3 = R

### Example
b = 1
R = 1
#
x = roots(c(1, b[1], - R))
x = as.vector(sapply(x, function(r) roots(c(1,0,0, -r))))
y = x*m3
sol = cbind(x, y)
# includes the Equal roots;
sol = rbind(sol, sol[,2:1])
sol

### Test
(x*y)^3 + b[1]*x^3
(x*y)^3 + b[1]*y^3

### Classic Polynomial:
round0.p(poly.calc(sol[1:6, 1]))
# P6 (using only 1 set of roots);

x^6 + b1*x^3 - R # = 0

