########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### Multiple/Mixed Leading Terms
###
### draft v.0.1a

### Polynomial Systems
### Heterogeneous Symmetric
### Leading Terms: Multiple + Mixed


### Types:

### Mixed: other
# Mm3.1.) x^3 + a*x^3*y = R; (TODO: P3 => P6)
# Mm3.2.) x^3 + a*x^3*y + b1*x*y = R;
# Mm3.3.) x^3 + a*x*y^3 = R; (TODO: P5 => P10)
#         https://www.youtube.com/watch?v=jRwiJTaAu4c
### Mixed 2 High-Power Terms:
# MT.1) a1*x^3*y + a2*x*y^3 + b*x = R; [P3 => interesting P6]


###############

###############
### History ###
###############


### v.0.1a:
# - started to move sections with Multiple/Mixed Leading terms:
#   from Poly.System.Hetero.Symmetric.R
#   to this file;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


##########################
##########################

##########################
### 2 Leading Terms:   ###
###   Simple + Mixed   ###
##########################

###############
### x^3 + x^3*y

# x^3 + a*x^3*y = R
# y^3 + a*y^3*x = R

### Extension:
# x^3 + a*x^3*y + b1*x*y = R

### Solution:

### Diff =>
# (x-y)*(x^2 + y^2 + x*y + a*x*y*(x+y)) = 0
# Case: x != 0 =>
# x^2 + y^2 + x*y + a*x*y*(x+y) = 0
# S^2 - x*y + a*x*y*S = 0
# x*y = - S^2 / (a*S - 1)

### Sum =>
# x^3 + y^3 + a*x*y*(x^2 + y^2) - 2*R = 0
# S^3 - 3*x*y*S + a*x*y*(S^2 - 2*x*y) - 2*R
# S^3 + 3*S^3 / (a*S - 1) - a*S^2*(S^2 - 2*x*y) / (a*S - 1) - 2*R
# a*S^4 - S^3 + 3*S^3 - a*S^2*(S^2 - 2*x*y) - 2*R*(a*S - 1)
# 2*S^3 + 2*a*S^2*x*y - 2*R*a*S + 2*R
# 2*S^3 - 2*a*S^4 / (a*S - 1) - 2*R*a*S + 2*R
# 2*a*S^4 - 2*S^3 - 2*a*S^4 - 2*R*a^2*S^2 + 2*R*a*S + 2*a*R*S - 2*R
# - 2*S^3 - 2*R*a^2*S^2 + 4*a*R*S - 2*R
# S^3 + R*a^2*S^2 - 2*a*R*S + R

### Extension: + b1*x*y
# (b1*a+1)*S^3 + (a^2*R-b1)*S^2 - 2*a*R*S + R

solve.x3y = function(a, R, b=0, type=1) {
	# type = 1: simple;
	if(type == 2) {
		coeff = c(b[1]*a[1]+1, a[1]^2*R - b[1], - 2*a[1]*R, R)
	} else if(type == 1) {
		coeff = c(1, R*a[1]^2, - 2*a[1]*R, R)
	} else if(type == 5) {
		coeff = c(a[1]^2, - a[1], -1, - a[1]^2*R, - 2*a[1]*R, - R)
	}
	x.s = roots(coeff)
	xy = if(type < 5) x.s^2 / (1 - a[1]*x.s) else x.s^2 / (a[1]*x.s + 1);
	x.d = sqrt(x.s^2 - 4*xy + 0i)
	x = (x.s + x.d)/2
	y = (x.s - x.d)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example
a = 3
R = 1
#
sol = solve.x3y(a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x^3*y
y^3 + a[1]*y^3*x

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))
err = 1 + 3*x - 2*x^3 - 3*x^4 + 9*x^5 + x^6
round0(err)

###
R = 1
l = sapply(-6:6, function(a) print(solve.x3y(a, R)$p))

### Extension

### Example 2: Extension
a = -3
b = c(1/2)
R = 2
#
sol = solve.x3y(a, R=R, b=b, type=2)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x^3*y + b[1]*x*y
y^3 + a[1]*y^3*x + b[1]*x*y



################
### Alternative:
### x^3 + x*y^3

# x^3 + a*x*y^3 = R
# y^3 + a*x^3*y = R

# Cheol-Hyun Cho:
# Fukaya category for Landau-Ginzburg orbifolds and Berglund-Hübsch
# homological mirror symmetry for curve singularities
# [from 53:30] https://www.youtube.com/watch?v=jRwiJTaAu4c

### TODO: add x*y & (x*y)^2 terms;

### Diff =>
# (x-y)*(x^2 + y^2 + x*y - a*x*y*(x+y)) = 0
# Case: x != 0 =>
# x^2 + y^2 + x*y - a*x*y*(x+y) = 0
# S^2 - x*y - a*x*y*S = 0
# x*y = S^2 / (a*S + 1)

### Sum =>
# x^3 + y^3 + a*x*y*(x^2 + y^2) - 2*R = 0
# S^3 - 3*x*y*S + a*x*y*(S^2 - 2*x*y) - 2*R
# S^3 - 3*S^3 / (a*S + 1) + a*S^2*(S^2 - 2*x*y) / (a*S + 1) - 2*R
# a*S^4 + S^3 - 3*S^3 + a*S^2*(S^2 - 2*x*y) - 2*R*(a*S + 1)
# a*S^4 - S^3 - a*S^2*x*y - a*R*S - R
# a*S^4 - S^3 - a*S^4/(a*S + 1) - a*R*S - R
# a^2*S^5 - a*S^4 - S^3 - a^2*R*S^2 - 2*a*R*S - R

### Example
a = 1
R = 1
#
sol = solve.x3y(a, R, type=5)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
x^3 + a[1]*x*y^3
y^3 + a[1]*y*x^3

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))
err = 1 - 3*x + 3*x^2 - 3*x^3 + 4*x^4 - 3*x^5 + 3*x^6 - 2*x^7 + x^8 - x^9 + x^10
round0(err)


############################
############################

############################
### 2 Mixed Highest Term ###
############################

###############################
### a1*x^j*y^k + a2*x^k*y^j ###
###############################

### a1*x^3*y + a2*x*y^3 + b*x

# a1*x^3*y + a2*x*y^3 + b1*x = R
# a1*x*y^3 + a2*x^3*y + b1*y = R

### Solution

### Diff =>
(a1-a2)*x*y*(x^2 - y^2) + b1*(x-y) # = 0
(x-y)*((a1-a2)*x*y*S + b1)#  = 0
# Case: x != y
(a1-a2)*x*y*S + b1 # = 0
# x*y*Z = -b1/(a1-a2);

### Sum =>
# (a1+a2)*x*y*(x^2 + y^2) + b1*(x+y) = 2*R
# (a1+a2)*x*y*(Z^2 - 2*x*y) + b1*Z - 2*R = 0
# (a1+a2)*x*y*Z*(Z^3 - 2*x*y*Z) + b1*Z^3 - 2*R*Z^2 = 0
# =>
# -(a1+a2)*b1/(a1-a2)*(Z^3 + 2*b1/(a1-a2)) + b1*Z^3 - 2*R*Z^2
# -(a1+a2)*b1*((a1-a2)*Z^3 + 2*b1) + b1*(a1-a2)^2*Z^3 - 2*(a1-a2)^2*R*Z^2
# -b1*(a1+a2)*(a1-a2)*Z^3 + b1*(a1-a2)^2*Z^3 - 2*(a1-a2)^2*R*Z^2 - 2*b1^2*(a1+a2)
# 2*b1*a2*(a1-a2)*Z^3 + 2*(a1-a2)^2*R*Z^2 + 2*b1^2*(a1+a2)
# b1*a2*(a1-a2)*Z^3 + (a1-a2)^2*R*Z^2 + b1^2*(a1+a2)

### Solver:
solve.ht2a = function(b, a, R) {
	x.sum = roots(c(b[1]*a[2]*(a[1]-a[2]), (a[1]-a[2])^2*R, 0, b[1]^2*(a[1]+a[2])))
	x.sum = x.sum[x.sum != 0]
	if(length(x.sum) == 0) { return(list(sol=NA, p=NA)); }
	xy = -b[1]/(a[1]-a[2]) / x.sum
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example 1:
b = 3
a = c(1,2)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
a[1]*x^3*y + a[2]*x*y^3 + b[1]*x
a[1]*x*y^3 + a[2]*x^3*y + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))


### Example 2:
b = 3
a = c(1,2)
R = 6
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2]
sol

### Test
a[1]*x^3*y + a[2]*x*y^3 + b[1]*x
a[1]*x*y^3 + a[2]*x^3*y + b[1]*y

### Classic Polynomial
round0.p(poly.calc(sol$sol[,1]))


###
a = c(1,2)
p = sapply(-6:6, function(r) print(solve.ht2a(6*r, a, 6)$p))
p = sapply(-6:6, function(r) print(solve.ht2a(6*r, c(1,1/2), 9)$p)) # interesting
p = sapply(-6:6, function(r) print(solve.ht2a(2*r, c(-1/3,2/3), 12)$p)) # interesting
p = sapply(-6:6, function(r) print(solve.ht2a(r, c(-1/3,2/3), 12)$p)) # interesting
# type: x - x^3 + x^5


a1 = a[1]; a2 = a[2]; b1 = b[1]; R = R[1];
#
(a1*x^4 + a2*x^4 + b1*x - R) *
	((- a2^3*b1^3) +
	(- 3*R^2*a1*a2^3 + 3*R^2*a1^2*a2^2 - R^2*a1^3*a2 + R^2*a2^4)*x +
	(- 3*R*a1^2*a2^2*b1 + 2*R*a1^3*a2*b1 + R*a2^4*b1)*x^2 +
	(a1*a2^3*b1^2 + 2*a1^2*a2^2*b1^2 - a1^3*a2*b1^2 - 2*a2^4*b1^2)*x^3 +
	(- 2*R*a1*a2^4 + 2*R*a1^3*a2^2 - R*a1^4*a2 + R*a2^5)*x^5 +
	(a1*a2^4*b1 + a1^2*a2^3*b1 - a1^3*a2^2*b1 - a2^5*b1)*x^6 )


### full: un-decomposed
(R*a2^3*b1^3) +
(3*R^3*a1*a2^3 - 3*R^3*a1^2*a2^2 + R^3*a1^3*a2 - R^3*a2^4 - a2^3*b1^4)*x +
(- 3*R^2*a1*a2^3*b1 + 6*R^2*a1^2*a2^2*b1 - 3*R^2*a1^3*a2*b1)*x^2 +
(- R*a1*a2^3*b1^2 - 5*R*a1^2*a2^2*b1^2 + 3*R*a1^3*a2*b1^2 + 3*R*a2^4*b1^2)*x^3 +
(2*a1^2*a2^2*b1^3 - a1^3*a2*b1^3 - 3*a2^4*b1^3)*x^4 +
(- 2*R*a1*a2^4*b1 - 4*R*a1^2*a2^3*b1 + 2*R*a1^3*a2^2*b1 + R*a1^4*a2*b1 + 3*R*a2^5*b1)*x^6 +
(4*a1^2*a2^3*b1^2 - a1^4*a2*b1^2 - 3*a2^5*b1^2)*x^7 +
(- R*a1*a2^5 - 2*R*a1^2*a2^4 + 2*R*a1^3*a2^3 + R*a1^4*a2^2 - R*a1^5*a2 + R*a2^6)*x^9 +
(2*a1^2*a2^4*b1 - a1^4*a2^2*b1 - a2^6*b1)*x^10


###################################
###################################

