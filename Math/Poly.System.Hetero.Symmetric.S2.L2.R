########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S2
### Heterogeneous Symmetric
### Multiple/Simple Leading Terms
###
### draft v.0.1a

### Polynomial Systems
### Heterogeneous Symmetric
### Leading Terms: Multiple Simple


### Types:

### 2 Leading Terms:
# B3.1) a1*x^3 + a2*y^3 + b*x = R; (P3 => P6)
# B3.2) a1*x^3 + a2*y^3 + b*x*y = R; (P3 => P6)
# B3.3) a1*x^3 + a2*y^3 + b2*x*y + b1*x = R; (P3 => P6)


###############

###############
### History ###
###############


### v.0.1a:
# - moved sections with Multiple Leading terms:
#   from Poly.System.Hetero.Symmetric.R
#   to this file;


####################
####################

### Helper Functions

source("Polynomials.Helper.R")


##########################
##########################

########################
### Leading Terms: 2 ###
########################

### Leading: 2 Simple
### a1*x^3 + a2*y^3 + b1*x

# a1*x^3 + a2*y^3 + b1*x = R
# a2*x^3 + a1*y^3 + b1*y = R

### Solution

### Diff =>
# x*y = S^2 + b1/(a1 - a2)

### Sum =>
2*(a1+a2)*S^3 + 3*b1*(a1+a2)/(a1 - a2)*S - b1*S + 2*R
### Eq:
(a1+a2)*S^3 + b1*(a1+2*a2)/(a1 - a2)*S + R # = 0


### Solver:

solve.S2Ht.LmP3 = function(b, a, R) {
	# TODO: Case a1 == a2;
	b2 = if(length(b) > 1) b[2] else 0; # Ext A1;
	div = a[1] - a[2];
	coeff = c(2*(a[1]+a[2]), 0, 3*b[1]*(a[1]+a[2])/div - b[1], 2*R[1])
	if(b2 != 0) coeff = coeff + c(0, 0, -2*b2/div, 0);
	x.sum = roots(coeff)
	xy = x.sum^2 + b[1]/div;
	x.diff = sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Examples:
# - has Fractions;

b = 1
a = c(1/2, 1/3)
R = 1
#
sol = solve.S2Ht.LmP3(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[1]*y

### Classic Polynomial:

R[1]^2*(a[1] - a[2])^3 + a[2]^2*b[1]^3 +
	- (a[1] - a[2])*b[1]*R[1] * (2*a[1]^2 - a[2]^2 - a[1]*a[2])*x +
	+ b[1]^2 * (a[1]^3 - a[2]^3)*x^2 +
	- 2*(a[1] - a[2])^3*(a[1] + a[2])*R[1]*x^3 +
	+ (a[1]^2 - a[2]^2)*b[1] * (2*a[1]^2 - a[2]^2 - a[1]*a[2])*x^4 +
	+ (a[1] - a[2])^3*(a[1] + a[2])^2*x^6


#########
### Ex 2:
b = -1
a = c(1/3, -2/3) # a2 = - 2*a1;
R = 1
#
sol = solve.S2Ht.LmP3(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[1]*y

### Classic Polynomial:
err = 5 + 3*x^2 + 6*x^3 + x^6
round0(err)


#########
### Ex 3:
b = 3
a = c(1/2, -1/4)
R = 1
#
sol = solve.S2Ht.LmP3(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[1]*y

err = 80 - 48*x + 48*x^2 - 8*x^3 + 12*x^4 + x^6
round0(err)


###############
### Extensions:
b = c(-1, -2)
a = c(1/3, -2/3)
R = 1
#
sol = solve.S2Ht.LmP3(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x + b[2]*(x+y)
a[2]*x^3 + a[1]*y^3 + b[1]*y + b[2]*(x+y)


#########################

### Leading Terms: > 1
### Variant:
### a1*x^3 + a2*y^3 + b1*x*y

# a1*x^3 + a2*y^3 + b1*x*y = R
# a2*x^3 + a1*y^3 + b1*x*y = R

### Solution

### Diff =>
# x*y = S^2

### Sum =>
(a1+a2)*S^3 - b1*S^2 + R


### Solver:
solve.htm = function(b, a, R) {
	x.sum = roots(c((a[1]+a[2]), - b[1], 0, R))
	xy = x.sum^2
	x.diff= sqrt(x.sum^2 - 4*xy + 0i)
	x = (x.sum + x.diff)/2
	y = (x.sum - x.diff)/2
	sol = cbind(x, y)
	sol = rbind(sol, sol[,2:1])
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}


### Example 1:
b = 3
a = c(1, 2)
R = 1
#
sol = solve.htm(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol


### Test
a[1]*x^3 + a[2]*y^3 + b[1]*x*y
a[2]*x^3 + a[1]*y^3 + b[1]*x*y

### Classic Polynomial
err = 1/9 + 1/3*x^2 - 2/3*x^3 + x^4 - x^5 + x^6
round0(err)


#########
### Ex 2:
b = 3
a = c(1/2, -1)
R = 1
#
sol = solve.htm(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Classic Polynomial
err = 4 + 12*x^2 + 4*x^3 + 36*x^4 + 6*x^5 + x^6
round0(err)


#########
### Ex 3:
b = -1 # variants: - 1/2; 1/2;
a = c(1, -1/2)
R = 1/2
#
sol = solve.htm(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Classic Polynomial
err = 1 - 2*x^2 - 2*x^3 + 4*x^4 + 2*x^5 + x^6
round0(err)


err = R^2 + R*b[1]*x^2 - 2*R*(a[1] + a[2])*x^3 + b[1]^2*x^4 - b[1]*(a[1] + a[2])*x^5 +
	+ (a[1] + a[2])^2*x^6
round0(err)



###################################

### Variant:
### a1*x^3 + a2*y^3 + b2*x*y + b1*x

# a1*x^3 + a2*y^3 + b2*x*y + b1*x = R
# a2*x^3 + a1*y^3 + b2*x*y + b1*y = R

### Solution

### Diff =>
# x*y = S^2 + b1/(a1 - a2)

### Sum =>
2*(a1-a2)*(a1+a2)*S^3 - 2*b2*(a1-a2)*S^2 + 3*b1*(a1+a2)*S - b1*(a1-a2)*S - 2*b1*b2 + 2*(a1-a2)*R


solve.ht2a = function(b, a, R) {
	a.s = a[1]+a[2]; a.d = a[1]-a[2];
	coeffs = c(2*a.d*a.s, -2*b[2]*a.d, 3*b[1]*a.s - b[1]*a.d, -2*b[1]*b[2] + 2*a.d*R)
	r.sum = roots(coeffs)
	xy = r.sum^2 + b[1]/(a[1] - a[2])
	r.diff = sqrt(r.sum^2 - 4*xy + 0i)
	x = (r.sum + r.diff)/2
	y = (r.sum - r.diff)/2
	sol = cbind(x, y) # TODO: add also x = y cases;
	sol = round0(rbind(sol, sol[,2:1]))
	p = round0.p(poly.calc(sol[,1]))
	return(list(sol=sol, p=p))
}

### Example:
a = c(2/3, 1/3)
b = c(3, -1)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[2]*x*y + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[2]*x*y + b[1]*y


### Example 2:
a = c(2/3, -1/3)
b = c(1, -3)
R = 1
#
sol = solve.ht2a(b, a, R)
x = sol$sol[,1]; y = sol$sol[,2];
sol

### Test
a[1]*x^3 + a[2]*y^3 + b[2]*x*y + b[1]*x
a[2]*x^3 + a[1]*y^3 + b[2]*x*y + b[1]*y


### Classic Polynomial:
a1 = a[1]; a2 = a[2]; b1 = b[1]; b2 = b[2];

(- 3*R^2*a1*a2^2 + 3*R^2*a1^2*a2 - R^2*a1^3 + R^2*a2^3 - a2^2*b1^3) +
	b[1]*(a2 - a1)*(R*a1*a2 - 2*R*a1^2 + R*a2^2 - 2*a2*b1*b2)*x +
	(a2 - a1)*(- 2*R*a1*a2*b2 + R*a1^2*b2 + R*a2^2*b2 + a1*a2*b1^2 - a2*b1*b2^2 +
		+ a1*b1*b2^2 + a1^2*b1^2 + a2^2*b1^2)*x^2 +
	(a1 - a2)^2 * (2*R*a1^2 - 2*R*a2^2 + a1*b1*b2 + 2*a2*b1*b2)*x^3 +
	(a1 - a2)^2 * (- 3*a1*a2*b1 - a1*b2^2 - 2*a1^2*b1 + a2*b2^2 - a2^2*b1)*x^4 +
	b[2]*(a[1] - a[2])^3*(a[1] + a[2])*x^5 +
	(a[2] - a[1])^3*(a[1] + a[2])^2*x^6

