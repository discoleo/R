########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Base Types
###
### draft v.0.1a-sol


### Asymmetric Polynomial Systems: 2 Variables
### Base Types

### Example:
# x^n + b*y = R1
# y^n + b*x = R2
# R1 != R2

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

###############
### Order 2 ###

# x^2 + b*y = R1
# y^2 + b*x = R2

### Sum =>
x^2 + y^2 + b*(x+y) - R1 - R2 # = 0
S^2 - 2*x*y + b*S - R1 - R2
# 2*x*y = S^2 + b*S - R1 - R2

### Sum(x*...) =>
x^3 + y^3 + 2*b*x*y - R1*x - R2*y # = 0
x^3 + y^3 + b*(S^2 + b*S - R1 - R2) - R1*x - R2*y # = 0
S^3 - 3*x*y*S + b*(S^2 + b*S - R1 - R2) - R1*x - R2*y # = 0
2*S^3 - 3*(S^2 + b*S - R1 - R2)*S + 2*b*(S^2 + b*S - R1 - R2) - 2*R1*x - 2*R2*y # = 0
-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2) - 2*R1*x - 2*R2*y # = 0
# 2*R1*x + 2*R2*y = -S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)

### Sum(y*...) =>
x^2*y + y^2*x + b*(x^2+y^2) - R1*y - R2*x # = 0
x*y*(x+y) + b*(S^2 - 2*x*y) - R1*y - R2*x # = 0
x*y*S + b*(S^2 - (S^2 + b*S - R1 - R2)) - R1*y - R2*x # = 0
x*y*S - b*(b*S - R1 - R2) - R1*y - R2*x # = 0
(S^2 + b*S - R1 - R2)*S - 2*b*(b*S - R1 - R2) - 2*R1*y - 2*R2*x # = 0
S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2) - 2*R1*y - 2*R2*x # = 0
# 2*R2*x + 2*R1*y = S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)

### Eq 2 / 3:
# 2*(R1^2 - R2^2)*x =
R1 * (-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2)) +
	- R2 * (S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2))
-(R1+R2)*S^3 - b*(R1+R2)*S^2 + (3*R1^2 + R2^2 + 4*R1*R2 + 2*b^2*R1+ 2*b^2*R2)*S - 2*b*(R1 + R2)^2
# 2*(R1^2 - R2^2)*y =
R1*(S^3 + b*S^2 - (R1 + R2 + 2*b^2)*S + 2*b*(R1 + R2)) +
	- R2*(-S^3 - b*S^2 + (3*R1 + 3*R2 + 2*b^2)*S - 2*b*(R1 + R2))
(R1+R2)*S^3 + b*(R1+R2)*S^2 - (R1^2 + 3*R2^2 + 4*R1*R2 + 2*b^2*R1 + 2*b^2*R2)*S + 2*b*(R1 + R2)^2

### Diff =>
(x-y)*(x+y) - b*(x-y) # = R1 - R2
(x-y)*S - b*(x-y) - R1 + R2 # = 0
(y - x)*(S - b) + R1 - R2 # = 0
2*(R1^2 - R2^2)*(y - x)*(S - b) + 2*(R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
# =>
((R1+R2)*S^3 + b*(R1+R2)*S^2 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S + 2*b*(R1 + R2)^2)*(S-b) +
	+ (R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
(R1+R2)*S^4 + b*(R1+R2)*S^3 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S^2 + 2*b*(R1 + R2)^2*S +
	- b*((R1+R2)*S^3 + b*(R1+R2)*S^2 - 2*(R1^2 + R2^2 + 2*R1*R2 + b^2*R1 + b^2*R2)*S + 2*b*(R1 + R2)^2) +
	+ (R1^3 - R1^2*R2 - R1*R2^2 + R2^3) # = 0
(R1+R2)*S^4 - (R1+R2)*(2*R1 + 2*R2 + 3*b^2)*S^2 +
	+ 2*b*(R1+R2)*(2*R1 + 2*R2 + b^2)*S +
	+ (R1+R2)*(R1^2+R2^2 - 2*R1*R2) - 2*b^2*(R1+R2)^2 # = 0
S^4 - (2*R1 + 2*R2 + 3*b^2)*S^2 + 2*b*(2*R1 + 2*R2 + b^2)*S +
	+ (R1^2+R2^2 - 2*R1*R2) - 2*b^2*(R1+R2) # = 0

### Diff(x*...) =>
x^3 - y^3 - R1*x + R2*y # = 0
(x-y)*(S^2 - x*y) - R1*x + R2*y # = 0
(S-b)*(x-y)*(S^2 - x*y) - (R1*x - R2*y)*(S - b) # = 0
(R2 - R1)*(S^2 - b*S + R1 + R2) + 2*(R1*x - R2*y)*(S - b) # = 0


### TODO:
# - solve;
# - evaluate higher powers;


### Solution:
solve.asym.S2P2 = function(R, b) {
	R12 = R[1] + R[2]
	coeff = c(1, 0, - (2*R[1] + 2*R[2] + 3*b[1]^2), 2*b[1]*(2*R[1] + 2*R[2] + b[1]^2),
		(R[1]^2+R[2]^2 - 2*R[1]*R[2]) - 2*b[1]^2*(R[1]+R[2]))
	S = roots(coeff)
	len = length(S)
	#
	div = 2*(R[1] - R[2])*R12
	x = -R12*S^3 - b*R12*S^2 + (2*R12^2 + R12*(R[1] - R[2]) + 2*b^2*R12)*S - 2*b*R12^2
	x = x / div;
	y = R12*S^3 + b*R12*S^2 - (2*R12^2 - R12*(R[1] - R[2]) + 2*b^2*R12)*S + 2*b*R12^2
	y = y / div;
	return(cbind(x=x, y=y))
}

### Examples:
R = c(1,2)
b = 1
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x


#########
### Ex 2:
# special Case:
# S^4 - 8*S^2 + 32;
R = c(1, -3)
b = 2
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 13 + 8*x - 2*x^2 + x^4
round0(err)


#########
### Ex 3:
# special Case: b0 == 0
# S^4 - 16*S^2 + 32*S;
R = c(3, -1)
b = 2
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 13 + 8*x - 6*x^2 + x^4
round0(err)


#########
### Ex 4:
# special Case: b0 == 0
# S^4 - 31*S^2 + 78*S;
R = c(4, -2)
b = 3
sol = solve.asym.S2P2(R, b)
x = sol[,1]; y = sol[,2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x

### Classic Polynomial:
round0.p(poly.calc(x))
round0.p(poly.calc(x + y))

err = 34 + 27*x - 8*x^2 + x^4
round0(err)


###########

### Debug
R = c(2,3)
b = 3
x = 0.9563705005
y = 0.3617851553
S = x+y; R1 = R[1]; R2 = R[2];

### Test
x^2 + b[1]*y
y^2 + b[1]*x
