########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Base Types
###
### draft v.0.1a


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

### TODO:
# - solve;
# - evaluate higher powers;


### Debug
R = c(2,3)
b = 3
x = 0.9563705005
y = 0.3617851553
S = x+y; R1 = R[1]; R2 = R[2];

### Test

x^2 + b[1]*y
y^2 + b[1]*x
