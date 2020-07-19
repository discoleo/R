
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: P2
### v.0.1-pre


####################
### Symmetric System

### n = 2
# (x + a)^2 + (y + a)^2 = R1
# x*y = R2

# x^2 + y^2 + 2*a*(x + y) + 2*a^2 - R1 = 0
# X = x + y =>
# X^2 + 2*a*X + 2*a^2 - R1 - 2*R2 = 0
# Step 1: Solve for X
# Step 2: x + y = X

### TODO:


###############

### n = 3
# (x + a)^3 + (y + a)^3 = R1
# x*y = R2

# x^3 + y^3 + 3*a*(x^2 + y^2) + 3*a^2*(x + y) + 2*a^3 - R1 = 0
# X = x + y =>
# X^3 - 3*R2*X + 3*a*X^2 - 6*a*R2 + 3*a^2*X + 2*a^3 - R1 = 0
# X^3 + 3*a*X^2 + 3*(a^2-R2)*X + 2*a^3 - R1 - 6*a*R2 = 0
# Step 1: Solve for X
# Step 2: x + y = X

# TODO: use exact solver for P3
library(polynom)

### Parameters
a = 2
R = c(1, 1)
### Solution

coeff = c(2*a^3 - R[1] - 6*a*R[2], 3*(a^2-R[2]), 3*a, 1)
X = solve(polynomial(coeff))
#
det = X^2 - 4*R[2]
x.m = sapply(det, function(det) if(Im(det) != 0 | Re(det) >= 0) sqrt(det) else complex(re=0, im=sqrt(-det)) )
x = (X + x.m)/2
y = (X - x.m)/2

### Test
(x + a)^3 + (y + a)^3
x*y

###############

### n = 4
# (x + a)^4 + (y + a)^4 = R1
# x*y = R2

### TODO:

