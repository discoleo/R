
########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### Heterogenous Symmetric
###
### draft v.0.1a

library(polynom)
library(pracma)


#############

### x^3 + b*y

# x^3 + b1*y = R
# y^3 + b1*x = R

# Diff =>
# x^3 - y^3 - b1*(x-y) = 0
# (x - y)*(x^2 + x*y + y^2 - b1) = 0
# => x = y *OR* x^2 + x*y + y^2 - b1 = 0;
# =>
# (x+y)^2 - x*y - b1 = 0
# x*y = (x+y)^2 - b1;

# Sum =>
# (x+y)^3 - 3*x*y*(x+y) + b1*(x+y) = 2*R
# s^3 - 3*(s^2 - b1)*s + b1*s - 2*R = 0
# s^3 - 2*b1*s + R = 0

### Step 2:
# Solve:
# x + y = s
# x*y = s^2 - b1

### Example
b = 3
R = 1
#
r.sum = roots(c(1,0, -2*b[1], R))
xy = r.sum^2 - b[1]
diff = sqrt(r.sum^2 - 4*xy + 0i)
x = (r.sum + diff)/2
y = (r.sum - diff)/2
sol = cbind(x, y)
sol # TODO: include also x = y cases

### Test
x^3 + b[1]*y
y^3 + b[1]*x
# Classic
err = x^6 - b[1]*x^4 - 2*R*x^3 + b[1]^2*x^2 + b[1]*R*x + R^2 - b[1]^3
round0(err)

### Classic
# b1*y = R - x^3
# =>
# (R - x^3)^3 / b1^3 + b1*x - R = 0
# (R - x^3)^3 + b1^4*x - R*b1^3
# (x^3 - R)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R - b1*x)^3 - b1^4*x + R*b1^3
# (x^3 + b1*x - R)^3 - 3*b1*x*(x^3 + b1*x - R)^2 + 3*b1^2*x^2*(x^3 + b1*x - R) - b1^3*x^3 - b1^4*x + R*b1^3
# let: p = (x^3 + b1*x - R)
# p^3 - 3*b1*x*p^2 + 3*b1^2*x^2*p - b1^3*p
# p*(p^2 - 3*b1*x*p + 3*b1^2*x^2 - b1^3)
# (x^3 + b1*x - R)*(x^6 - b1*x^4 - 2*R*x^3 + b1^2*x^2 + b1*R*x + R^2 - b1^3)



#############

### x^4 + b*y

# x^4 + b1*y = R
# y^4 + b1*x = R

# Diff =>
# x^4 - y^4 - b1*(x-y) = 0
# (x - y)*((x+y)*(x^2+y^2) - b1) = 0
# => x = y *OR* (x+y)*(x^2+y^2) - b1 = 0;
# =>
# (x+y)^2 - 2*x*y = b1 / s
# x*y = (s^2 - b1/s)/2

### Sum =>
# x^4 + y^4 + b1*s - 2*R = 0
# s^4 - 4*x*y*s^2 + 2*(x*y)^2 + b1*s - 2*R
# s^4 - 2*(s^4 - b1*s) + (s^2 - b1/s)^2/2 + b1*s - 2*R
# -s^4 + 2*b1*s + (s^2 - b1/s)^2/2 + b1*s - 2*R
# s^4 - 3*b1*s - (s^2 - b1/s)^2/2 + 2*R
# 2*s^6 - 6*b1*s^3 - (s^3 - b1)^2 + 4*R*s^2
# s^6 - 4*b1*s^3 + 4*R*s^2 - b1^2


### Example
b = 2
R = 1
#
r.sum = roots(c(1,0,0, -4*b[1], 4*R, 0, - b[1]^2))
xy = (r.sum^2 - b[1]/r.sum)/2
r.diff = sqrt(r.sum^2 - 4*xy + 0i)
x = (r.sum + r.diff)/2
y = (r.sum - r.diff)/2
sol = cbind(x, y) # TODO: add also x = y cases;
sol

### Test
x^4 + b[1]*y
y^4 + b[1]*x

